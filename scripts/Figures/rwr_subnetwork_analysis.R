# ============================================================
# Script : reseau HPE (seeds/signal) + recuperation de variants
#          supplementaires (voie signal / voie diagnostique)
# ------------------------------------------------------------
# Objectif final :
#   - full_cohort_network : reseau (edges) signal+seed, avec le
#     statut d'evidence variant par gene (node1)
#   - node_recap           : table par gene, pensee pour lecture
#     directe par le biologiste (statut, preuves, patients
#     concernes) -> export article/tables/
#
# Trois niveaux d'evidence variant, du plus strict au plus
# recupere :
#   "seed"            -> variant strict initial (filtrage
#                         diagnostique standard, cf. patient_data.R)
#   "variant_refined" -> variant recupere avec PRED_value_reworked
#                         >= 0.7 (pLI ignore), sur un gene SIGNAL
#                         de la cohorte HPE (logique reseau)
#   "variant_diag"    -> variant recupere avec filtres qualite
#                         assouplis (seq_depth >= 10, INHOUSE_FREQ
#                         <= 0.001, DB_freq <= 0.01,
#                         PRED_value_reworked >= 0.3) ET gene
#                         diagnostique connu pour CE patient
#                         (clinic_clean) -> logique clinique
# ============================================================

library(Matrix)
library(tidyverse)
library(parallel)
library(tibble)
library(igraph)

setwd(paste(here::here(), "projects", sep = "/"))


source("article/scripts/custom_fct.R")
source("article/scripts/RECENT_patient_gene_pvalue.R") # process_group() / process_patient()

dir.create("article/networks", recursive = TRUE, showWarnings = FALSE)
dir.create("article/tables", recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 0. Chargement des objets de base
# ============================================================

load("article/transit_data/data_3p.RData")
load("article/transit_data/RWR/HPE_RWR_results_.RData")


load("article/transit_data/patient_variants.RData")
HPE_patient_variants <- patient_variants
rm(patient_variants)

load("article/transit_data/permutation_bank_cds_weighted.RData")
bank <- bank_blend

coex_genes <- rownames(data_3p$coex_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

patient_ids <- read.csv("article/transit_data/patient_ids.csv", stringsAsFactors = FALSE)[[1]]
names(HPE_patient_variants) <- patient_ids
colnames(HPE_results) <- patient_ids # appariement explicite patient <-> colonne RWR

# ============================================================
# Numerotation anonyme des patients (1, 2, 3, ...), utilisee a la
# place de l'identifiant patient dans les colonnes de recap
# destinees au biologiste (variant_diag_patients / variant_refined_patients).
# L'ordre suit celui de patient_ids.csv (stable, reproductible).
# Cette table de correspondance est sauvegardee a part, dans
# transit_data (PAS dans article/tables/), pour ne pas exposer les
# identifiants patients dans les tables destinees a l'article.
# ============================================================
patient_number_map <- setNames(seq_along(patient_ids), patient_ids)
save(patient_number_map, file = "article/transit_data/patient_number_map.RData")

HPE_seeds <- lapply(HPE_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(HPE_seeds) <- patient_ids

all_seed_genes <- unique(unlist(HPE_seeds))

seed_freq_cohort <- tibble(gene = unlist(HPE_seeds)) %>%
    dplyr::count(gene, name = "seed_n")

# clinic_clean : genes diagnostiques par patient (colonne "patient" = id, "genes" = symbols "/"-separes)
clinic_clean <- read.csv("/home/jules/Documents/phd/Data/vcf_HPE/clinic_clean.csv", header = TRUE, stringsAsFactors = FALSE)

get_diag_genes <- function(pid) {
    row <- filter(clinic_clean, patient == pid)
    if (nrow(row) == 0 || is.na(row$genes[1]) || row$genes[1] == "") {
        return(character(0))
    }
    strsplit(row$genes[1], "/")[[1]] %>% trimws()
}
diag_genes_by_patient <- setNames(lapply(patient_ids, get_diag_genes), patient_ids)

HPE_seeds %>%
    unlist() %>%
    unique() %>%
    sort()

# ============================================================
# 1. Genes signal, cohorte HPE complete (RWR + permutation bank)
# ============================================================

cohort_result <- process_profiles(
    patient_ids      = names(HPE_seeds),
    RWR_results      = HPE_results,
    seeds            = HPE_seeds,
    permutation_bank = bank
)
signal_genes <- cohort_result$gene[cohort_result$signal_final]


# ============================================================
# 2. Score d'effet variant (identique a patient_data.R)
# ============================================================

score_variant_effect <- function(PRED_value, PRED_effect, CLNSIG, PRED_splice) {
    high_impact <- c("frameshift deletion", "frameshift insertion", "stopgain")
    moderate_high <- c("stoploss", "startloss")
    moderate <- c("nonframeshift deletion", "nonframeshift insertion", "nonframeshift substitution")

    effect_score <- dplyr::case_when(
        PRED_value == 1 & PRED_effect %in% high_impact ~ 0.9,
        PRED_value == 1 & PRED_effect %in% moderate_high ~ 0.8,
        PRED_value == 1 & PRED_effect %in% moderate ~ 0,
        PRED_value < 1 ~ pmin(pmax(PRED_value, 0), 1),
        TRUE ~ 0.3
    )

    splice_no <- c("0/1", "0/2", "0/3", "0/4")
    splice_low <- c("1/2", "1/3", "1/4", "2/4")
    splice_moderate <- c("2/2", "2/3", "3/4")
    splice_high <- c("3/3", "4/4")

    effect_score <- pmax(effect_score, dplyr::case_when(
        PRED_effect %in% splice_high ~ 0.9,
        PRED_effect %in% splice_moderate ~ 0.8,
        PRED_effect %in% splice_low ~ 0.5,
        PRED_effect %in% splice_no ~ 0,
        TRUE ~ 0.3
    ))

    is_clinvar_patho <- CLNSIG %in% c("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic")
    pmax(effect_score, ifelse(is_clinvar_patho, 0.95, 0))
}


# ============================================================
# 2bis. Known HPE genes (Orphanet, Holoprosencephaly gene list,
#       accessed 2026-10-15)
# ------------------------------------------------------------
#  CRIPTO is the current HGNC symbol of TDGF1: both are kept so
#  the match works whatever the symbol version of the network.
# ============================================================

HPE_established_genes <- c(
    "SHH", "ZIC2", "SIX3", "TGIF1", "CDON", "DISP1", "DLL1", "FGF8", "FGFR1", "STIL",
    "KMT2D", "PPP1R12A", "RAD21", "SMC1A", "SMC3", "STAG2", "CNOT1", "CENPF", "DHCR7",
    "GLI2", "SUFU", "GAS1", "TDGF1", "FOXH1", "NODAL", "BOC", "PTCH1"
)

hpe_gene_tier <- function(gene) {
    ifelse(gene %in% HPE_established_genes, "HPE_genes", NA_character_)
}
# ============================================================
# 3. Recuperation de variants supplementaires par patient
#    (voie A "signal" / voie B "diag")
# ============================================================

patient_files_all <- list.files(path = "/home/jules/Documents/phd/Data/vcf_HPE/RESULTS/", full.names = TRUE)
names(patient_files_all) <- patient_ids # meme ordre que patient_ids (cf. patient_data.R)

variant_key_cols <- c("Chr", "Start", "End", "Ref", "Alt", "Gene.refGeneWithVer")

process_patient_relaxed <- function(pid) {
    patient_data <- read.table(patient_files_all[[pid]], header = TRUE, sep = "\t", fill = TRUE, quote = "")
    patient_data$seq_depth <- patient_data[, 17] %>% as.numeric()
    patient_data$pLI <- cleaning_pLI(patient_data$pLI)
    patient_data <- patient_data %>%
        mutate(Gene.refGeneWithVer = strsplit(as.character(Gene.refGeneWithVer), ";")) %>%
        unnest(Gene.refGeneWithVer) %>%
        as.data.frame()

    patient_data <- filter(
        patient_data,
        nchar(Ref) < 30,
        nchar(Alt) < 30,
        seq_depth >= 10,
        INHOUSE_FREQ <= 0.001,
        DB_freq <= 0.01
    )

    # PRED_value_reworked calcule une seule fois, reutilise par les 2 voies
    patient_data$PRED_value_reworked <- seq_len(nrow(patient_data)) %>% sapply(function(i) {
        score_variant_effect(
            patient_data[i, "PRED_value"], patient_data[i, "PRED_effect"],
            patient_data[i, "CLNSIG"], patient_data[i, "PRED_splice"]
        )
    })

    strict_df <- HPE_patient_variants[[pid]]
    diag_genes <- diag_genes_by_patient[[pid]]

    # ---- Voie A : "variant_refined" ------------------------------------
    # pLI ignore, PRED_value_reworked >= 0.7, restreint aux genes SIGNAL
    dataA <- filter(
        patient_data,
        DB_freq <= 0.001,
        seq_depth > 20,
        INHOUSE_FREQ <= 0.0001,
        nchar(Ref) < 30,
        nchar(Alt) < 30,
        PRED_value_reworked >= 0.7,
        Gene.refGeneWithVer %in% signal_genes
    )
    new_refined <- if (nrow(dataA) > 0) dplyr::anti_join(dataA, strict_df, by = variant_key_cols) else dataA

    # ---- Voie B : "variant_diag" ----------------------------------------
    # filtres qualite assouplis + gene diagnostique connu pour ce patient
    dataB <- filter(
        patient_data,
        nchar(Ref) < 30,
        nchar(Alt) < 30,
        seq_depth >= 10,
        INHOUSE_FREQ <= 0.001,
        DB_freq <= 0.01,
        PRED_value_reworked >= 0.3,
        Gene.refGeneWithVer %in% diag_genes
    )
    new_diag <- if (nrow(dataB) > 0) dplyr::anti_join(dataB, strict_df, by = variant_key_cols) else dataB

    list(new_refined = new_refined, new_diag = new_diag)
}

recovered <- mclapply(patient_ids, process_patient_relaxed, mc.cores = 15)
names(recovered) <- patient_ids

new_refined_by_patient <- lapply(recovered, `[[`, "new_refined")
new_diag_by_patient <- lapply(recovered, `[[`, "new_diag")


# ============================================================
# 4. patient_variants_new : strict + voie A + voie B, avec
#    colonnes seed / variant_type
# ------------------------------------------------------------
#  Chevauchement voie A / voie B pour un meme variant (gene a la
#  fois signal ET diagnostique) : priorite a "variant_diag"
#  (justification clinique directe), le variant n'est compte
#  qu'une fois.
# ============================================================

build_refined_patient_df <- function(pid) {
    strict_df <- HPE_patient_variants[[pid]] %>%
        mutate(seed = TRUE, variant_type = "seed")

    refined_df <- new_refined_by_patient[[pid]]
    diag_df <- new_diag_by_patient[[pid]]

    if (nrow(refined_df) > 0 && nrow(diag_df) > 0) {
        refined_df <- dplyr::anti_join(refined_df, diag_df, by = variant_key_cols)
    }
    if (nrow(refined_df) > 0) {
        refined_df <- refined_df %>%
            mutate(seed = Gene.refGeneWithVer %in% HPE_seeds[[pid]], variant_type = "variant_refined")
    }
    if (nrow(diag_df) > 0) {
        diag_df <- diag_df %>%
            mutate(seed = Gene.refGeneWithVer %in% HPE_seeds[[pid]], variant_type = "variant_diag")
    }

    bind_rows(strict_df, refined_df, diag_df)
}

patient_variants_new <- setNames(lapply(patient_ids, build_refined_patient_df), patient_ids)
save(patient_variants_new, file = "article/transit_data/patient_variants_new.RData")


# ============================================================
# 5. Reseau complet (degre) et filtrage des seeds : connectees
#    a au moins un gene signal (voisins dans le reseau complet)
# ============================================================

g_full <- graph_from_adjacency_matrix(
    data_3p$coex_case,
    mode = "undirected", weighted = TRUE, diag = FALSE
)
V(g_full)$name <- V(g_full)$name %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

full_degree <- data.frame(gene = V(g_full)$name, degree = igraph::degree(g_full))

filter_seeds_connected_to_signal <- function(seed_genes, signal_genes, full_graph) {
    anchor_genes <- intersect(unique(signal_genes), V(full_graph)$name)
    if (length(anchor_genes) == 0) {
        return(character(0))
    }
    neighbor_genes <- unique(unlist(lapply(anchor_genes, function(g) neighbors(full_graph, g)$name)))
    intersect(seed_genes, neighbor_genes)
}

connected_seed_genes <- filter_seeds_connected_to_signal(all_seed_genes, signal_genes, g_full)


# ============================================================
# 6. Construction du reseau (dataframe node1/node2)
# ------------------------------------------------------------
# Sous-graphe induit par l'union des genes signal et des seeds
# directement connectees a un gene signal (cf. section 5).
#
# Le reseau est duplique (edges + edges_mirror) pour que CHAQUE
# gene apparaisse comme node1 avec l'integralite de ses aretes
# (annotations completes cote node1) -> a filtrer duplicate==TRUE
# une fois l'import Cytoscape fait.
#
# Colonne "link" -- PROPRIETE DE L'ARETE (identique sur la ligne
# et sa version miroir) :
#   "link"      : cor >= 0.4, quel que soit le reste (priorite absolue)
#   "link_weak" : cor < 0.4 MAIS egale a la correlation maximale
#                 d'AU MOINS UN des deux noeuds de l'arete (le
#                 noeud considere son extremite comme sa meilleure
#                 connexion, meme si l'autre extremite a mieux
#                 ailleurs). Si un noeud n'a qu'une seule arete,
#                 celle-ci est forcement son max -> "link_weak"
#                 (ou "link" si >= 0.4), donc duplicate et
#                 non-duplicate portent la meme valeur.
#   "none"      : ni l'un ni l'autre, pour les deux noeuds.
# ============================================================

build_full_cohort_network <- function(signal_genes, seed_genes, full_graph,
                                      full_degree_df, seed_freq_cohort_df, all_seed_genes) {
    group_genes <- unique(c(signal_genes, seed_genes))
    group_genes <- intersect(group_genes, V(full_graph)$name)
    if (length(group_genes) < 2) {
        message("Moins de 2 genes dans le reseau, impossible de construire le graphe")
        return(NULL)
    }

    subg <- induced_subgraph(full_graph, vids = group_genes)

    gene_attrs <- data.frame(gene = group_genes) %>%
        mutate(signal = gene %in% signal_genes, seed = gene %in% all_seed_genes) %>%
        left_join(full_degree_df, by = "gene") %>%
        left_join(seed_freq_cohort_df, by = "gene") %>%
        mutate(
            seed_n = replace_na(seed_n, 0),
            condition = case_when(
                seed & signal ~ "seed+signal",
                seed ~ "seed",
                signal ~ "signal",
                TRUE ~ NA_character_
            )
        )

    # Une ligne par arete NON orientee (pas encore de duplication)
    edges_raw <- igraph::as_data_frame(subg, what = "edges") %>%
        dplyr::select(from, to, cor = weight)
    if (nrow(edges_raw) == 0) {
        message("Aucune arete dans le sous-reseau")
        return(NULL)
    }

    # Correlation max par noeud, tous voisins confondus (le noeud
    # peut apparaitre comme "from" ou "to" selon les aretes)
    node_max_cor <- bind_rows(
        edges_raw %>% dplyr::select(node = from, cor),
        edges_raw %>% dplyr::select(node = to, cor)
    ) %>%
        group_by(node) %>%
        summarise(max_cor = max(cor), .groups = "drop")

    # Statut "link" calcule UNE FOIS par arete non orientee, en
    # regardant le max de CHAQUE extremite -> propriete de l'arete
    edges_raw <- edges_raw %>%
        left_join(node_max_cor, by = c("from" = "node")) %>%
        dplyr::rename(max_cor_from = max_cor) %>%
        left_join(node_max_cor, by = c("to" = "node")) %>%
        dplyr::rename(max_cor_to = max_cor) %>%
        mutate(
            link = case_when(
                cor >= 0.4 ~ "link",
                cor == max_cor_from | cor == max_cor_to ~ "link_weak",
                TRUE ~ "none"
            )
        ) %>%
        dplyr::select(from, to, cor, link)

    edges <- edges_raw %>%
        dplyr::select(node1 = from, node2 = to, cor, link) %>%
        left_join(gene_attrs, by = c("node1" = "gene")) %>%
        mutate(duplicate = FALSE)

    edges_mirror <- edges_raw %>%
        dplyr::select(node1 = to, node2 = from, cor, link) %>%
        left_join(gene_attrs, by = c("node1" = "gene")) %>%
        mutate(duplicate = TRUE)

    bind_rows(edges, edges_mirror) %>%
        dplyr::select(node1, node2, signal, seed, seed_n, condition, degree, cor, link, duplicate)
}

full_cohort_network <- build_full_cohort_network(
    signal_genes = signal_genes,
    seed_genes = connected_seed_genes,
    full_graph = g_full,
    full_degree_df = full_degree,
    seed_freq_cohort_df = seed_freq_cohort,
    all_seed_genes = all_seed_genes
)

nrow(full_cohort_network)

# ============================================================
# 7. Evidence variant par gene (a partir de patient_variants_new)
# ------------------------------------------------------------
#  Pour chaque gene du reseau, priorite d'affichage :
#     variant (seed strict) > variant_diag > variant_refined > none
#  + details utiles au biologiste : nb et liste des patients
#    apportant une preuve supplementaire (refined / diag)
#
#  Les patients sont identifies par un NUMERO (1, 2, 3, ...),
#  via patient_number_map (cf. section 0), plutot que par leur
#  identifiant d'origine, dans les colonnes de recap destinees
#  au biologiste (variant_diag_patients / variant_refined_patients).
# ============================================================

all_new_variants <- purrr::imap_dfr(patient_variants_new, function(df, pid) {
    df %>% mutate(patient = pid)
}) %>%
    filter(variant_type != "seed") %>% # uniquement les preuves supplementaires
    mutate(patient = patient_number_map[patient]) # identifiant -> numero anonyme

variant_evidence_by_gene <- all_new_variants %>%
    group_by(gene = Gene.refGeneWithVer, variant_type) %>%
    summarise(patients = paste(sort(unique(patient)), collapse = ";"), n_patients = n_distinct(patient), .groups = "drop") %>%
    pivot_wider(
        names_from  = variant_type,
        values_from = c(patients, n_patients),
        names_glue  = "{variant_type}_{.value}"
    )

# noms de colonnes attendus si absents (aucun gene concerne par une voie)
for (col in c(
    "variant_refined_n_patients", "variant_refined_patients",
    "variant_diag_n_patients", "variant_diag_patients"
)) {
    if (!col %in% names(variant_evidence_by_gene)) {
        variant_evidence_by_gene[[col]] <- if (grepl("n_patients", col)) 0L else NA_character_
    }
}

variant_status_from_gene <- function(gene) {
    case_when(
        gene %in% all_seed_genes ~ "variant",
        gene %in% (all_new_variants %>% filter(variant_type == "variant_diag") %>% pull(Gene.refGeneWithVer)) ~ "variant_diag",
        gene %in% (all_new_variants %>% filter(variant_type == "variant_refined") %>% pull(Gene.refGeneWithVer)) ~ "variant_refined",
        TRUE ~ "none"
    )
}

full_cohort_network <- full_cohort_network %>%
    mutate(variant_status = variant_status_from_gene(node1))


# ============================================================
# 7bis. Annotations supplementaires par gene
# ------------------------------------------------------------
#  - label                    : alias de node1 pour le mapping
#                                Label Cytoscape, VIDE ("") des
#                                que variant_status == "none"
#                                (allege l'affichage des noeuds
#                                sans aucune preuve variant)
#  - HPE_genes                : tier HPE diagnostique connu
#                                (litterature, cf section 2bis) :
#                                "major" / "minor" / NA
#  - signal_variant_status     : combinaison condition (signal
#                                ou non) x variant_status, avec
#                                un libelle de niveau d'interet :
#                                  seed                  -> variant seed strict
#                                                          (repris par les filtres d'origine)
#                                  variant_diag_excluded -> variant diagnostique
#                                                          connu mais non retenu au
#                                                          filtrage strict (interet moyen,
#                                                          gene deja connu)
#                                  variant_refined       -> variant recupere post-RWR
#                                                          sur gene signal (interet faible)
#                                  no_variant             -> aucune preuve variant
#                                                          (interet tres faible)
#                                prefixe signal_ / non_signal_ selon condition
#  - isolated_low_confidence  : gene variant_refined/none, sans
#                                aucune arete forte ("link") ET
#                                sans arete (link ou link_weak)
#                                vers un gene seed+signal -> aucune
#                                ancre credible dans le reseau
# ------------------------------------------------------------
#  Note : node1 == g couvre deja TOUS les voisins de g (edges +
#  edges_mirror), pas besoin de refaire un induced_subgraph. Le
#  join sur node2 recupere la condition du voisin, absente
#  jusqu'ici du tableau (seul node1 etait annote).
# ============================================================

condition_by_gene <- full_cohort_network %>%
    dplyr::distinct(node1, condition) %>%
    dplyr::rename(node2_condition = condition)

node_flags <- full_cohort_network %>%
    dplyr::left_join(condition_by_gene, by = c("node2" = "node1")) %>%
    dplyr::group_by(node1) %>%
    dplyr::summarise(
        has_strong_link = any(link == "link"),
        linked_to_seed_signal = any(link %in% c("link", "link_weak") & node2_condition == "seed+signal"),
        .groups = "drop"
    )

signal_variant_status_from <- function(condition, variant_status) {
    signal_flag <- condition %in% c("signal", "seed+signal")
    variant_label <- dplyr::case_when(
        variant_status == "variant" ~ "seed",
        variant_status == "variant_diag" ~ "variant_diag_excluded",
        variant_status == "variant_refined" ~ "variant_refined",
        variant_status == "none" ~ "no_variant",
        TRUE ~ NA_character_
    )
    paste0(ifelse(signal_flag, "signal", "non_signal"), "_", variant_label)
}

# ============================================================
# 7ter. Fragments du reseau (identification automatique)
# ------------------------------------------------------------
#  Reproduit l'exploration Cytoscape, de facon automatique :
#   1. retrait des genes isolated_low_confidence ;
#   2. aretes "link" / "link_weak" uniquement (une seule ligne par
#      arete : duplicate == FALSE) ;
#   3. composantes connexes du graphe restant.
#
#  Etiquettes :
#   "main"     : plus grosse composante (contient normalement les
#                genes seed+signal deja explores)
#   "fragN"    : autres composantes de >= 2 genes contenant au moins
#                un gene avec preuve variant, numerotees par taille
#                decroissante (puis nb de genes variant, puis ordre
#                alphabetique du premier gene -> numerotation stable)
#   "isolated" : gene avec preuve variant reste seul (aucune arete
#                link/link_weak vers un gene conserve)
#
#  network_fragment ne concerne QUE les genes portant une preuve
#  variant (variant_status != "none"). Les genes sans preuve variant,
#  les genes du fragment principal et les genes
#  isolated_low_confidence recoivent NA (meme regle qu'avant).
# ============================================================

node_status <- full_cohort_network %>%
    dplyr::distinct(node1, condition, variant_status) %>%
    dplyr::left_join(node_flags, by = "node1") %>%
    dplyr::mutate(
        isolated_low_confidence = variant_status %in% c("variant_refined", "none") &
            !has_strong_link & !linked_to_seed_signal
    )

identify_network_fragments <- function(network_df, node_status, link_types = c("link", "link_weak")) {
    kept_nodes <- node_status %>%
        dplyr::filter(!isolated_low_confidence) %>%
        dplyr::pull(node1)

    kept_edges <- network_df %>%
        dplyr::filter(!duplicate, link %in% link_types, node1 %in% kept_nodes, node2 %in% kept_nodes) %>%
        dplyr::select(node1, node2)

    g_kept <- igraph::graph_from_data_frame(kept_edges, directed = FALSE, vertices = data.frame(name = kept_nodes))
    comp <- igraph::components(g_kept)

    comp_df <- tibble(gene = names(comp$membership), component = as.integer(comp$membership)) %>%
        dplyr::left_join(node_status %>% dplyr::select(gene = node1, condition, variant_status), by = "gene") %>%
        dplyr::group_by(component) %>%
        dplyr::mutate(
            component_size = dplyr::n(),
            n_variant_genes = sum(variant_status != "none"),
            n_seed_signal = sum(condition == "seed+signal")
        ) %>%
        dplyr::ungroup()

    # Composante principale = la plus grosse
    main_component <- comp_df %>%
        dplyr::distinct(component, component_size) %>%
        dplyr::arrange(dplyr::desc(component_size)) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::pull(component)

    # Fragments secondaires (>= 2 genes, au moins un gene variant)
    frag_ids <- comp_df %>%
        dplyr::filter(component != main_component, component_size >= 2, n_variant_genes >= 1) %>%
        dplyr::group_by(component, component_size, n_variant_genes) %>%
        dplyr::summarise(first_gene = min(gene), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(component_size), dplyr::desc(n_variant_genes), first_gene) %>%
        dplyr::mutate(fragment_label = paste0("frag", dplyr::row_number())) %>%
        dplyr::select(component, fragment_label)

    by_gene <- comp_df %>%
        dplyr::left_join(frag_ids, by = "component") %>%
        dplyr::mutate(
            fragment_label = dplyr::case_when(
                component == main_component ~ "main",
                !is.na(fragment_label) ~ fragment_label,
                component_size == 1 & variant_status != "none" ~ "isolated",
                TRUE ~ NA_character_ # fragment sans aucun gene variant
            ),
            network_fragment = dplyr::if_else(
                variant_status == "none" | fragment_label == "main",
                NA_character_, fragment_label
            )
        )

    summary <- by_gene %>%
        dplyr::filter(!is.na(fragment_label)) %>%
        dplyr::group_by(fragment = fragment_label) %>%
        dplyr::summarise(
            n_genes = dplyr::n(),
            n_variant_genes = sum(variant_status != "none"),
            n_seed_signal = sum(condition == "seed+signal"),
            variant_genes = paste(sort(gene[variant_status != "none"]), collapse = ";"),
            .groups = "drop"
        ) %>%
        dplyr::arrange(fragment != "main", fragment == "isolated", readr::parse_number(fragment))

    # Controle : les genes seed+signal devraient tous etre dans le fragment principal
    seed_signal_out <- by_gene %>%
        dplyr::filter(condition == "seed+signal", fragment_label != "main" | is.na(fragment_label)) %>%
        dplyr::pull(gene)
    if (length(seed_signal_out) > 0) {
        message("Genes seed+signal hors du fragment principal : ", paste(seed_signal_out, collapse = ", "))
    }

    list(by_gene = by_gene, summary = summary)
}

fragments <- identify_network_fragments(full_cohort_network, node_status)
print(fragments$summary, n = Inf)

dir.create("article/tables/Figure_4", recursive = TRUE, showWarnings = FALSE)
write.csv(fragments$summary, "article/tables/Figure_4/HPE_network_fragments_summary.csv", row.names = FALSE)

full_cohort_network <- full_cohort_network %>%
    dplyr::left_join(node_status %>% dplyr::select(node1, isolated_low_confidence), by = "node1") %>%
    dplyr::left_join(fragments$by_gene %>% dplyr::select(node1 = gene, network_fragment), by = "node1") %>%
    dplyr::mutate(
        HPE_genes = hpe_gene_tier(node1),
        signal_variant_status = signal_variant_status_from(condition, variant_status),
        label = ifelse(variant_status == "none", "", node1)
    ) %>%
    dplyr::select(
        node1, node2, label, signal, seed, seed_n, condition, degree, cor, link, duplicate,
        variant_status, signal_variant_status, HPE_genes, isolated_low_confidence, network_fragment
    )


# ============================================================
# 8. node_recap : table par gene lisible par le biologiste
# ------------------------------------------------------------
#  Une ligne par gene du reseau, avec :
#   - gene              : symbole
#   - condition         : statut HPE reseau (seed / signal / seed+signal)
#   - degree            : degre dans le reseau complet (hub-ness)
#   - seed_n            : nb de patients ou le gene est seed (strict)
#   - variant_status    : meilleur niveau d'evidence disponible
#                          (variant > variant_diag > variant_refined > none)
#   - HPE_genes            : tier HPE diagnostique connu
#                          (litterature, major/minor)
#   - signal_variant_status : combinaison condition x variant_status
#                          (cf. section 7bis)
#   - isolated_low_confidence : aucune ancre credible dans le reseau
#                          (cf. section 7bis)
#   - variant_diag_n_patients / variant_diag_patients       : preuve
#     clinique supplementaire (gene diagnostique connu, filtres assouplis).
#     variant_diag_patients contient des NUMEROS de patients
#     (cf. patient_number_map, section 0), pas les identifiants d'origine.
#   - variant_refined_n_patients / variant_refined_patients : preuve
#     reseau supplementaire (gene signal, PRED_value_reworked >= 0.7).
#     variant_refined_patients contient des NUMEROS de patients,
#     meme logique que ci-dessus.
# ============================================================

node_recap <- full_cohort_network %>%
    dplyr::select(
        node1, label, signal, seed, seed_n, condition, degree, variant_status,
        signal_variant_status, HPE_genes, isolated_low_confidence, network_fragment
    ) %>%
    dplyr::distinct(node1, .keep_all = TRUE) %>%
    dplyr::rename(gene = node1) %>%
    left_join(variant_evidence_by_gene, by = "gene") %>%
    mutate(
        variant_diag_n_patients    = replace_na(variant_diag_n_patients, 0L),
        variant_refined_n_patients = replace_na(variant_refined_n_patients, 0L)
    ) %>%
    dplyr::select(
        gene, label, condition, degree, seed_n, variant_status,
        signal_variant_status, HPE_genes, isolated_low_confidence, network_fragment,
        variant_diag_n_patients, variant_diag_patients,
        variant_refined_n_patients, variant_refined_patients
    ) %>%
    dplyr::arrange(desc(variant_status == "variant"), desc(seed_n), desc(degree))


# ============================================================
# 9. Sauvegarde
# ============================================================

write.csv(
    full_cohort_network,
    "article/tables/Figure_4/HPE_full_cohort_signal_seed_network_refined.csv",
    row.names = FALSE, quote = FALSE
)

write.csv(
    node_recap,
    "article/tables/Figure_4/HPE_full_cohort_signal_seed_node_recap_refined.csv",
    row.names = FALSE, quote = FALSE
)

library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "node_recap")

writeData(wb, "node_recap", node_recap, headerStyle = createStyle(textDecoration = "bold"))
freezePane(wb, "node_recap", firstRow = TRUE)
setColWidths(wb, "node_recap", cols = seq_len(ncol(node_recap)), widths = "auto")

saveWorkbook(
    wb,
    "article/tables/Figure_4/HPE_full_cohort_signal_seed_node_recap_refined.xlsx",
    overwrite = TRUE
)





hpe_seed_genes <- HPE_seeds %>%
    unlist() %>%
    unique()

hpe_check <- tibble(gene = HPE_established_genes) %>%
    mutate(
        in_network = gene %in% coex_genes,
        signal     = gene %in% signal_genes,
        seed       = gene %in% hpe_seed_genes
    ) %>%
    arrange(desc(in_network), desc(signal), desc(seed), gene)

print(hpe_check, n = Inf)

# Comptes
hpe_check %>%
    summarise(
        n_genes = n(),
        in_network = sum(in_network),
        signal = sum(signal),
        seed = sum(seed),
        seed_in_net = sum(seed & in_network),
        signal_and_seed = sum(signal & seed),
        signal_no_seed = sum(signal & !seed)
    )

# Listes utiles pour le texte
hpe_check %>%
    filter(!in_network) %>%
    pull(gene) # absents du réseau
hpe_check %>%
    filter(signal & seed) %>%
    pull(gene) # références Fig 4
hpe_check %>%
    filter(signal & !seed) %>%
    pull(gene) # signal par propagation seule
