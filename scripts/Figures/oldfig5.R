# Loading packages and functions
library(Matrix)
library(tidyverse)
library(parallel)
library(tibble)
library(igraph)
library(EGAD)
library(ggplot2)
library(cluster)
library(stats)
library(RANKS)
library(circlize)
library(ComplexHeatmap)
library(plotly)
library(WebGestaltR)
library(scales)

# Setting working directory
rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-4) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()

# Loading custom functions
source("article/scripts/custom_fct.R")
source("article/scripts/RECENT_patient_gene_pvalue.R") # process_group() / process_patient()
source("article/scripts/RECENT_rwr_normalization.R") # build_permutation_bank() / normalize_rwr_score()
source("article/scripts/RECENT_webgestalt_ORA.R")

load("article/transit_data/data_3p.RData")

load("article/transit_data/RWR_results_case.RData")
HPE_results <- as.matrix(RWR_results_z)
rm(RWR_results_z)

load("article/transit_data/patient_variants.RData")
HPE_patient_variants <- patient_variants
rm(patient_variants)

# Precomputed permutation bank (build_permutation_bank(), cf. RECENT_results_analyses.R)
# -> reused as-is, no RWR recomputation
load("article/transit_data/permutation_bank.RData") # -> `bank`

coex_genes <- rownames(data_3p$coex_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

context_genes <- rownames(data_3p$norm_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

HPE_seeds <- lapply(HPE_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(HPE_seeds) <- colnames(HPE_results)

# groupes fonctionnels (identiques a Figure_2.R)
BMP_TGFbeta <- c("FBN3", "FBN2", "NR3C1", "PCSK5", "PLOD1", "COL2A1", "LTBP1", "BMP2", "SMAD3", "TGIF1")
SHH <- c("BOC", "SHH", "GLI2", "LRP2", "CEP128", "CROCC", "GLI3", "SCUBE2", "ZIC2", "CREBBP", "DNHD1", "LAMB1", "PTCH1", "SIX3", "FOXA2")
NOTCH <- c("MKI67", "CIT", "ASPH", "LFNG", "RBL1", "SALL2", "ADAM19", "ASCL1", "NOTCH1", "NOTCH2", "NOTCH3")
FGF <- c("APLP2", "AXL", "ETV4", "FGFR4", "FGF8")
WNT_PCP <- c("MACF1", "CACHD1", "ITGA5", "LGR5", "FAT1", "TNFRSF19", "CELSR1", "CELSR2", "APC", "FAT3", "NCOA6", "APC2", "AXIN2", "LEF1", "PTK7", "TCF7L2")

pathway_groups <- list(
    BMP_TGFbeta = BMP_TGFbeta,
    SHH         = SHH,
    NOTCH       = NOTCH,
    FGF         = FGF,
    WNT_PCP     = WNT_PCP
)

get_patients_with_seed_in <- function(gene_set, seeds_list) {
    names(seeds_list)[sapply(seeds_list, function(seeds) any(seeds %in% gene_set))]
}

patient_groups <- lapply(pathway_groups, get_patients_with_seed_in, seeds_list = HPE_seeds)


# ============================================================
# 1. Genes signal par groupe (identique a Figure_2.R, step 2)
# ============================================================

pathway_group_results <- lapply(patient_groups, function(patient_ids) {
    if (length(patient_ids) == 0) {
        message("Groupe vide, ignore")
        return(NULL)
    }
    process_group(
        patient_ids      = patient_ids,
        RWR_results      = HPE_results,
        seeds            = HPE_seeds,
        permutation_bank = bank
    )
})

pathway_signal_genes <- lapply(pathway_group_results, function(df) {
    if (is.null(df)) {
        return(character(0))
    }
    df$gene[df$signal_final]
})

sapply(pathway_signal_genes, length)


# ============================================================
# 2. Extraction des seeds par groupe
# ------------------------------------------------------------
# Union des seeds (variants) des patients appartenant a chaque
# groupe de voie de signalisation.
# ============================================================

pathway_seed_genes <- lapply(patient_groups, function(patient_ids) {
    if (length(patient_ids) == 0) {
        return(character(0))
    }
    unique(unlist(HPE_seeds[patient_ids]))
})

sapply(pathway_seed_genes, length)


# ============================================================
# 3. Frequences des seeds
# ------------------------------------------------------------
#  - seed_n_cohort : nb de patients (cohorte HPE complete, tous
#                     groupes confondus) ou le gene est une seed
#  - seed_n_group   : nb de patients DU GROUPE ou le gene est
#                     une seed (calcule separement par groupe,
#                     a partir des memes patient_ids que ceux
#                     utilises pour pathway_seed_genes)
# ============================================================

count_seed_freq <- function(patient_ids, seeds_list, colname) {
    if (length(patient_ids) == 0) {
        return(data.frame(gene = character(0), n = integer(0)) %>% rename(!!colname := n))
    }
    seeds_list[patient_ids] %>%
        unlist() %>%
        table() %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        rename(gene = 1, n = Freq) %>%
        mutate(gene = as.character(gene)) %>%
        rename(!!colname := n)
}

# Frequence cohorte complete (tous les patients HPE, tous groupes confondus)
seed_freq_cohort <- count_seed_freq(names(HPE_seeds), HPE_seeds, "seed_n_cohort")

# Frequence par groupe
pathway_seed_freq_group <- lapply(names(patient_groups), function(pw) {
    count_seed_freq(patient_groups[[pw]], HPE_seeds, "seed_n_group")
})
names(pathway_seed_freq_group) <- names(patient_groups)


# ============================================================
# 4. Attributs genes (a l'echelle du reseau complet)
# ------------------------------------------------------------
#  - degree    : degre dans le reseau RWR complet (coex_case)
#  - HPE_genes : liste des genes diagnostiques HPE connus
# ============================================================

g_full <- graph_from_adjacency_matrix(
    data_3p$coex_case,
    mode     = "undirected",
    weighted = TRUE,
    diag     = FALSE
)
V(g_full)$name <- V(g_full)$name %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

full_degree <- data.frame(
    gene   = V(g_full)$name,
    degree = igraph::degree(g_full)
)

HPE_diag_genes <- c(
    "SHH", "GLI2", "FGF8", "ZIC2", "SIX3", "TGIF1", "PTCH1", "GAS1", "FGFR1", "BOC", "SCUBE2", "LRP2"
)


# ============================================================
# 4bis. Annotation HPO (Human Phenotype Ontology)
# ------------------------------------------------------------
# Meme source d'annotation que celle utilisee par run_ora_hpo()
# (RECENT_webgestalt_ORA.R) : enrichDatabase =
# "phenotype_Human_Phenotype_Ontology", enrichDatabaseType =
# "genesymbol", organism = "hsapiens".
#
# loadGeneSet() renvoie :
#  - $geneSet    : table AU FORMAT LONG (une ligne par gene) :
#                   colonnes "geneSet" (id HP:xxxxxxx), "description"
#                   (nom du terme), "gene". ATTENTION : malgre
#                   enrichDatabaseType = "genesymbol" (qui ne
#                   controle que le typage du gene list d'INTERET,
#                   pas celui de l'annotation source), la colonne
#                   "gene" de l'annotation HPO brute est en
#                   ENTREZID -> conversion ENTREZID -> SYMBOL
#                   obligatoire via gene_converter() avant tout
#                   matching avec les symboles utilises partout
#                   ailleurs dans le pipeline (coex_genes, seeds,
#                   signal genes, etc.)
#  - $geneSetDes : table geneSet (HP:xxxxxxx) -> description,
#                   une ligne par terme (dedupliquee). Utilisee ici
#                   pour resoudre un nom de terme fourni manuellement
#                   en identifiant HP:xxxxxxx.
# ============================================================

hpo_full <- WebGestaltR::loadGeneSet(
    organism = "hsapiens",
    enrichDatabase = "phenotype_Human_Phenotype_Ontology",
    enrichDatabaseType = "genesymbol"
)
hpo_geneset_des <- hpo_full$geneSetDes # geneSet (HP id), description (nom du terme, dedupliquee)

# Conversion ENTREZID -> SYMBOL de la colonne "gene" (annotation brute)
hpo_geneset <- hpo_full$geneSet %>%
    mutate(gene_symbol = gene_converter(as.character(gene), from = "ENTREZID", to = "SYMBOL")) %>%
    filter(!is.na(gene_symbol))

# Resout des noms de termes (ou des identifiants HP:xxxxxxx deja
# fournis tels quels) en identifiants HP:xxxxxxx, via geneSetDes
map_hpo_terms_to_ids <- function(term_names, geneset_des_df) {
    ids <- character(0)
    for (term in term_names) {
        term_trim <- trimws(term)
        if (grepl("^HP:[0-9]{7}$", term_trim)) {
            ids <- c(ids, term_trim)
            next
        }
        match_id <- geneset_des_df %>%
            dplyr::filter(tolower(trimws(description)) == tolower(term_trim)) %>%
            dplyr::pull(geneSet)
        if (length(match_id) == 0) {
            # repli : recherche partielle si aucune correspondance exacte
            match_id <- geneset_des_df %>%
                dplyr::filter(grepl(tolower(term_trim), tolower(description), fixed = TRUE)) %>%
                dplyr::pull(geneSet)
        }
        if (length(match_id) == 0) {
            message("Terme HPO non trouve (verifier l'orthographe) : ", term)
        }
        ids <- c(ids, match_id)
    }
    unique(ids)
}

get_genes_from_hpo_terms <- function(hpo_terms, geneset_df, geneset_des_df) {
    if (length(hpo_terms) == 0) {
        return(character(0))
    }
    hpo_ids <- map_hpo_terms_to_ids(hpo_terms, geneset_des_df)
    if (length(hpo_ids) == 0) {
        return(character(0))
    }
    geneset_df %>%
        dplyr::filter(geneSet %in% hpo_ids) %>%
        dplyr::pull(gene_symbol) %>%
        unique()
}

group_hpo_terms <- list(
    BMP_TGFbeta = c(
        "Pituitary dwarfism",
        "Ectopic anterior pituitary gland",
        "Anterior pituitary agenesis",
        "Median cleft lip and palate",
        "Hypoplastic iris stroma",
        "Basal encephalocele",
        "Ethmoidal encephalocele",
        "Abnormal hard palate morphology",
        "Midnasal stenosis",
        "Abnormality of the midnasal cavity",
        "Abnormality of the nasal cavity",
        "Median cleft palate"
    ),
    SHH = c(
        "Pituitary dwarfism",
        "Ectopic anterior pituitary gland",
        "Median cleft lip and palate",
        "Median cleft palate",
        "Median cleft lip",
        "Cleft lip"
    ),
    NOTCH = c(
        "Pituitary dwarfism",
        "Ectopic anterior pituitary gland",
        "Anterior pituitary agenesis",
        "Median cleft lip and palate",
        "Basal encephalocele",
        "Ethmoidal encephalocele",
        "Abnormal hard palate morphology",
        "Midnasal stenosis",
        "Abnormality of the midnasal cavity",
        "Abnormality of the nasal cavity",
        "Median cleft palate",
        "Cyclopia"
    ),
    FGF = c(
        "Pituitary dwarfism",
        "Ectopic anterior pituitary gland",
        "Anterior pituitary agenesis",
        "Median cleft lip and palate",
        "Hypoplastic iris stroma",
        "Basal encephalocele",
        "Ethmoidal encephalocele",
        "Midnasal stenosis",
        "Abnormality of the midnasal cavity",
        "Abnormality of the nasal cavity",
        "Median cleft palate",
        "Cyclopia",
        "Median cleft lip"
    ),
    WNT_PCP = c(
        "Abnormal hard palate morphology",
        "Median cleft lip",
        "Abnormal hypothalamus physiology",
        "Expressive language delay",
        "Abnormality of the diencephalon"
    )
)

pathway_hpo_annotated_genes <- lapply(names(group_hpo_terms), function(pw) {
    get_genes_from_hpo_terms(group_hpo_terms[[pw]], hpo_geneset, hpo_geneset_des)
})
names(pathway_hpo_annotated_genes) <- names(group_hpo_terms)

sapply(pathway_hpo_annotated_genes, length)


# ============================================================
# 5. Filtrage des seeds : uniquement celles directement connectees
#    a au moins un gene signal DU GROUPE OU a au moins un gene
#    HPE diagnostique (voisins dans le reseau RWR complet)
# ============================================================

filter_seeds_connected_to_signal_or_hpe <- function(seed_genes, signal_genes, hpe_genes, full_graph) {
    anchor_genes <- intersect(unique(c(signal_genes, hpe_genes)), V(full_graph)$name)
    if (length(anchor_genes) == 0) {
        return(character(0))
    }
    neighbor_genes <- unique(unlist(
        lapply(anchor_genes, function(g) neighbors(full_graph, g)$name)
    ))
    intersect(seed_genes, neighbor_genes)
}

pathway_connected_seed_genes <- lapply(names(pathway_seed_genes), function(pw) {
    filter_seeds_connected_to_signal_or_hpe(
        seed_genes   = pathway_seed_genes[[pw]],
        signal_genes = pathway_signal_genes[[pw]],
        hpe_genes    = HPE_diag_genes,
        full_graph   = g_full
    )
})
names(pathway_connected_seed_genes) <- names(pathway_seed_genes)

sapply(pathway_connected_seed_genes, length)


# ============================================================
# 6. Construction du reseau (dataframe node1/node2) par groupe
# ------------------------------------------------------------
# Le reseau du groupe est le sous-graphe du reseau RWR complet
# induit par l'union des genes signal et des seeds directement
# connectees a un gene signal OU a un gene HPE diagnostique.
# Pour chaque gene du groupe (node1), une ligne par voisin
# (node2) dans ce sous-graphe.
#
# node1          : gene primaire (attributs de la ligne)
# node2          : gene voisin de node1 dans le sous-reseau
# HPE            : node1 est un gene diagnostique HPE (T/F)
# seed           : node1 est une seed (directement connectee a
#                  un gene signal ou un gene HPE) dans ce groupe (T/F)
# seed_n_group   : nb de patients DU GROUPE ou node1 est une seed
# seed_n_cohort  : nb de patients DE LA COHORTE COMPLETE ou node1
#                  est une seed
# signal         : node1 est un gene signal dans ce groupe (T/F)
# HPO_annotated  : node1 est annote par au moins un des termes HPO
#                  fournis pour ce groupe (group_hpo_terms) (T/F)
# degree         : degre de node1 DANS LE RESEAU RWR COMPLET
# category       : 1 = gene HPE
#                  2 = seed non-HPE, directement liee a un signal
#                  3 = gene signal (non HPE, non seed)
#                  4 = le reste
# cor            : poids de l'arete node1-node2 dans le reseau RWR
#                  complet (coefficient de coexpression / correlation)
# corr_qual      : qualite de la correlation, evaluee PAR NODE1
#                  (sur |cor|) :
#                  2 = |cor| >= 0.5
#                  1 = |cor| <  0.5 MAIS c'est la correlation la
#                      plus elevee parmi les voisins de ce node1
#                  0 = sinon
# duplicate      : FALSE pour les aretes "originales", TRUE pour
#                  les aretes miroir (node1/node2 inverses),
#                  ajoutees pour que chaque gene apparaisse comme
#                  node1 -> a masquer dans Cytoscape si besoin
# ============================================================

build_group_network <- function(group_name,
                                signal_genes,
                                seed_genes,
                                hpo_annotated_genes,
                                full_graph,
                                full_degree_df,
                                seed_freq_cohort_df,
                                seed_freq_group_df,
                                hpe_diag_genes) {
    group_genes <- unique(c(signal_genes, seed_genes))
    group_genes <- intersect(group_genes, V(full_graph)$name)

    if (length(group_genes) < 2) {
        message("Groupe ", group_name, " : moins de 2 genes dans le reseau, ignore")
        return(NULL)
    }

    # Sous-graphe induit par les genes signal + seeds connectees
    subg <- induced_subgraph(full_graph, vids = group_genes)

    # Table d'attributs par gene (utilisee pour annoter node1,
    # que ce soit dans la version originale ou inversee)
    gene_attrs <- data.frame(gene = group_genes) %>%
        mutate(
            HPE           = gene %in% hpe_diag_genes,
            seed          = gene %in% seed_genes,
            signal        = gene %in% signal_genes,
            HPO_annotated = gene %in% hpo_annotated_genes
        ) %>%
        left_join(full_degree_df, by = "gene") %>%
        left_join(seed_freq_group_df, by = "gene") %>%
        left_join(seed_freq_cohort_df, by = "gene") %>%
        mutate(
            seed_n_group = replace_na(seed_n_group, 0),
            seed_n_cohort = replace_na(seed_n_cohort, 0),
            category = case_when(
                HPE ~ 1L,
                seed & !HPE ~ 2L,
                signal ~ 3L,
                TRUE ~ 4L
            )
        )

    # Liste d'aretes du sous-graphe (une ligne par arete, non orientee)
    # -> "weight" est l'attribut d'arete du graphe complet (coex_case),
    #    conserve tel quel par induced_subgraph()
    edges <- igraph::as_data_frame(subg, what = "edges") %>%
        dplyr::select(node1 = from, node2 = to, cor = weight)

    if (nrow(edges) == 0) {
        message("Groupe ", group_name, " : aucune arete dans le sous-reseau")
        return(NULL)
    }

    edges <- edges %>%
        left_join(gene_attrs, by = c("node1" = "gene")) %>%
        mutate(duplicate = FALSE)

    # Version miroir : node1 <-> node2 inverses, ATTRIBUTS
    # RECALCULES pour le nouveau node1 (ex-node2). cor reste la
    # meme valeur car la correlation node1-node2 est symetrique.
    edges_mirror <- igraph::as_data_frame(subg, what = "edges") %>%
        dplyr::select(node1 = to, node2 = from, cor = weight) %>%
        left_join(gene_attrs, by = c("node1" = "gene")) %>%
        mutate(duplicate = TRUE)

    bind_rows(edges, edges_mirror) %>%
        dplyr::group_by(node1) %>%
        dplyr::mutate(
            max_abs_cor_node1 = max(abs(cor)),
            corr_qual = dplyr::case_when(
                abs(cor) >= 0.5 ~ 2L,
                abs(cor) == max_abs_cor_node1 ~ 1L,
                TRUE ~ 0L
            )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(node1, node2, HPE, seed, seed_n_group, seed_n_cohort, signal, HPO_annotated, degree, category, cor, corr_qual, duplicate)
}

pathway_networks <- lapply(names(pathway_signal_genes), function(pw) {
    build_group_network(
        group_name = pw,
        signal_genes = pathway_signal_genes[[pw]],
        seed_genes = pathway_connected_seed_genes[[pw]],
        hpo_annotated_genes = pathway_hpo_annotated_genes[[pw]],
        full_graph = g_full,
        full_degree_df = full_degree,
        seed_freq_cohort_df = seed_freq_cohort,
        seed_freq_group_df = pathway_seed_freq_group[[pw]],
        hpe_diag_genes = HPE_diag_genes
    )
})
names(pathway_networks) <- names(pathway_signal_genes)

sapply(pathway_networks, function(df) if (is.null(df)) 0 else nrow(df))


# ============================================================
# 7. Sauvegarde (un fichier par groupe, import direct Cytoscape)
# ------------------------------------------------------------
# quote = FALSE : pas de guillemets autour des symboles de genes
# ============================================================
dir.create("article/networks", recursive = TRUE, showWarnings = FALSE)

invisible(lapply(names(pathway_networks), function(pw) {
    df <- pathway_networks[[pw]]
    if (is.null(df)) {
        return(invisible(NULL))
    }
    write.csv(
        df,
        file.path("article/networks", paste0(pw, "_signal_seed_network.csv")),
        row.names = FALSE,
        quote     = FALSE
    )
}))
