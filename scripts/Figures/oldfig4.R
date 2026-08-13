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

HPE_seeds <- lapply(HPE_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(HPE_seeds) <- colnames(HPE_results)


# ============================================================
# 1. Genes signal, cohorte HPE complete (tous patients)
# ============================================================

cohort_result <- process_group(
    patient_ids      = names(HPE_seeds),
    RWR_results      = HPE_results,
    seeds            = HPE_seeds,
    permutation_bank = bank
)

signal_genes <- cohort_result$gene[cohort_result$signal_final]


# ============================================================
# 2. Frequence des seeds, cohorte HPE complete
# ------------------------------------------------------------
#  - seed_n : nb de patients (cohorte HPE complete) ou le gene
#             est une seed
# ============================================================

seed_freq_cohort <- tibble(gene = unlist(HPE_seeds)) %>%
    dplyr::count(gene, name = "seed_n")

all_seed_genes <- unique(unlist(HPE_seeds))


# ============================================================
# 3. Reseau complet (degre) et filtrage des seeds : uniquement
#    celles directement connectees a au moins un gene signal
#    (voisins dans le reseau RWR complet)
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

filter_seeds_connected_to_signal <- function(seed_genes, signal_genes, full_graph) {
    anchor_genes <- intersect(unique(signal_genes), V(full_graph)$name)
    if (length(anchor_genes) == 0) {
        return(character(0))
    }
    neighbor_genes <- unique(unlist(
        lapply(anchor_genes, function(g) neighbors(full_graph, g)$name)
    ))
    intersect(seed_genes, neighbor_genes)
}

connected_seed_genes <- filter_seeds_connected_to_signal(
    seed_genes   = all_seed_genes,
    signal_genes = signal_genes,
    full_graph   = g_full
)


# ============================================================
# 4. Construction du reseau (dataframe node1/node2)
# ------------------------------------------------------------
# Le reseau est le sous-graphe du reseau RWR complet induit par
# l'union des genes signal et des seeds directement connectees
# a un gene signal.
# Pour chaque gene (node1), une ligne par voisin (node2) dans ce
# sous-graphe.
#
# node1     : gene primaire (attributs de la ligne)
# node2     : gene voisin de node1 dans le sous-reseau
# signal    : node1 est un gene signal de la cohorte HPE (T/F)
# seed      : node1 est une seed chez au moins un patient HPE (T/F)
# seed_n    : nb de patients HPE ou node1 est une seed
# condition : statut de node1 -> "seed", "signal" ou "seed+signal"
#             (un gene peut etre a la fois seed et signal)
# degree    : degre de node1 DANS LE RESEAU RWR COMPLET
# cor       : poids de l'arete node1-node2 dans le reseau RWR
#             complet (coefficient de coexpression / correlation)
# duplicate : FALSE pour les aretes "originales", TRUE pour les
#             aretes miroir (node1/node2 inverses), ajoutees pour
#             que chaque gene apparaisse comme node1 (annotations
#             completes) -> a masquer dans Cytoscape une fois
#             l'import fait (filtrer duplicate == TRUE)
# ============================================================

build_full_cohort_network <- function(signal_genes,
                                      seed_genes,
                                      full_graph,
                                      full_degree_df,
                                      seed_freq_cohort_df,
                                      all_seed_genes) {
    group_genes <- unique(c(signal_genes, seed_genes))
    group_genes <- intersect(group_genes, V(full_graph)$name)

    if (length(group_genes) < 2) {
        message("Moins de 2 genes dans le reseau, impossible de construire le graphe")
        return(NULL)
    }

    subg <- induced_subgraph(full_graph, vids = group_genes)

    gene_attrs <- data.frame(gene = group_genes) %>%
        mutate(
            signal = gene %in% signal_genes,
            seed   = gene %in% all_seed_genes
        ) %>%
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

    edges <- igraph::as_data_frame(subg, what = "edges") %>%
        dplyr::select(node1 = from, node2 = to, cor = weight)

    if (nrow(edges) == 0) {
        message("Aucune arete dans le sous-reseau")
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
        dplyr::select(node1, node2, signal, seed, seed_n, condition, degree, cor, duplicate)
}

full_cohort_network <- build_full_cohort_network(
    signal_genes         = signal_genes,
    seed_genes           = connected_seed_genes,
    full_graph           = g_full,
    full_degree_df       = full_degree,
    seed_freq_cohort_df  = seed_freq_cohort,
    all_seed_genes       = all_seed_genes
)

nrow(full_cohort_network)


# ============================================================
# 5. Table recap (niveau noeud, une ligne par gene)
# ------------------------------------------------------------
# Reprend toutes les infos de node1 dans full_cohort_network,
# dedupliquees par gene (independamment de duplicate et des
# aretes) -> une ligne par gene du reseau.
# ============================================================

node_recap <- full_cohort_network %>%
    dplyr::select(node1, signal, seed, seed_n, condition, degree) %>%
    dplyr::distinct(node1, .keep_all = TRUE) %>%
    dplyr::rename(gene = node1) %>%
    dplyr::arrange(desc(seed_n), desc(degree))

nrow(node_recap)


# ============================================================
# 6. Sauvegarde (import direct Cytoscape + table recap)
# ------------------------------------------------------------
# quote = FALSE : pas de guillemets autour des symboles de genes
# ============================================================
dir.create("article/networks", recursive = TRUE, showWarnings = FALSE)
dir.create("article/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(
    full_cohort_network,
    file.path("article/networks", "HPE_full_cohort_signal_seed_network.csv"),
    row.names = FALSE,
    quote     = FALSE
)

write.csv(
    node_recap,
    file.path("article/tables", "HPE_full_cohort_signal_seed_node_recap.csv"),
    row.names = FALSE,
    quote     = FALSE
)
