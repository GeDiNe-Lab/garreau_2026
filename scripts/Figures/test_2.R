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

load("article/transit_data/iop_patient_variants.RData")
IOP_patient_variants <- patient_variants
rm(patient_variants)
load("article/transit_data/patient_variants.RData")
HPE_patient_variants <- patient_variants
rm(patient_variants)

# --- Coex network degree centrality ---
g <- graph_from_adjacency_matrix(
    data_3p$coex_case,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
)

centrality_matrix <- data.frame(
    node = V(g)$name %>% gene_converter(from = "ENSEMBL", to = "SYMBOL"),
    degree = igraph::degree(g)
)
net_genes <- centrality_matrix$node

compute_seed_freq <- function(patient_variants, net_genes) {
    all_genes <- patient_variants %>%
        lapply(function(patient_df) {
            intersect(unique(patient_df$Gene.refGeneWithVer), net_genes)
        }) %>%
        unlist(use.names = FALSE)

    as.data.frame(table(all_genes), stringsAsFactors = FALSE) %>%
        dplyr::rename(node = all_genes, Freq = Freq)
}

HPE_seeds_freq <- compute_seed_freq(HPE_patient_variants, net_genes)
IOP_seeds_freq <- compute_seed_freq(IOP_patient_variants, net_genes)


# --- Extract per-gene pLI from patient variant tables ---
# pLI is a gene-level annotation (gnomAD) carried per-variant in each patient's
# table, so it should be constant across patients/cohorts for a given gene.
# We pool variants from both cohorts, restrict to network genes, and keep one
# pLI value per gene (first non-NA). Genes with inconsistent pLI values across
# rows are flagged so this assumption can be checked.
extract_gene_pLI <- function(patient_variants_list, net_genes) {
    patient_variants_list <- lapply(patient_variants_list, function(patient_df) {
        patient_df$pLI <- cleaning_pLI(patient_df$pLI)
        patient_df
    })

    bind_rows(patient_variants_list) %>%
        filter(Gene.refGeneWithVer %in% net_genes) %>%
        dplyr::select(node = Gene.refGeneWithVer, pLI) %>%
        filter(!is.na(pLI))
}

gene_pLI_HPE <- extract_gene_pLI(HPE_patient_variants, net_genes)
gene_pLI_IOP <- extract_gene_pLI(IOP_patient_variants, net_genes)

gene_pLI_all <- bind_rows(gene_pLI_HPE, gene_pLI_IOP)

# Check for genes with inconsistent pLI values across variants/patients
pLI_inconsistencies <- gene_pLI_all %>%
    distinct(node, pLI) %>%
    group_by(node) %>%
    filter(n() > 1) %>%
    ungroup()

if (nrow(pLI_inconsistencies) > 0) {
    warning(
        n_distinct(pLI_inconsistencies$node),
        " gene(s) have inconsistent pLI values across variants; the first",
        " non-NA value encountered will be used. See `pLI_inconsistencies`."
    )
}

gene_pLI <- gene_pLI_all %>%
    distinct(node, .keep_all = TRUE)

# --- Merge into one row per gene, with HPE_Freq / IOP_Freq columns ---
# A gene that only carries variants in one cohort gets Freq = 0 for the other,
# instead of being duplicated across two rows.
seed_centrality_df <- centrality_matrix %>%
    filter(node %in% union(HPE_seeds_freq$node, IOP_seeds_freq$node)) %>%
    left_join(HPE_seeds_freq %>% dplyr::rename(HPE_Freq = Freq), by = "node") %>%
    left_join(IOP_seeds_freq %>% dplyr::rename(IOP_Freq = Freq), by = "node") %>%
    mutate(
        HPE_Freq = replace_na(HPE_Freq, 0L),
        IOP_Freq = replace_na(IOP_Freq, 0L)
    ) %>%
    left_join(gene_pLI, by = "node") %>%
    relocate(node, HPE_Freq, IOP_Freq, degree, pLI)

# --- Cohort summary: unique seeds vs number of patients ---
n_patients_HPE <- length(HPE_patient_variants)
n_patients_IOP <- length(IOP_patient_variants)
n_seeds_HPE <- sum(seed_centrality_df$HPE_Freq > 0)
n_seeds_IOP <- sum(seed_centrality_df$IOP_Freq > 0)

cohort_summary <- data.frame(
    cohort = c("HPE", "IOP"),
    n_patients = c(n_patients_HPE, n_patients_IOP),
    n_unique_seeds = c(n_seeds_HPE, n_seeds_IOP)
) %>%
    mutate(seeds_per_patient = round(n_unique_seeds / n_patients, 2))
cohort_summary

subtitle_text <- paste0(
    "HPE: ", n_seeds_HPE, " unique seeds / ", n_patients_HPE, " patients   |   ",
    "IOP: ", n_seeds_IOP, " unique seeds / ", n_patients_IOP, " patients"
)

# --- Long format (plotting only) ---
# seed_centrality_df stays one-row-per-gene; this reshaped copy is only used
# to feed the cohort-coloured scatterplots below.
seed_centrality_long <- seed_centrality_df %>%
    pivot_longer(
        cols = c(HPE_Freq, IOP_Freq),
        names_to = "cohort",
        values_to = "Freq"
    ) %>%
    mutate(cohort = str_remove(cohort, "_Freq")) %>%
    filter(Freq > 0)

# --- Dotplot: seed frequency (x) vs network degree centrality (y) ---
p_seed_centrality <- ggplot(seed_centrality_long, aes(x = Freq, y = degree, color = cohort)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c(HPE = "#D55E00", IOP = "#0072B2")) +
    labs(
        title = "Seed gene frequency vs network degree centrality",
        subtitle = subtitle_text,
        x = "Seed frequency (number of patients carrying the variant)",
        y = "Degree centrality (coex network)",
        color = "Cohort"
    ) +
    custom_theme()



# --- Dotplot: seed frequency (x) vs network degree centrality (y), colored by pLI ---
p_seed_centrality_pLI <- ggplot(seed_centrality_long, aes(x = Freq, y = degree, color = pLI, shape = cohort)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_viridis_c(na.value = "grey80") +
    labs(
        title = "Seed gene frequency vs network degree centrality, coloured by pLI",
        subtitle = subtitle_text,
        x = "Seed frequency (number of patients carrying the variant)",
        y = "Degree centrality (coex network)",
        color = "pLI",
        shape = "Cohort"
    ) +
    custom_theme()


# --- Fonction: récapitulatif des variants pour un gène donné, dans les deux cohortes ---
get_gene_variants <- function(gene, HPE_patient_variants, IOP_patient_variants) {
    extract_variants <- function(patient_variants, cohort_label) {
        names(patient_variants) <- paste0("Patient_", seq_along(patient_variants))
        patient_variants %>%
            lapply(function(patient_df) {
                filter(patient_df, Gene.refGeneWithVer == gene)
            }) %>%
            bind_rows(.id = "patient") %>%
            mutate(cohort = cohort_label)
    }

    bind_rows(
        extract_variants(HPE_patient_variants, "HPE"),
        extract_variants(IOP_patient_variants, "IOP")
    ) %>%
        relocate(cohort, patient, .before = everything())
}


seed_centrality_df %>% View()

View(get_gene_variants("LRP2", HPE_patient_variants, IOP_patient_variants))








# --- Number of seed genes per patient, distribution ---
seeds_per_patient <- bind_rows(
    tibble(
        n_seeds = sapply(HPE_patient_variants, function(df) {
            length(intersect(unique(df$Gene.refGeneWithVer), net_genes))
        }),
        cohort = "HPE"
    ),
    tibble(
        n_seeds = sapply(IOP_patient_variants, function(df) {
            length(intersect(unique(df$Gene.refGeneWithVer), net_genes))
        }),
        cohort = "IOP"
    )
)
seeds_per_patient %>% View()
ggplot(seeds_per_patient, aes(x = n_seeds, fill = cohort)) +
    geom_histogram(position = "identity", alpha = 0.6, binwidth = 2, color = "black") +
    scale_fill_manual(values = c(HPE = "#D55E00", IOP = "#0072B2")) +
    labs(
        title = "Distribution of seed genes per patient",
        x = "Number of seed genes",
        y = "Number of patients",
        fill = "Cohort"
    ) +
    custom_theme()





# ============================================================
# GO-BP / HPO ORA on seed genes, per cohort
# ============================================================
# Background (referenceGene) = all network genes, consistent with how
# seed_centrality_df was built (seeds are already restricted to net_genes).

HPE_seed_genes <- seed_centrality_df %>%
    filter(HPE_Freq > 0) %>%
    pull(node)

IOP_seed_genes <- seed_centrality_df %>%
    filter(IOP_Freq > 0) %>%
    pull(node)

# --- HPE ---
HPE_seeds_GO <- run_ora_gobp(
    genelist = HPE_seed_genes,
    referenceGene = net_genes,
    interestGeneType = "genesymbol",
    outputDirectory = "article/enrichment",
    projectName = "HPE_seed_centrality_GO",
    fdrThr = 0.05,
    minNum = 5,
    maxNum = 2000,
    noRedundant = TRUE,
    is.output = TRUE
)

HPE_seeds_HPO <- run_ora_hpo(
    genelist = HPE_seed_genes,
    referenceGene = net_genes,
    interestGeneType = "genesymbol",
    outputDirectory = "article/enrichment",
    projectName = "HPE_seed_centrality_HPO",
    fdrThr = 0.05,
    minNum = 5,
    maxNum = 2000,
    is.output = TRUE
)

HPE_GO_plot <- plot_ora_dotplot(HPE_seeds_GO, n = 20, title = "HPE seeds - Top 20 GO-BP")
HPE_HPO_plot <- plot_ora_dotplot(HPE_seeds_HPO, n = 20, title = "HPE seeds - Top 20 HPO")

print(HPE_GO_plot)
print(HPE_HPO_plot)

ggsave("article/figures/HPE_seed_centrality_GO.png", HPE_GO_plot, width = 10, height = 6, dpi = 300)
ggsave("article/figures/HPE_seed_centrality_HPO.png", HPE_HPO_plot, width = 10, height = 6, dpi = 300)

# --- IOP ---
IOP_seeds_GO <- run_ora_gobp(
    genelist = IOP_seed_genes,
    referenceGene = net_genes,
    interestGeneType = "genesymbol",
    outputDirectory = "article/enrichment",
    projectName = "IOP_seed_centrality_GO",
    fdrThr = 0.05,
    minNum = 5,
    maxNum = 2000,
    noRedundant = TRUE,
    is.output = TRUE
)

IOP_seeds_HPO <- run_ora_hpo(
    genelist = IOP_seed_genes,
    referenceGene = net_genes,
    interestGeneType = "genesymbol",
    outputDirectory = "article/enrichment",
    projectName = "IOP_seed_centrality_HPO",
    fdrThr = 0.05,
    minNum = 5,
    maxNum = 2000,
    is.output = TRUE
)

IOP_GO_plot <- plot_ora_dotplot(IOP_seeds_GO, n = 20, title = "IOP seeds - Top 20 GO-BP")
IOP_HPO_plot <- plot_ora_dotplot(IOP_seeds_HPO, n = 20, title = "IOP seeds - Top 20 HPO")

print(IOP_GO_plot)
print(IOP_HPO_plot)

ggsave("article/figures/IOP_seed_centrality_GO.png", IOP_GO_plot, width = 10, height = 6, dpi = 300)
ggsave("article/figures/IOP_seed_centrality_HPO.png", IOP_HPO_plot, width = 10, height = 6, dpi = 300)
