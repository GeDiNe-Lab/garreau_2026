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

# Setting working directory
rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-3) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()

# Loading custom functions
source("article/scripts/custom_fct.R")

load("article/transit_data/RWR_results_case.RData")
HPE_results <- as.matrix(RWR_results_z)
rm(RWR_results_z)
load("article/transit_data/RWR_iron_results_case.RData")
FER_results <- as.matrix(RWR_results_z)
rm(RWR_results_z)
load("article/transit_data/data_3p.RData")
load("article/transit_data/iron_patient_variants.RData")
FER_patient_variants <- patient_variants
rm(patient_variants)
load("article/transit_data/patient_variants.RData")
HPE_patient_variants <- patient_variants
rm(patient_variants)
load("article/transit_data/permutation_scores.RData")


coex_genes <- rownames(data_3p$coex_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")
HPE_seeds <- lapply(HPE_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(HPE_seeds) <- colnames(HPE_results)

RWR_score_matrix <- HPE_results
seed_list <- HPE_seeds

RWR_score_matrix_corr <- RWR_score_matrix

for (p in colnames(RWR_score_matrix_corr)) {
    patient_scores <- RWR_score_matrix_corr[, p]
    patient_seeds <- seed_list[[p]]

    if (is.null(patient_seeds) || length(patient_seeds) == 0) next

    seed_idx <- rownames(RWR_score_matrix_corr) %in% patient_seeds
    non_seed_median <- median(patient_scores[!seed_idx], na.rm = TRUE)

    RWR_score_matrix_corr[seed_idx, p] <- non_seed_median
}

# 1. Discard patients with NA scores (whole column NA)
valid_patients <- colnames(RWR_score_matrix_corr)[!is.na(colMeans(RWR_score_matrix_corr, na.rm = FALSE))]
RWR_score_matrix_f <- RWR_score_matrix_corr[, valid_patients, drop = FALSE]
seed_list_f <- seed_list[valid_patients]

# 2. Mean score per gene across valid patients, pick top 50
mean_scores <- rowMeans(RWR_score_matrix_f, na.rm = TRUE)
top50_genes <- sort(mean_scores, decreasing = TRUE) %>%
    head(50) %>%
    names()

# 3. Build long-format df: one row per gene x patient, for the top 50 genes only
top50_long_df <- RWR_score_matrix_f[top50_genes, , drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = "patient", values_to = "score")

# 4. Number of patients each gene is a seed for (0 if none), among valid patients
seed_counts <- sapply(top50_genes, function(g) {
    sum(sapply(seed_list_f, function(seeds) g %in% seeds))
})
seed_counts_df <- data.frame(gene = top50_genes, n_patient_seed = seed_counts)

top50_long_df <- top50_long_df %>%
    left_join(seed_counts_df, by = "gene")

# 5. Order genes by mean score (descending) for the plot
top35_long_df <- top35_long_df %>%
    mutate(gene = fct_reorder(gene, mean_scores[as.character(gene)]))

# 5b. Flag points that are outliers relative to their gene's box (default ggplot/Tukey rule: outside 1.5*IQR from Q1/Q3)
top35_long_df <- top35_long_df %>%
    group_by(gene) %>%
    mutate(
        q1 = quantile(score, 0.25, na.rm = TRUE),
        q3 = quantile(score, 0.75, na.rm = TRUE),
        iqr = q3 - q1,
        lower_whisker = q1 - 1.5 * iqr,
        upper_whisker = q3 + 1.5 * iqr,
        is_outlier = score < lower_whisker | score > upper_whisker
    ) %>%
    ungroup()

outliers_df <- top35_long_df %>% filter(is_outlier)

label_df <- top35_long_df %>%
    group_by(gene, n_patient_seed) %>%
    summarise(y_pos = max(score, na.rm = TRUE), .groups = "drop")

# 6. Box plot + jitter only for outlier points
ggplot(top35_long_df, aes(y = gene, x = score)) +
    geom_boxplot(fill = "#2E86AB", alpha = 0.4, outlier.shape = NA, width = 0.6) +
    geom_jitter(
        data = outliers_df,
        width = 0.15, alpha = 0.6, size = 1.4, color = "#A23B72"
    ) +
    geom_text(
        data = label_df,
        aes(y = gene, x = y_pos, label = n_patient_seed), size = 3, fontface = "bold", inherit.aes = FALSE
    ) +
    labs(
        title = "Distribution of RWR scores across patients — top 50 genes",
        subtitle = "Seed scores replaced by non-seed median per patient · jitter shows outliers only · number next to box = seed count",
        y = "Gene",
        x = "RWR score"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 9),
        panel.grid.minor = element_blank()
    )



ggsave("article/figures/top50_genes_score_boxplot_corrected.png", width = 12, height = 7, dpi = 300)
