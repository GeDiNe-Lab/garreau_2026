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

HPE_results <- HPE_results[, which(!is.na(colMeans(HPE_results)))]

c(1:ncol(HPE_results)) %>%
    lapply(function(i) {
        x <- HPE_results[, i]
        true_mean <- mean(x[which(x < 5)])
        x[which(x >= 5)] <- true_mean
        return(x)
    }) %>%
    Reduce(cbind, .) -> HPE_score_matrix
colnames(HPE_score_matrix) <- colnames(HPE_results)
rownames(HPE_score_matrix) <- rownames(HPE_results)

FER_results <- FER_results[, which(!is.na(colMeans(FER_results)))]

c(1:ncol(FER_results)) %>%
    lapply(function(i) {
        x <- FER_results[, i]
        true_mean <- mean(x[which(x < 5)])
        x[which(x >= 5)] <- true_mean
        return(x)
    }) %>%
    Reduce(cbind, .) -> FER_score_matrix
colnames(FER_score_matrix) <- colnames(FER_results)
rownames(FER_score_matrix) <- rownames(FER_results)

hpe_means <- rowMeans(HPE_score_matrix) %>%
    sort(decreasing = TRUE)
plot(hpe_means)

# filter HPE_score_matrix and FER_score_matrix for HPE_score_matrix top 1000 genes
top_genes <- names(hpe_means)[1:1000]
HPE_score_matrix_f <- HPE_score_matrix[top_genes, ]
FER_score_matrix_f <- FER_score_matrix[top_genes, ]

# Wilcoxon test for genes in HPE compared to genes in FER
common_genes <- intersect(rownames(HPE_score_matrix_f), rownames(FER_score_matrix_f))
wilcox_results <- data.frame(
    gene = common_genes,
    p.value = NA_real_,
    statistic = NA_real_,
    HPE_median = NA_real_,
    FER_median = NA_real_,
    adj.p.value = NA_real_,
    stringsAsFactors = FALSE
)

for (i in seq_along(common_genes)) {
    g <- common_genes[i]
    hpe <- HPE_score_matrix_f[g, ]
    fer <- FER_score_matrix_f[g, ]
    if (all(is.na(hpe)) || all(is.na(fer))) {
        next
    }
    test <- wilcox.test(hpe, fer, paired = FALSE, alternative = "greater", exact = FALSE)
    wilcox_results$p.value[i] <- test$p.value
    wilcox_results$statistic[i] <- as.numeric(test$statistic)
    wilcox_results$HPE_median[i] <- median(hpe, na.rm = TRUE)
    wilcox_results$FER_median[i] <- median(fer, na.rm = TRUE)
}

wilcox_results$adj.p.value <- p.adjust(wilcox_results$p.value, method = "fdr")
# add fold change column
wilcox_results <- wilcox_results %>%
    mutate(fold_change = HPE_median / FER_median)
wilcox_results <- wilcox_results %>% arrange(p.value)
View(wilcox_results)

wilcox_results$adj.p.value %>% sort()
# Output directory
out_dir <- "article/figures/top1000_HPE_distribs"

for (dpy_gene in rownames(HPE_score_matrix_f)) {
    dpy_hpe <- HPE_score_matrix_f[dpy_gene, ]
    dpy_fer <- FER_score_matrix_f[dpy_gene, ]

    dpy_data <- data.frame(
        score = c(dpy_hpe, dpy_fer),
        group = c(rep("HPE", length(dpy_hpe)), rep("FER", length(dpy_fer)))
    )

    p <- ggplot(dpy_data, aes(x = score, y = group)) +
        geom_boxplot(width = 0.3, outlier.shape = NA) +
        geom_jitter(height = 0.15, alpha = 0.75) +
        ggridges::geom_density_ridges(
            aes(fill = group),
            alpha = 0.2,
            scale = 0.9,
            rel_min_height = 0.01,
            position = position_nudge(y = -0.3)
        ) +
        labs(title = paste("Distribution of scores for", dpy_gene), x = "Score", y = NULL) +
        theme_light() +
        theme(legend.position = "none")

    ggsave(
        filename = file.path(out_dir, paste0(dpy_gene, ".png")),
        plot = p,
        width = 7, height = 4, dpi = 150
    )
}
