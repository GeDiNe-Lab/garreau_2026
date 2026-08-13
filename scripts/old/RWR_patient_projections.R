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

genes <- rownames(data_3p$coex_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

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
HPE_score_matrix %>% head()

# PCA of patients (columns of score matrix)
pca_patients <- prcomp(t(HPE_score_matrix))

# Percent variance explained
var_exp <- (pca_patients$sdev^2) / sum(pca_patients$sdev^2) * 100

# Prepare data frame for plotting
pca_df <- as.data.frame(pca_patients$x)
pca_df$sample <- rownames(pca_df)

# Basic scatter plot of PC1 vs PC2
pc_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
    geom_point(size = 2) +
    geom_text(hjust = 0.5, vjust = -1, size = 3) +
    xlab(paste0("PC1 (", round(var_exp[1], 1), "%)")) +
    ylab(paste0("PC2 (", round(var_exp[2], 1), "%)")) +
    theme_minimal()
pc_plot

# Ensure output directory exists and save plot
out_dir <- "article/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
ggsave(file.path(out_dir, "HPE_patients_PCA.png"), pc_plot, width = 7, height = 5)

print(pc_plot)
