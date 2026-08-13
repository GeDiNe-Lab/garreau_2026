# Loading packages and functions
library(Matrix)
library(tidyverse)
library(parallel)
library(tibble)
library(ggplot2)
library(tidyr)
library(scales)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(purrr)

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

load("article/transit_data/patient_variants.RData")
load("article/transit_data/data_3p.RData")

patient_variants[[1]] %>% head()

combined_patient_variants <- bind_rows(patient_variants, .id = "patient")
combined_patient_variants$patient <- as.integer(combined_patient_variants$patient)

View(combined_patient_variants)

# Create dataframe with gene symbol, patient number and counts per patient
df_gene_patient_counts <- combined_patient_variants %>%
    group_by(Gene.refGeneWithVer, patient) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(Gene.refGeneWithVer, patient)

# create piechart of genes for count of variant
df_gene_counts <- combined_patient_variants %>%
    group_by(Gene.refGeneWithVer) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(desc(count)) %>%
    filter(Gene.refGeneWithVer %in% gene_converter(rownames(data_3p$coex_case), "ENSEMBL", "SYMBOL"))

df_gene_counts$Gene.refGeneWithVer %>% write(file = "article/transit_data/patient_seeds.csv")

# add a column ENSEMBL ID to the df_gene_counts
df_gene_counts <- df_gene_counts %>%
    mutate(ENSEMBL = gene_converter(Gene.refGeneWithVer, "SYMBOL", "ENSEMBL"))


# Compute number of unique genes per patient
unique_gene_counts <- patient_variants %>%
    sapply(function(patient) {
        filter(patient, Gene.refGeneWithVer %in% gene_converter(rownames(data_3p$coex_case), "ENSEMBL", "SYMBOL"))$Gene.refGeneWithVer %>%
            unique() %>%
            length()
    })
unique_gene_counts

# Convert to data frame for ggplot
df_genes <- data.frame(
    patient = as.character(c(1:199)),
    n_genes = unique_gene_counts
)

View(df_genes)

p <- ggplot(df_genes, aes(x = n_genes)) +
    geom_histogram(binwidth = 2, fill = "steelblue", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(
        x = "Number of seed per patient",
        y = "Count of patients",
        title = "Distribution of RWR seed genes across patients"
    )
p


# Compute number of unique genes per patient
unique_gene_counts <- patient_variants %>%
    sapply(function(patient) {
        filter(patient, Gene.refGeneWithVer %in% gene_converter(rownames(data_3p$coex_case), "ENSEMBL", "SYMBOL"))$Gene.refGeneWithVer %>%
            unique() %>%
            length()
    })
unique_gene_counts

# Convert to data frame for ggplot
df_genes <- data.frame(
    patient = as.character(c(1:199)),
    n_genes = unique_gene_counts
)

View(df_genes)

p <- ggplot(df_genes, aes(x = n_genes)) +
    geom_histogram(binwidth = 2, fill = "steelblue", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(
        x = "Number of seed per patient",
        y = "Count of patients",
        title = "Distribution of RWR seed genes across patients"
    )
p
