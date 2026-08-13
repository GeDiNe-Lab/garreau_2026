# Loading packages and functions
library(Matrix)
library(tidyverse)
library(parallel)
library(tibble)
library(EGAD)

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

load("article/transit_data/data_3p.RData")

hpo <- read.csv("article/phenotype_to_genes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

raw_RWR_net <- abs(data_3p$coex_consensus) %>% as.matrix()
#  Removing lower correlations interactions
raw_RWR_net[abs(raw_RWR_net) < 0.7] <- 0
#  Removing genes with no interactions
keep <- which(rowSums(raw_RWR_net != 0) > 1)
RWR_net <- raw_RWR_net[keep, keep]
rownames(RWR_net) <- gene_converter(rownames(RWR_net), "ENSEMBL", "SYMBOL")
colnames(RWR_net) <- gene_converter(colnames(RWR_net), "ENSEMBL", "SYMBOL")

load("article/transit_data/hpo_matrix.RData")

hpo_matrix <- t(hpo_matrix)

neighvote_results <- run_GBA(RWR_net, hpo_matrix)

neighvote_df <- neighvote_results[[1]] %>% as.data.frame()

hpo_names <- hpo$hpo_name[match(rownames(neighvote_df), hpo$hpo_id)] %>% as.vector()

neighvote_df$name <- hpo_names

ORA_results <- read.table("article/transit_data/hpo_term_clusters.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ORA_results %>% head()
neighvote_df_ora <- neighvote_df[ORA_results$HPO_ID, ]
neighvote_df_ora %>% View()
