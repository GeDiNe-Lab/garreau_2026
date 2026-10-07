# Loading packages and functions
library(Matrix)
library(tidyverse)

setwd(paste(here::here(), "projects", sep = "/"))


# Loading custom functions
source("article/scripts/custom_fct.R")

# loading single cell data from Zeng et al from week3,week4 and week5
sc_counts <- readMM("panasenkava_2024/scdata/week345_wholebody.mtx") %>% t()

# keeping only cells related to Central Nervous System (SNC)
cell_ids <- read.table("panasenkava_2024/scdata/indices_week345_wholebody.csv", header = TRUE, sep = ",")
genes <- read.table("panasenkava_2024/scdata/genes_week345_wholebody.tsv", header = FALSE)$V1
rownames(sc_counts) <- genes

formated_ids <- c(1:nrow(cell_ids)) %>% sapply(function(i) {
    c(str_split(cell_ids[i, ]$index, "-")[[1]][1:2], cell_ids[i, ]$week_stage) %>%
        paste(collapse = "-") %>%
        return()
})
colnames(sc_counts) <- formated_ids

# keeping only cells tagged as Forebrain progenitor (FB) and/or Brain Neureon depending on the week
authors_meta <- read.csv("panasenkava_2024/scdata/week345_wholebody_metadata.csv")
authors_meta_f <- filter(
    authors_meta,
    (celltype_region == "FB" & week_stage == "W3-1") |
        (celltype_region %in% c("Brain Neuron", "FB") & week_stage == "W4-1") |
        (celltype_region == "Brain Neuron" & week_stage == "W5-1")
)

sc_counts_f1 <- sc_counts[, authors_meta_f$barcode]

nrow(sc_counts_f1)

# counts: genes x cells
n_cells <- ncol(sc_counts_f1)
pct_detected <- Matrix::rowSums(sc_counts_f1 > 0) / n_cells # fraction of cells expressing the gene
total_counts <- Matrix::rowSums(sc_counts_f1) # total UMI per gene

keep <- pct_detected >= 0.005 & total_counts >= 25
sc_counts_f <- sc_counts_f1[keep, ]
dim(sc_counts_f)
0.005 * 4659

sc_counts_f <- sc_counts_f1[which(rowSums(sc_counts_f1) > 25), ]
nrow(sc_counts_f)
sc_genes <- rownames(sc_counts_f)

save(sc_genes, file = "article/transit_data/sc_genes.RData")
