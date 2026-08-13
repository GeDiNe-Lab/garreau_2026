# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(tibble)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(readr)
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
load("article/transit_data/data_3p.RData")
load("article/transit_data/patient_variants.RData")

raw_RWR_net <- abs(data_3p$coex_consensus) %>% as.matrix()
#  Removing lower correlations interactions
raw_RWR_net[abs(raw_RWR_net) < 0.1] <- 0
#  Removing genes with no interactions
keep <- which(rowSums(raw_RWR_net != 0) > 1)
RWR_net <- raw_RWR_net[keep, keep]

nrow(data_3p$coex_consensus)
nrow(RWR_net)
write.csv(colnames(RWR_net), file = "article/transit_data/network_genes.csv")

# convert adjacency matrix to edge list dataframe (exclude zeros and self-loops, keep correlation)
idx <- which(upper.tri(RWR_net) & RWR_net != 0, arr.ind = TRUE)

network_df <- data.frame(
    gene1 = rownames(RWR_net)[idx[, 1]],
    gene2 = colnames(RWR_net)[idx[, 2]],
    correlation = RWR_net[idx]
)

library(igraph)

# build undirected graph from edge list
g <- graph_from_data_frame(network_df[, c("gene1", "gene2")], directed = FALSE)

# find largest connected component
cmp <- components(g)
largest_id <- which.max(cmp$csize)
largest_nodes <- V(g)$name[cmp$membership == largest_id]

# keep only edges where both endpoints belong to the largest component
network_df <- network_df %>%
    filter(gene1 %in% largest_nodes & gene2 %in% largest_nodes)
# optional: sort by absolute correlation descending
network_df <- network_df[order(-abs(network_df$correlation)), ]

network_df$link <- rep("base", nrow(network_df))


gene1_notin_gene2 <- which(network_df$gene1 %in% setdiff(unique(network_df$gene1), unique(network_df$gene2)))
gene1_notin_gene2
gene2_notin_gene1 <- which(network_df$gene2 %in% setdiff(unique(network_df$gene2), unique(network_df$gene1)))
gene2_notin_gene1


rev_net_1 <- data.frame(
    gene1 = network_df[gene1_notin_gene2, "gene2"],
    gene2 = network_df[gene1_notin_gene2, "gene1"],
    correlation = network_df[gene1_notin_gene2, "correlation"],
    link = rep("dup", length(network_df[gene1_notin_gene2, "correlation"]))
)
rev_net_2 <- data.frame(
    gene1 = network_df[gene2_notin_gene1, "gene2"],
    gene2 = network_df[gene2_notin_gene1, "gene1"],
    correlation = network_df[gene2_notin_gene1, "correlation"],
    link = rep("dup", length(network_df[gene2_notin_gene1, "correlation"]))
)

network_df_full <- rbind(network_df, rev_net_1, rev_net_2)
View(network_df_full)
head(network_df_full)


# compute shortest-path lengths from each gene1 to the target ENSG
target <- "ENSG00000164690"
View(network_df_full)

# build graph from the full edge list (undirected, unweighted)
full_g <- graph_from_data_frame(network_df_full[, c("gene1", "gene2")], directed = FALSE)

# if target not present, assign NA; otherwise compute distances
if (!(target %in% V(full_g)$name)) {
    network_df_full$dist_to_ENSG00000164690 <- NA_integer_
} else {
    dmat <- distances(full_g, v = network_df_full$gene1, to = target, weights = NA)
    dvec <- as.integer(dmat[, 1])
    dvec[is.infinite(dvec)] <- NA_integer_
    network_df_full$dist_to_ENSG00000164690 <- dvec
}


# Folder containing your CSV files
folder <- "article/RWR_patient_results/"

# List all CSV files
csv_files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)

# Read and extract gene lists
gene_lists <- csv_files %>%
    set_names() %>% # keep filenames for tracking
    map(function(file) {
        df <- read_csv(file, show_col_types = FALSE)

        # assuming the gene column is the ONLY column
        # or replace with df$GeneColumnName
        gene_col <- df[[1]] %>% unique()

        data.frame(
            gene = gene_col,
            file = basename(file),
            stringsAsFactors = FALSE
        )
    })

# Combine all CSV contents
all_genes <- bind_rows(gene_lists)

gene_iterations <- all_genes$gene %>%
    table() %>%
    as.data.frame()
rownames(gene_iterations) <- gene_iterations$.
gene_iterations %>% View()

network_df_full$symbol1 <- gene_converter(network_df_full$gene1, "ENSEMBL", "SYMBOL")
network_df_full$symbol2 <- gene_converter(network_df_full$gene2, "ENSEMBL", "SYMBOL")


rwr_occurences <- network_df_full$symbol1 %>%
    mclapply(function(symbol) {
        if (symbol %in% rownames(gene_iterations)) {
            return(gene_iterations[symbol, "Freq"])
        } else {
            return(0)
        }
    }, mc.cores = 10) %>%
    unlist() %>%
    as.vector()
rwr_occurences
network_df_full$rwr_occurences <- rwr_occurences

full_g <- graph_from_data_frame(network_df_full[, c("gene1", "gene2")], directed = FALSE)

### Degree for all genes
deg <- degree(full_g, mode = "all")
deg_df <- tibble(gene1 = names(deg), degree_gene1 = as.integer(deg))

### Mean correlation for each gene1
mean_corr_df <- network_df_full %>%
    group_by(gene1) %>%
    summarise(mean_corr_gene1 = mean(abs(correlation), na.rm = TRUE))

### Merge into network_df_full
network_df_full <- network_df_full %>%
    left_join(deg_df, by = "gene1") %>%
    left_join(mean_corr_df, by = "gene1")

patient_variants[[1]] %>% head()
patient_genes <- lapply(patient_variants, function(df) unique(df$Gene.refGeneWithVer))
names(patient_genes) <- names(patient_variants) # optional but useful

# all unique genes seen in any patient
all_seen_genes <- unique(unlist(patient_genes))

# number of patients
n_patients <- length(patient_genes)

# count in how many patients each gene appears
patient_gene_counts <- table(unlist(
    lapply(patient_genes, unique)
)) %>% as.data.frame()

colnames(patient_gene_counts) <- c("gene", "patients_with_variant")
rownames(patient_gene_counts) <- patient_gene_counts$gene
patient_gene_counts

network_df_full$patient_variant_occurrence <- sapply(
    network_df_full$symbol1,
    function(g) {
        if (g %in% patient_gene_counts$gene) {
            patient_gene_counts[g, "patients_with_variant"]
        } else {
            0
        }
    }
)


View(network_df_full)

write.csv(network_df_full, "article/transit_data/network_observation.csv", row.names = FALSE, quote = FALSE)
network_df_full
