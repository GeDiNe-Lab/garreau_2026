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


load("article/transit_data/clinic_f_cluster.RData")
clinic_f %>% View()

HPE_top_genes_list


## ============================================================
##  Gene Enrichment Analysis: Gene Ontology + Human Phenotype Ontology
##  Dependencies: clusterProfiler, org.Hs.eg.db, enrichplot,
##                phenOnto / hpoglue  (HPO via ontologyIndex + DOSE)
## ============================================================

# ── 0. Install / load packages ────────────────────────────────────────────────
packages_bioc <- c(
    "clusterProfiler", "org.Hs.eg.db", "enrichplot",
    "AnnotationDbi", "DOSE", "ReactomePA"
)
packages_cran <- c(
    "ontologyIndex", "ontologySimilarity",
    "ggplot2", "dplyr", "stringr", "patchwork"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

for (pkg in packages_bioc) {
    if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}
for (pkg in packages_cran) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(AnnotationDbi)
    library(DOSE)
    library(ggplot2)
    library(dplyr)
    library(stringr)
    library(patchwork)
    library(ontologyIndex)
})


# ── 1. Your gene list ────────────────────────────────────────────────────────
#  Replace with your own genes (HGNC symbols, Entrez IDs, or Ensembl IDs).
#  Here we use an example list of 40 genes involved in DNA damage / cancer.

clinic_f
# based on cluster and patient columns, we can get the top genes 200 genes for each cluster
#  using mean score


HPE_top_genes_cluster <- lapply(c(1:4), function(i) {
    patient_scores <- HPE_score_matrix[, clinic_f$index[which(clinic_f$Cluster == i)]] %>%
        rowMeans(na.rm = TRUE)
    top_genes <- names(sort(patient_scores, decreasing = TRUE))[1:200]
    return(top_genes)
})
HPE_top_genes_cluster

# Background / universe: all genes in the annotation DB (recommended)
all_entrez <- keys(org.Hs.eg.db, keytype = "ENTREZID")


# ── 2. ID conversion  (Symbol → Entrez) ──────────────────────────────────────
id_map <- bitr(gene_symbols,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
)

cat(sprintf(
    "Input: %d symbols  |  Mapped: %d Entrez IDs\n",
    length(gene_symbols), nrow(id_map)
))

entrez_ids <- id_map$ENTREZID


# ── 3. Gene Ontology enrichment ───────────────────────────────────────────────
run_go <- function(genes, ont = "BP", universe = all_entrez) {
    enrichGO(
        gene          = genes,
        universe      = universe,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = ont, # "BP", "MF", or "CC"
        pAdjustMethod = "fdr",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05,
        readable      = TRUE # convert Entrez → Symbol in results
    )
}

go_bp <- run_go(entrez_ids, ont = "BP")

go_mf <- run_go(entrez_ids, ont = "MF")
go_cc <- run_go(entrez_ids, ont = "CC")

cat(sprintf("\nGO-BP significant terms : %d\n", nrow(go_bp)))
cat(sprintf("GO-MF significant terms : %d\n", nrow(go_mf)))
cat(sprintf("GO-CC significant terms : %d\n", nrow(go_cc)))


# ── 4. HPO enrichment ─────────────────────────────────────────────────────────
#  Strategy: download the official HPO gene-phenotype annotations from JAX,
#  build a term2gene / term2name data.frame, then use enricher().

hpo_annot_url <- paste0(
    "https://github.com/obophenotype/human-phenotype-ontology/",
    "releases/latest/download/genes_to_phenotype.txt"
)

hpo_obo_url <- paste0(
    "https://github.com/obophenotype/human-phenotype-ontology/",
    "releases/latest/download/hp.obo"
)

hpo_annot_file <- "genes_to_phenotype.txt"
hpo_obo_file <- "hp.obo"

# Download once, cache locally
if (!file.exists(hpo_annot_file)) {
    message("Downloading HPO gene-phenotype annotations …")
    download.file(hpo_annot_url, hpo_annot_file, quiet = TRUE)
}
if (!file.exists(hpo_obo_file)) {
    message("Downloading hp.obo …")
    download.file(hpo_obo_url, hpo_obo_file, quiet = TRUE)
}

# Parse gene-phenotype table
#  Columns (no header): ncbi_gene_id | gene_symbol | hpo_id | hpo_name | ...
hpo_raw <- read.delim(hpo_annot_file,
    comment.char = "#", header = FALSE,
    col.names = c(
        "ncbi_gene_id", "gene_symbol",
        "hpo_id", "hpo_name",
        "frequency", "disease_id"
    ),
    fill = TRUE, quote = ""
)

# Build TERM2GENE and TERM2NAME for enricher()
term2gene <- hpo_raw %>%
    select(hpo_id, ncbi_gene_id) %>%
    mutate(ncbi_gene_id = as.character(ncbi_gene_id)) %>%
    distinct()

term2name <- hpo_raw %>%
    select(hpo_id, hpo_name) %>%
    distinct()

# HPO universe = all genes in the annotation file
hpo_universe <- unique(as.character(hpo_raw$ncbi_gene_id))

# Run enrichment
hpo_res <- enricher(
    gene          = entrez_ids,
    universe      = hpo_universe,
    TERM2GENE     = term2gene,
    TERM2NAME     = term2name,
    pAdjustMethod = "FDR",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
)

cat(sprintf("HPO significant terms    : %d\n", nrow(hpo_res)))


enrichment <- HPE_top_genes_cluster %>% lapply(function(cluster) {
    gene_symbols <- cluster
    entrez_ids <- bitr(gene_symbols,
        fromType = "SYMBOL",
        toType   = "ENTREZID",
        OrgDb    = org.Hs.eg.db
    )$ENTREZID
    go_bp <- run_go(entrez_ids, ont = "BP")
    hpo_res <- enricher(
        gene          = entrez_ids,
        universe      = hpo_universe,
        TERM2GENE     = term2gene,
        TERM2NAME     = term2name,
        pAdjustMethod = "fdr",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.05
    )
    return(list(
        hpo_res = hpo_res,
        go_bp = go_bp
    ))
})
names(enrichment) <- paste0("Cluster_", 1:4)
enrichment$Cluster_2$hpo_res %>%
    as.data.frame() %>%
    View()

# ── 5. Results tables ─────────────────────────────────────────────────────────
print_top <- function(res, label, n = 15) {
    cat(sprintf("\n── Top %d %s terms ──────────────────────────────────\n", n, label))
    df <- as.data.frame(res)
    if (nrow(df) == 0) {
        cat("  (no significant terms)\n")
        return(invisible(NULL))
    }
    df %>%
        arrange(p.adjust) %>%
        head(n) %>%
        select(ID, Description, GeneRatio, BgRatio, p.adjust, Count) %>%
        print(row.names = FALSE)
}

print_top(go_bp, "GO-BP")
print_top(hpo_res, "HPO")













# ── 6. Visualisations ─────────────────────────────────────────────────────────
# helper: safe dotplot wrapper
safe_dotplot <- function(res, title, showCategory = 15) {
    if (is.null(res) || nrow(as.data.frame(res)) == 0) {
        return(NULL)
    }
    dotplot(res, showCategory = showCategory, font.size = 9) +
        ggtitle(title) +
        theme(plot.title = element_text(face = "bold", size = 11))
}

p_gobp <- safe_dotplot(go_bp, "GO – Biological Process")
p_gomf <- safe_dotplot(go_mf, "GO – Molecular Function")
p_gocc <- safe_dotplot(go_cc, "GO – Cellular Component")
p_hpo <- safe_dotplot(hpo_res, "HPO – Human Phenotype Ontology")

## 6a. Combined dotplots (GO panels)
pdf("go_dotplots.pdf", width = 16, height = 7)
print(p_gobp + p_gomf + p_gocc)
dev.off()
message("Saved: go_dotplots.pdf")

## 6b. HPO dotplot
if (!is.null(p_hpo)) {
    pdf("hpo_dotplot.pdf", width = 9, height = 7)
    print(p_hpo)
    dev.off()
    message("Saved: hpo_dotplot.pdf")
}

## 6c. Enrichment map (GO-BP) — clusters related terms
if (nrow(as.data.frame(go_bp)) >= 5) {
    go_bp2 <- pairwise_termsim(go_bp)
    pdf("go_bp_emapplot.pdf", width = 10, height = 8)
    print(emapplot(go_bp2, showCategory = 30, cex_label_category = 0.7) +
        ggtitle("GO-BP enrichment map"))
    dev.off()
    message("Saved: go_bp_emapplot.pdf")
}

## 6d. Upset plot — gene overlap across GO-BP terms
if (nrow(as.data.frame(go_bp)) >= 3) {
    pdf("go_bp_upsetplot.pdf", width = 10, height = 5)
    print(upsetplot(go_bp, n = 12) + ggtitle("GO-BP gene overlap (upset)"))
    dev.off()
    message("Saved: go_bp_upsetplot.pdf")
}

## 6e. Bar chart for HPO
if (!is.null(hpo_res) && nrow(as.data.frame(hpo_res)) > 0) {
    hpo_df <- as.data.frame(hpo_res) %>%
        arrange(p.adjust) %>%
        head(20) %>%
        mutate(
            Description = str_trunc(Description, 55),
            Description = factor(Description, levels = rev(Description))
        )

    p_hpo_bar <- ggplot(hpo_df, aes(x = Count, y = Description, fill = p.adjust)) +
        geom_col() +
        scale_fill_gradient(
            low = "#e63946", high = "#a8dadc",
            name = "adj. p-value"
        ) +
        labs(title = "Top HPO phenotype terms", x = "Gene count", y = NULL) +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(face = "bold"))

    pdf("hpo_barplot.pdf", width = 9, height = 7)
    print(p_hpo_bar)
    dev.off()
    message("Saved: hpo_barplot.pdf")
}


# ── 7. Export results to CSV ──────────────────────────────────────────────────
write_res <- function(res, fname) {
    df <- as.data.frame(res)
    if (nrow(df) > 0) {
        write.csv(df, fname, row.names = FALSE)
        message(sprintf("Saved: %s  (%d terms)", fname, nrow(df)))
    }
}

write_res(go_bp, "results_GO_BP.csv")
write_res(go_mf, "results_GO_MF.csv")
write_res(go_cc, "results_GO_CC.csv")
write_res(hpo_res, "results_HPO.csv")


# ── 8. Session info ───────────────────────────────────────────────────────────
cat("\n── Session info ─────────────────────────────────────────\n")
sessionInfo()
