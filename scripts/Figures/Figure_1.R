# Loading packages and functions
library(Matrix)
library(tidyverse)
library(parallel)
library(tibble)
library(igraph)
library(ggplot2)
library(stats)

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
source("article/scripts/RECENT_webgestalt_ORA.R")

load("article/transit_data/data_3p.RData")

g <- graph_from_adjacency_matrix(
    data_3p$coex_case,
    mode = "undirected", # or "directed", "upper", "lower", "max", "min", "plus"
    weighted = TRUE, # set to TRUE to use matrix values as edge weights
    diag = FALSE # FALSE to ignore self-loops on the diagonal
)

centrality_matrix <- data.frame(
    node = V(g)$name %>% gene_converter(from = "ENSEMBL", to = "SYMBOL"),
    degree = igraph::degree(g),
    degree_rank = percent_rank(igraph::degree(g)) * 100,
    betweenness = betweenness(g),
    closeness = closeness(g),
    eigenvector = eigen_centrality(g)$vector
)


# Figure 1 : iPSC anterior neuroectoderm transcriptomic data enable the construction of a highly contextual co-expression network
# Figure 1A : network degree distribution
# Histogram of network degree centrality
degree_centrality_histogram <- ggplot(centrality_matrix, aes(x = degree)) +
    geom_histogram(bins = 30, fill = "#A23B72", color = "white", linewidth = 0.5) +
    labs(
        title = "Distribution of Network Degree Centrality",
        x = "Degree",
        y = "Frequency"
    ) +
    custom_theme() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
    )

ggsave("article/figures/Figure_1/degree_centrality_histogram.png", plot = degree_centrality_histogram, width = 10, height = 6, dpi = 300)


# Figure 1:B : top 50 genes by degree centrality
top_genes_degree <- centrality_matrix %>%
    arrange(desc(degree)) %>%
    slice_head(n = 50)
top50_genes_by_degree <- ggplot(top_genes_degree, aes(x = degree, y = reorder(node, degree), color = degree)) +
    geom_point(size = 3) +
    scale_color_gradientn(
        colours = c("orange", "black", "purple")
    ) +
    labs(
        title = "Top 50 Genes by Degree Centrality",
        x = "Degree",
        y = "Gene",
        color = "Degree"
    ) +
    custom_theme() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
    )
ggsave("article/figures/Figure_1/top50_genes_by_degree.png", plot = top50_genes_by_degree, width = 10, height = 6, dpi = 300)

# Figure 1:C: Enrichment analysis of the top 500 genes by degree centrality
top_genes_500 <- centrality_matrix %>%
    arrange(desc(degree)) %>%
    slice_head(n = 500) %>%
    pull(node)

GO_enrichment_results <- run_ora_gobp(
    genelist = top_genes_500,
    interestGeneType = "genesymbol",
    fdrThr = 0.05,
    minNum = 5,
    maxNum = 2000,
    noRedundant = TRUE
)
GO_plot <- plot_ora_dotplot(
    ora_results = GO_enrichment_results,
    n = 20
)
ggsave("article/figures/Figure_1/GO_BP_enrichment_top500_network.png", plot = GO_plot, width = 10, height = 6, dpi = 300)

HPO_enrichment_results <- run_ora_hpo(
    genelist = top_genes_500,
    interestGeneType = "genesymbol",
    fdrThr = 0.05,
    minNum = 5,
    maxNum = 2000
)
HPO_plot <- plot_ora_dotplot(
    ora_results = HPO_enrichment_results,
    n = 20
)
ggsave("article/figures/Figure_1/HPO_enrichment_top500_network.png", plot = HPO_plot, width = 10, height = 6, dpi = 300)

reactome_enrichment_results <- run_ora_reactome(
    genelist = top_genes_500,
    interestGeneType = "genesymbol",
    fdrThr = 0.05,
    minNum = 5,
    maxNum = 2000
)
reactome_plot <- plot_ora_dotplot(
    ora_results = reactome_enrichment_results,
    n = 20
)
ggsave("article/figures/Figure_1/reactome_enrichment_top500_network.png", plot = reactome_plot, width = 10, height = 6, dpi = 300)

# ---- Figure 1 summary dataframes ----

gene_network_summary <- centrality_matrix %>%
    mutate(
        top50 = node %in% top_genes_degree$node,
        top500 = node %in% top_genes_500
    ) %>%
    arrange(desc(degree))

GO_enrichment_summary <- GO_enrichment_results %>%
    dplyr::select(geneSet, description, size, overlap, enrichmentRatio, pValue, FDR, userId) %>%
    arrange(-enrichmentRatio)

HPO_enrichment_summary <- HPO_enrichment_results %>%
    dplyr::select(geneSet, description, size, overlap, enrichmentRatio, pValue, FDR, userId) %>%
    arrange(-enrichmentRatio)

reactome_enrichment_summary <- reactome_enrichment_results %>%
    dplyr::select(geneSet, description, size, overlap, enrichmentRatio, pValue, FDR, userId) %>%
    arrange(-enrichmentRatio)

write_csv(gene_network_summary, "article/tables/Figure_1/gene_network_summary.csv")
write_csv(GO_enrichment_summary, "article/tables/Figure_1/GO_BP_enrichment_summary.csv")
write_csv(HPO_enrichment_summary, "article/tables/Figure_1/HPO_enrichment_summary.csv")
write_csv(reactome_enrichment_summary, "article/tables/Figure_1/reactome_enrichment_summary.csv")

library(openxlsx)

write_single_xlsx <- function(df, path, sheet_name = "Sheet1") {
    wb <- createWorkbook()
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df, headerStyle = createStyle(textDecoration = "bold"))
    freezePane(wb, sheet_name, firstRow = TRUE)
    setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
    saveWorkbook(wb, path, overwrite = TRUE)
}

write_single_xlsx(gene_network_summary, "article/tables/Figure_1/gene_network_summary.xlsx", "gene_network_summary")
write_single_xlsx(GO_enrichment_summary, "article/tables/Figure_1/GO_BP_enrichment_summary.xlsx", "GO_BP_enrichment")
write_single_xlsx(HPO_enrichment_summary, "article/tables/Figure_1/HPO_enrichment_summary.xlsx", "HPO_enrichment")
write_single_xlsx(reactome_enrichment_summary, "article/tables/Figure_1/reactome_enrichment_summary.xlsx", "reactome_enrichment")
