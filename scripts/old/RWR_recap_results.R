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

coex_genes <- rownames(data_3p$coex_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")
HPE_seeds <- lapply(HPE_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
HPE_seeds %>%
    unlist() %>%
    table() %>%
    as.data.frame() %>%
    arrange(desc(Freq)) %>%
    filter(Freq > 2) %>%
    write.csv("article/tables/HPE_seeds_freq.csv", row.names = FALSE)

# set seeds score in patient to 0
HPE_results_noseed <- c(1:ncol(HPE_results)) %>%
    lapply(function(i) {
        x <- HPE_results[, i]
        x[which(rownames(HPE_results) %in% HPE_seeds[[i]])] <- min(HPE_results)
        return(x)
    }) %>%
    Reduce(cbind, .) -> HPE_results_noseed
colnames(HPE_results_noseed) <- colnames(HPE_results)

SHH_patient <- which(sapply(HPE_patient_variants, function(df) {
    "SHH" %in% df$Gene.refGeneWithVer
}))

# get SHH patient top 200 genes in term of mean score
SHH_patient_top200 <- HPE_results_noseed[, SHH_patient] %>%
    rowMeans() %>%
    sort(decreasing = TRUE) %>%
    head(200) %>%
    names()
SHH_patient_top200 %>% write.csv("article/tables/SHH_patient_top200_genes.csv", row.names = FALSE)
View(SHH_patient_top200)

SIX3_patient <- which(sapply(HPE_patient_variants, function(df) {
    "SIX3" %in% df$Gene.refGeneWithVer
}))

# get SIX3 patient top 200 genes in term of mean score
SIX3_patient_top200 <- HPE_results_noseed[, SIX3_patient] %>%
    rowMeans() %>%
    sort(decreasing = TRUE) %>%
    head(200) %>%
    names()
SIX3_patient_top200 %>% write.csv("article/tables/SIX3_patient_top200_genes.csv", row.names = FALSE)
View(SIX3_patient_top200)

ZIC2_patient <- which(sapply(HPE_patient_variants, function(df) {
    "ZIC2" %in% df$Gene.refGeneWithVer
}))

# get ZIC2 patient top 200 genes in term of mean score
ZIC2_patient_top200 <- HPE_results_noseed[, ZIC2_patient] %>%
    rowMeans() %>%
    sort(decreasing = TRUE) %>%
    head(200) %>%
    names()
ZIC2_patient_top200 %>% write.csv("article/tables/ZIC2_patient_top200_genes.csv", row.names = FALSE)
View(ZIC2_patient_top200)

FGF8_patient <- which(sapply(HPE_patient_variants, function(df) {
    "FGF8" %in% df$Gene.refGeneWithVer
}))

# get FGF8 patient top 200 genes in term of mean score
FGF8_patient_top200 <- HPE_results_noseed[, FGF8_patient] %>%
    rowMeans() %>%
    sort(decreasing = TRUE) %>%
    head(200) %>%
    names()
FGF8_patient_top200 %>% write.csv("article/tables/FGF8_patient_top200_genes.csv", row.names = FALSE)
View(FGF8_patient_top200)




#  HPE results mean quantile 0.9 value
mean_quantile <- HPE_results[, 2:199] %>%
    apply(2, function(x) {
        quantile(x, 0.95)
    }) %>%
    mean()
mean_quantile
# max 3086
patient <- 100
HPE_seed_counts[[patient]]
plot(HPE_results[which(!rownames(HPE_results) %in% HPE_seed_counts[[patient]]), patient] %>% sort(decreasing = TRUE))
abline(h = HPE_results[, patient] %>% quantile(0.95), col = "red", lty = 2)
#  value where the patient score plateau define as difference between consecutive scores < quantile 0.95
rownames(HPE_results)[which(HPE_results[which(!rownames(HPE_results) %in% HPE_seed_counts[[patient]]), patient] > HPE_results[, patient] %>% quantile(0.95))]
# get top quantile 0.95 (whithout taking seed genes in account from HPE_seed_counts) for each patient

patient_q95 <- lapply(2:ncol(HPE_results), function(patient) {
    rownames(HPE_results)[which(HPE_results[which(!rownames(HPE_results) %in% HPE_seed_counts[[patient]]), patient] > HPE_results[, patient] %>% quantile(0.95))]
})
patient_q95 %>%
    lapply(length) %>%
    unlist()
patient_q95 %>%
    unlist() %>%
    table() %>%
    as.data.frame() %>%
    View()
HPE_seed_counts <- lapply(HPE_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
# Histogram of HPE seed counts
ggplot(data.frame(counts = HPE_seed_counts), aes(x = counts)) +
    geom_histogram(bins = 20, fill = "#2E86AB", color = "white", linewidth = 0.5) +
    labs(
        title = "Distribution of HPE Seed Counts",
        x = "Number of Seeds per Patient",
        y = "Frequency"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
    )
ggsave("article/figures/HPE_seed_counts_histogram.png", width = 10, height = 6, dpi = 300)

FER_seed_counts <- lapply(FER_patient_variants, function(df) {
    length(filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer)
}) %>% unlist()

# Histogram of FER seed counts
ggplot(data.frame(counts = FER_seed_counts), aes(x = counts)) +
    geom_histogram(bins = 20, fill = "#2E86AB", color = "white", linewidth = 0.5) +
    labs(
        title = "Distribution of IRON Seed Counts",
        x = "Number of Seeds per Patient",
        y = "Frequency"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
    )

ggsave("article/figures/FER_seed_counts_histogram.png", width = 10, height = 6, dpi = 300)



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

# Write top 200 genes by degree to CSV
centrality_matrix %>%
    arrange(desc(degree)) %>%
    slice_head(n = 200) %>%
    write.csv("article/figures/top200_genes_by_degree.csv", row.names = FALSE)
nrow(data_3p$coex_case)

# plot the top 25 genes by degree as a dotplot, with degree on x axis and gene name on y axis, color by degree
top_genes_degree <- centrality_matrix %>%
    arrange(desc(degree)) %>%
    slice_head(n = 50)
ggplot(top_genes_degree, aes(x = degree, y = reorder(node, degree), color = degree)) +
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
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
    )
ggsave("article/figures/top25_genes_by_degree.png", width = 10, height = 6, dpi = 300)

# Histogram of network degree centrality
ggplot(centrality_matrix, aes(x = degree)) +
    geom_histogram(bins = 30, fill = "#A23B72", color = "white", linewidth = 0.5) +
    labs(
        title = "Distribution of Network Degree Centrality",
        x = "Degree",
        y = "Frequency"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
    )

ggsave("article/figures/degree_centrality_histogram.png", width = 10, height = 6, dpi = 300)



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

HPE_long_df <- HPE_score_matrix %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = "patient", values_to = "score")

HPE_long_df <- HPE_long_df %>%
    left_join(
        centrality_matrix %>% dplyr::select(node, degree_rank),
        by = c("gene" = "node")
    )
HPE_long_df$score_scale <- 1 / -log10(HPE_long_df$score)
HPE_long_df_f <- filter(HPE_long_df, score > 0.3)

HPE_gene_order <- HPE_long_df %>%
    group_by(gene) %>%
    summarise(mean_score = mean(score, na.rm = TRUE)) %>%
    arrange(desc(mean_score)) # or arrange(mean_score) for ascending

HPE_long_df <- HPE_long_df %>%
    mutate(gene = factor(gene, levels = HPE_gene_order$gene)) %>%
    arrange(gene)

FER_long_df <- FER_score_matrix %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = "patient", values_to = "score")

FER_long_df <- FER_long_df %>%
    left_join(
        centrality_matrix %>% dplyr::select(node, degree_rank),
        by = c("gene" = "node")
    )
FER_long_df$score_scale <- 1 / -log10(FER_long_df$score)
FER_long_df_f <- filter(FER_long_df, score > 0.3)

FER_gene_order <- FER_long_df %>%
    group_by(gene) %>%
    summarise(mean_score = mean(score, na.rm = TRUE)) %>%
    arrange(desc(mean_score)) # or arrange(mean_score) for ascending

FER_long_df <- FER_long_df %>%
    mutate(gene = factor(gene, levels = FER_gene_order$gene)) %>%
    arrange(gene)

# Reorder HPE_long_df based on that gene order
reorder_HPE_long_df <- HPE_long_df %>%
    mutate(gene = factor(gene, levels = FER_gene_order$gene)) %>%
    arrange(gene)

# Reorder FER_long_df based on that gene order
reorder_FER_long_df <- FER_long_df %>%
    mutate(gene = factor(gene, levels = HPE_gene_order$gene)) %>%
    arrange(gene)

# save tables for gene by gene distribution analyses
save(HPE_long_df, file = "article/transit_data/HPE_RWR_distrib_df.RData")
save(FER_long_df, file = "article/transit_data/FER_RWR_distrib_df.RData")

HPE_profile <- ggplot(data = HPE_long_df, aes(x = gene, y = score, color = degree_rank)) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.5) +
    scale_color_gradientn(
        colours = c("orange", "black", "purple")
    ) +
    # ylim(-0.1, 1) +
    custom_theme() +
    geom_hline(yintercept = 0.02) +
    theme(panel.grid.major.x = element_blank(), axis.text.x = element_blank()) +
    labs(
        title = "RWR gene ranks across HPE patients, ordered by mean score",
        x = "Gene",
        y = "RWR normalized Zscore",
        color = "Degree rank (%)"
    )
ggsave("article/figures/RWR_HPE_profile.png", HPE_profile, width = 12, height = 6)

FER_profile <- ggplot(data = FER_long_df, aes(x = gene, y = score, color = degree_rank)) +
    geom_jitter(width = 0.2, alpha = 0.7, size = 0.5) +
    scale_color_gradientn(
        colours = c("orange", "black", "purple")
    ) +
    custom_theme() +
    theme(panel.grid.major.x = element_blank(), axis.text.x = element_blank()) +
    labs(
        title = "RWR gene ranks across Iron patients, ordered by mean score",
        x = "Gene",
        y = "RWR normalized Zscore",
        color = "Degree rank (%)"
    )
ggsave("article/figures/RWR_FER_profile.png", FER_profile, width = 12, height = 6)


HPE_score_matrix %>%
    rowMeans() %>%
    sort(decreasing = TRUE) %>%
    head(200) %>%
    names() %>%
    write.csv("article/tables/RWR_HPE_top200_genes.csv", row.names = FALSE)

FER_score_matrix %>%
    rowMeans() %>%
    sort(decreasing = TRUE) %>%
    head(200) %>%
    names() %>%
    write.csv("article/tables/RWR_FER_top200_genes.csv", row.names = FALSE)

# combine HPE and FER score matrices

#  add suffixes to patient names to distinguish them
colnames(HPE_score_matrix_common) <- paste0(colnames(HPE_score_matrix_common), "_HPE")
colnames(FER_score_matrix_common) <- paste0(colnames(FER_score_matrix_common), "_FER")

combined_score_matrix <- cbind(HPE_score_matrix_common, FER_score_matrix_common)

# --- 1. Normalize columns to proper probability distributions ---
normalize_columns <- function(mat) {
    apply(mat, 2, function(x) {
        x <- pmax(x, 0) # ensure non-negative
        x / sum(x) # normalize to sum = 1
    })
}

# --- 2. KL divergence (base 2, safe with zeros) ---
kl_divergence <- function(p, q) {
    # Only compute where both > 0
    mask <- p > 0 & q > 0
    sum(p[mask] * log2(p[mask] / q[mask]))
}

# --- 3. Jensen-Shannon Divergence (scalar, for 2 vectors) ---
js_divergence <- function(p, q) {
    m <- (p + q) / 2
    0.5 * kl_divergence(p, m) + 0.5 * kl_divergence(q, m)
}

# --- 4. Jensen-Shannon Distance (square root → true metric) ---
js_distance <- function(p, q) {
    sqrt(js_divergence(p, q))
}

# --- 5. Full pairwise JSD matrix across profiles (columns) ---
jsd_matrix <- function(mat) {
    mat <- normalize_columns(mat)
    n <- ncol(mat)
    nms <- colnames(mat)

    D <- matrix(0, nrow = n, ncol = n, dimnames = list(nms, nms))

    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            d <- js_distance(mat[, i], mat[, j])
            D[i, j] <- d
            D[j, i] <- d
        }
    }
    D
}

rank_matrix <- apply(combined_score_matrix, 2, rank)
D_rank <- jsd_matrix(rank_matrix)
# Réduire la dimensionnalité avant de clusterer
variance_per_gene <- apply(combined_score_matrix, 1, var)

# Visualiser la distribution des variances
hist(log1p(variance_per_gene),
    main = "Distribution des variances par gène",
    xlab = "log(variance + 1)"
)

# Filtrer — seuil à ajuster selon tes données
min_var <- quantile(variance_per_gene, 0.5) # ou plus strict : quantile(variance_per_gene, 0.1)
rwr_filtered <- combined_score_matrix[variance_per_gene > min_var, ]

cat("Gènes avant filtrage :", nrow(combined_score_matrix), "\n")
cat("Gènes après filtrage :", nrow(rwr_filtered), "\n")

# PCA sur matrice filtrée
pca <- prcomp(t(rwr_filtered), scale. = TRUE)

# Visualiser la variance expliquée
plot(pca)

# Clusterer sur les premiers PCs
D_pca <- dist(pca$x[, 1:5]) # ajuster le nombre de PCs
combined_clustering <- hclust(D_pca, method = "ward.D2") # Ward OK sur espace euclidien

combined_clusters <- cutree(combined_clustering, k = 5)
table(combined_clusters)
combined_score_matrix %>% colnames()
combined_clusters

# using complesxHeamap to plot the jaccard similarity matrix with the clusters as annotation, use a color gradient from white to black for the jaccard similarity, and show the patient names on the x and y axis
# Row/col annotation based on clusters
# Cluster order: reorder patients so clusters are contiguous
combined_patient_order <- order(combined_clusters)
combined_cluster_colors <- c("1" = "#E64B35", "2" = "#4DBBD5", "3" = "#00A087", "4" = "#3C5488", "5" = "#F39B7F")

combined_ann <- HeatmapAnnotation(
    Cluster = factor(combined_clusters[combined_patient_order]),
    col = list(Cluster = combined_cluster_colors),
    annotation_legend_param = list(
        Cluster = list(title = "Cluster", nrow = 1)
    ),
    show_annotation_name = TRUE
)
combined_row_ann <- rowAnnotation(
    Cluster = factor(combined_clusters[combined_patient_order]),
    col = list(Cluster = combined_cluster_colors),
    show_annotation_name = FALSE,
    show_legend = FALSE # avoid duplicate legend
)

# Color scale for Jaccard similarity
combined_col_fun <- colorRamp2(
    c(0, 0.1, 0.2, 0.3, 0.4),
    c(
        "white",
        "yellow",
        "orange",
        "purple",
        "black"
    )
)

#  show annotation for FER/HPE status
combined_status <- sapply(colnames(combined_score_matrix), function(name) {
    if (grepl("_HPE$", name)) {
        return("HPE")
    } else if (grepl("_FER$", name)) {
        return("FER")
    } else {
        return(NA)
    }
})
combined_status_colors <- c("HPE" = "#2E86AB", "FER" = "#A23B72")
combined_status_ann <- HeatmapAnnotation(
    Status = factor(combined_status[combined_patient_order]),
    col = list(Status = combined_status_colors),
    annotation_legend_param = list(
        Status = list(title = "Patient Group", nrow = 1)
    ),
    show_annotation_name = TRUE
)

# plot the heatmap
combined_jaccard_HM <- Heatmap(
    combined_jaccard_matrix[combined_patient_order, combined_patient_order],
    name = "Jaccard\nSimilarity",
    col = combined_col_fun,

    # Use pre-defined clustering — disable internal clustering
    cluster_rows = FALSE,
    cluster_columns = FALSE,

    # Annotations
    top_annotation = combined_ann,
    left_annotation = combined_row_ann,
    bottom_annotation = combined_status_ann,

    # Cluster boundary lines
    row_split = factor(combined_clusters[combined_patient_order]),
    column_split = factor(combined_clusters[combined_patient_order]),
    row_gap = unit(2, "mm"),
    column_gap = unit(2, "mm"),

    # Labels
    row_title = "Patients",
    column_title = "HPE and FER Patients — Jaccard Similarity of Top 200 RWR Genes",
    column_title_side = "top",
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),

    # Hide individual patient labels if many patients
    show_row_names = FALSE,
    show_column_names = FALSE,
    # Cell appearance
    rect_gp = gpar(col = NA), # no cell borders (cleaner for large matrices)

    heatmap_legend_param = list(
        direction = "vertical",
        legend_height = unit(4, "cm")
    )
)
combined_jaccard_HM




# jaccard distance of top 200 genes of each patient
HPE_top_genes_list <- lapply(1:ncol(HPE_score_matrix), function(i) {
    patient_scores <- HPE_score_matrix[, i]
    top_genes <- names(sort(patient_scores, decreasing = TRUE))[1:200]
    return(top_genes)
})
names(HPE_top_genes_list) <- colnames(HPE_score_matrix)

HPE_jaccard_matrix <- matrix(0, nrow = ncol(HPE_score_matrix), ncol = ncol(HPE_score_matrix))
for (i in 1:ncol(HPE_score_matrix)) {
    for (j in 1:ncol(HPE_score_matrix)) {
        HPE_jaccard_matrix[i, j] <- length(intersect(HPE_top_genes_list[[i]], HPE_top_genes_list[[j]])) / length(union(HPE_top_genes_list[[i]], HPE_top_genes_list[[j]]))
    }
}

as.dist(1 - HPE_jaccard_matrix) %>%
    hclust() -> clustering

clusters <- cutree(clustering, k = 4)
table(clusters)

# get the 200 top mean genes for each cluster
cluster_top_genes <- lapply(1:4, function(cluster) {
    cluster_patients <- which(clusters == cluster)
    cluster_mean_scores <- rowMeans(HPE_score_matrix[, cluster_patients, drop = FALSE])
    top_genes <- names(sort(cluster_mean_scores, decreasing = TRUE))[1:200]
    return(top_genes)
})

# print top 10 genes for each cluster
for (i in 1:4) {
    cat(paste0("Cluster ", i, " top 10 genes:\n"))
    print(head(cluster_top_genes[[i]], 10))
    cat("\n")
}

# write them in files
for (i in 1:4) {
    write.csv(cluster_top_genes[[i]], paste0("article/tables/RWR_HPE_cluster_", i, "_top200_genes.csv"), row.names = FALSE)
}

# plot HPE jaccard matrix as heatmap, use previously define clustering use the complexHeatmap package, show annotation for clusters
# Row/col annotation based on clusters
# Cluster order: reorder patients so clusters are contiguous
patient_order <- order(clusters)
patient_order
clusters

cluster_colors <- c("1" = "#E64B35", "2" = "#4DBBD5", "3" = "#00A087", "4" = "#3C5488")

ann <- HeatmapAnnotation(
    Cluster = factor(clusters[patient_order]),
    col = list(Cluster = cluster_colors),
    annotation_legend_param = list(
        Cluster = list(title = "Cluster", nrow = 1)
    ),
    show_annotation_name = TRUE
)
row_ann <- rowAnnotation(
    Cluster = factor(clusters[patient_order]),
    col = list(Cluster = cluster_colors),
    show_annotation_name = FALSE,
    show_legend = FALSE # avoid duplicate legend
)


# Color scale for Jaccard similarity
col_fun <- colorRamp2(
    c(0, 0.1, 0.2, 0.3, 0.4),
    c(
        "white",
        "yellow",
        "orange",
        "purple",
        "black"
    )
)


jaccard_HM <- Heatmap(
    HPE_jaccard_matrix[patient_order, patient_order],
    name = "Jaccard\nSimilarity",
    col = col_fun,

    # Use pre-defined clustering — disable internal clustering
    cluster_rows = FALSE,
    cluster_columns = FALSE,

    # Annotations
    top_annotation = ann,
    left_annotation = row_ann,

    # Cluster boundary lines
    row_split = factor(clusters[patient_order]),
    column_split = factor(clusters[patient_order]),
    row_gap = unit(2, "mm"),
    column_gap = unit(2, "mm"),

    # Labels
    row_title = "Patients",
    column_title = "HPE Patients — Jaccard Similarity of Top 200 RWR Genes",
    column_title_side = "top",
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),

    # Hide individual patient labels if many patients
    show_row_names = FALSE,
    show_column_names = FALSE,
    # Cell appearance
    rect_gp = gpar(col = NA), # no cell borders (cleaner for large matrices)

    heatmap_legend_param = list(
        direction = "vertical",
        legend_height = unit(4, "cm")
    )
)

png_save(
    plot = jaccard_HM,
    path = "article/figures/RWR_HPE_jaccard_heatmap.png",
    width = 2550,
    height = 2000
)
clusters
# add cluster inf to clinic_f data

names(clusters) <- colnames(HPE_score_matrix)

load("article/transit_data/clinic_f.RData")
clinic_f$Cluster <- factor(clusters[clinic_f$index])
View(clinic_f)
save(clinic_f, file = "article/transit_data/clinic_f_cluster.RData")

HPE_seeds <- lapply(HPE_patient_variants, function(patient_df) {
    patient_genes <- intersect(rownames(HPE_score_matrix), patient_df$Gene.refGeneWithVer)
    if (length(patient_genes) == 0) {
        return(NA)
    }
    return(patient_genes)
})
names(HPE_seeds) <- paste0("Patient_", seq_along(HPE_seeds))


# rowMeans(HPE_score_matrix) %>%
#     sort(decreasing = TRUE) %>%
#     head(200) %>%
#     names() %>%
#     write.csv("article/tables/RWR_HPE_top200_genes.csv", row.names = FALSE)

# rowMeans(FER_score_matrix) %>%
#     sort(decreasing = TRUE) %>%
#     head(200) %>%
#     names() %>%
#     write.csv("article/tables/RWR_FER_top200_genes.csv", row.names = FALSE)
FER_score_matrix
FER_score_matrix_red <- 1 / -log10(FER_score_matrix)
FER_score_matrix_red
FER_score_matrix_red[which(FER_score_matrix_red <= 0.25)] <- 0

names(FER_patient_variants) <- paste0("Patient_", seq_along(FER_patient_variants))

colnames(FER_score_matrix) %>%
    lapply(function(patient) {
        patient_variant_genes <- unique(FER_patient_variants[[patient]]$Gene.refGeneWithVer)
        patient_seeds <- intersect(patient_variant_genes, genes)
        patient_seed_ratio <- length(patient_seeds) / length(patient_variant_genes)
        patient_top50_genes <- names(which(FER_score_matrix[, patient] >= sort(FER_score_matrix[, patient], decreasing = TRUE)[50]))
        patient_recall50_genes <- intersect(patient_top50_genes, patient_seeds)
        patient_top200_genes <- names(which(FER_score_matrix[, patient] >= sort(FER_score_matrix[, patient], decreasing = TRUE)[200]))
        patient_recall200_genes <- intersect(patient_top200_genes, patient_seeds)
        patient_top500_genes <- names(which(FER_score_matrix[, patient] >= sort(FER_score_matrix[, patient], decreasing = TRUE)[500]))
        patient_recall500_genes <- intersect(patient_top500_genes, patient_seeds)
        return(c(
            patient,
            length(patient_variant_genes),
            length(patient_seeds),
            patient_seed_ratio,
            length(patient_recall50_genes),
            length(patient_recall200_genes),
            length(patient_recall500_genes),
            paste(patient_recall50_genes, collapse = "/"),
            paste(patient_recall200_genes, collapse = "/"),
            paste(patient_recall500_genes, collapse = "/"),
            paste(patient_variant_genes, collapse = "/"),
            paste(patient_seeds, collapse = "/"),
            paste(patient_top50_genes, collapse = "/"),
            paste(patient_top200_genes, collapse = "/"),
            paste(patient_top500_genes, collapse = "/")
        ))
    }) %>%
    Reduce(rbind, .) %>%
    as.data.frame() -> FER_recap
colnames(FER_recap) <- c(
    "Patient",
    "Variant_gene_count",
    "Seed_gene_count",
    "Seed_ratio",
    "Recall50_count",
    "Recall200_count",
    "Recall500_count",
    "Recall50_genes",
    "Recall200_genes",
    "Recall500_genes",
    "Variant_genes",
    "Seed_genes",
    "Top50_genes",
    "Top200_genes",
    "Top500_genes"
)
FER_recap$Seed_ratio <- as.numeric(FER_recap$Seed_ratio)
View(FER_recap)
write.csv(FER_recap, "article/tables/RWR_FER_recap.csv")

names(HPE_patient_variants) <- paste0("Patient_", seq_along(HPE_patient_variants))
HPE_patient_variants
names(HPE_patient_variants)

colnames(HPE_score_matrix) %>%
    lapply(function(patient) {
        patient_variant_genes <- unique(HPE_patient_variants[[patient]]$Gene.refGeneWithVer)
        patient_seeds <- intersect(patient_variant_genes, genes)
        patient_seed_ratio <- length(patient_seeds) / length(patient_variant_genes)
        patient_top50_genes <- names(which(HPE_score_matrix[, patient] >= sort(HPE_score_matrix[, patient], decreasing = TRUE)[50]))
        patient_recall50_genes <- intersect(patient_top50_genes, patient_seeds)
        patient_top200_genes <- names(which(HPE_score_matrix[, patient] >= sort(HPE_score_matrix[, patient], decreasing = TRUE)[200]))
        patient_recall200_genes <- intersect(patient_top200_genes, patient_seeds)
        patient_top500_genes <- names(which(HPE_score_matrix[, patient] >= sort(HPE_score_matrix[, patient], decreasing = TRUE)[500]))
        patient_recall500_genes <- intersect(patient_top500_genes, patient_seeds)
        return(c(
            patient,
            as.numeric(length(patient_variant_genes)),
            as.numeric(length(patient_seeds)),
            patient_seed_ratio,
            as.numeric(length(patient_recall50_genes)),
            as.numeric(length(patient_recall200_genes)),
            as.numeric(length(patient_recall500_genes)),
            paste(patient_recall50_genes, collapse = "/"),
            paste(patient_recall200_genes, collapse = "/"),
            paste(patient_recall500_genes, collapse = "/"),
            paste(patient_variant_genes, collapse = "/"),
            paste(patient_seeds, collapse = "/"),
            paste(patient_top50_genes, collapse = "/"),
            paste(patient_top200_genes, collapse = "/"),
            paste(patient_top500_genes, collapse = "/")
        ))
    }) %>%
    Reduce(rbind, .) %>%
    as.data.frame() -> HPE_recap
colnames(HPE_recap) <- c(
    "Patient",
    "Variant_gene_count",
    "Seed_gene_count",
    "Seed_ratio",
    "Recall50_count",
    "Recall200_count",
    "Recall500_count",
    "Recall50_genes",
    "Recall200_genes",
    "Recall500_genes",
    "Variant_genes",
    "Seed_genes",
    "Top50_genes",
    "Top200_genes",
    "Top500_genes"
)
HPE_recap$Seed_ratio <- as.numeric(HPE_recap$Seed_ratio)
HPE_recap$Variant_gene_count <- as.numeric(HPE_recap$Variant_gene_count)
HPE_recap$Seed_gene_count <- as.numeric(HPE_recap$Seed_gene_count)
HPE_recap$Recall50_count <- as.numeric(HPE_recap$Recall50_count)
HPE_recap$Recall200_count <- as.numeric(HPE_recap$Recall200_count)
HPE_recap$Recall500_count <- as.numeric(HPE_recap$Recall500_count)
View(HPE_recap)
write.csv(HPE_recap, "article/tables/RWR_HPE_recap.csv")

# # Seed ratio distribution
# seed_ratio_df <- rbind(
#     HPE_recap %>% dplyr::select(Patient, Seed_ratio) %>% mutate(Group = "HPE"),
#     FER_recap %>% dplyr::select(Patient, Seed_ratio) %>% mutate(Group = "FER")
# )
# seed_ratio_plot <- ggplot(data = seed_ratio_df, aes(x = Group, y = Seed_ratio, color = Group)) +
#     geom_boxplot(alpha = 0.5, outlier.shape = NA) +
#     geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
#     custom_theme() +
#     ylim(0, 1) +
#     labs(
#         title = "Seed ratio distribution across patient groups",
#         x = "Group",
#         y = "Seed ratio",
#         color = "Group"
#     ) +
#     theme(legend.position = "none")

# seed_ratio_plot
# ggsave("article/figures/RWR_seed_ratio_distribution.png", seed_ratio_plot, width = 10, height = 6)
# data_3p$coex_case %>% nrow()
