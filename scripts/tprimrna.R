# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(tibble)

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

rawcounts <- readcounts("panasenkava_2024/data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("panasenkava_2024/data/meta.csv", sep = ",", header = TRUE)

meta_case <- filter(rawmeta, CRISPR %in% c("no", "control"), type %in% c("ventral", "cyclo"), sample != "L9C1_2", diff %in% c("diff9", "diff12"))
rawcounts_case <- rawcounts[, meta_case$sample]
counts_case <- rawcounts_case[which(rowSums(rawcounts_case) >= 25), ]

meta_control <- filter(rawmeta, CRISPR %in% c("no", "control"), type == "ventral", sample != "L9C1_2", day %in% c("day10", "day12"))
rawcounts_control <- rawcounts[, meta_control$sample]
counts_control <- rawcounts_control[which(rowSums(rawcounts_control) >= 25), ]

load("article/transit_data/sc_genes.RData")

norm_case <- DESeqDataSetFromMatrix(
    countData = counts_case[which(gene_converter(rownames(counts_case), "ENSEMBL", "SYMBOL") %in% sc_genes == TRUE), ],
    colData = meta_case,
    design = ~ line * cyclo_dose_qual
) %>%
    vst(blind = FALSE) %>%
    assay() %>%
    limma::removeBatchEffect(batch = meta_case$line)

DEGs_case <- DESeqDataSetFromMatrix(
    countData = counts_case[which(gene_converter(rownames(counts_case), "ENSEMBL", "SYMBOL") %in% sc_genes == TRUE), ],
    colData = meta_case,
    design = ~ line * cyclo_dose_qual
) %>%
    DESeq() %>%
    results(alpha = 0.05, contrast = c("cyclo_dose_qual", "high", "none")) %>%
    as.data.frame() %>%
    na.omit() %>%
    filter(padj < 0.05)
DEGs_case$gene <- rownames(DEGs_case) %>% gene_converter("ENSEMBL", "SYMBOL")

norm_control() <- DESeqDataSetFromMatrix(
    countData = counts_control[which(gene_converter(rownames(counts_control), "ENSEMBL", "SYMBOL") %in% sc_genes == TRUE), ],
    colData = meta_control,
    design = ~ line * day
) %>%
    vst(blind = FALSE) %>%
    assay() %>%
    limma::removeBatchEffect(batch = meta_control$line)


sft_case <- pickSoftThreshold(t(norm_case),
    powerVector = c(1:20),
    verbose = 5
)

sft_control <- pickSoftThreshold(t(norm_control),
    powerVector = c(1:20),
    verbose = 5
)

filter(sft_case$fitIndices, SFT.R.sq > 0.8, mean.k. < 100)
filter(sft_control$fitIndices, SFT.R.sq > 0.8, mean.k. < 100)

coex_case <- cor(t(norm_case))**4
coex_control <- cor(t(norm_control))**6

raw_case_net <- abs(coex_case) %>% as.matrix()
raw_case_net[abs(raw_case_net) < 0.1] <- 0
case_keep <- which(rowSums(raw_case_net != 0) > 1)
case_net <- raw_case_net[case_keep, case_keep]

nrow(case_net)

raw_control_net <- abs(coex_control) %>% as.matrix()
raw_control_net[abs(raw_control_net) < 0.1] <- 0
control_keep <- which(rowSums(raw_control_net != 0) > 1)
control_net <- raw_control_net[control_keep, control_keep]

rownames(norm_case) %>%
    gene_converter("ENSEMBL", "SYMBOL") %>%
    write.csv("article/tables/context_genes.txt", row.names = FALSE, quote = FALSE)

data_3p <- list(
    norm_case = norm_case,
    norm_control = norm_control,
    meta_case = meta_case,
    meta_control = meta_control,
    coex_case = case_net,
    coex_control = control_net,
    DEGs_case = DEGs_case
)
save(data_3p, file = "article/transit_data/data_3p.RData")
View(data_3p$DEGs_case)
case_net[which(case_net != 0)] %>%
    as.vector() %>%
    hist(breaks = 100)













data_3p$meta_LON %>% View()
data_3p$meta_LON[which(data_3p$meta_LON$type == "dorsal"), ]$cyclo_dose_qual <- "dorsal"
data_3p$meta_WTC[which(data_3p$meta_WTC$type == "dorsal"), ]$cyclo_dose_qual <- "dorsal"

genes_plot <- c("ENSG00000133124", "ENSG00000139174", "ENSG00000142627", "ENSG00000007372")

# COMBINING BOTH CELL LINES
# combine LON and WTC for the same genes and plot together
lon_mat <- data_3p$norm_LON[genes_plot, , drop = FALSE]
rownames(lon_mat) <- gene_converter(rownames(lon_mat), "ENSEMBL", "SYMBOL")
lon_df <- as.data.frame(t(lon_mat)) %>%
    rownames_to_column("sample") %>%
    pivot_longer(cols = -sample, names_to = "gene", values_to = "count") %>%
    left_join(data_3p$meta_LON, by = "sample") %>%
    mutate(cell_line = "LON")

wtc_mat <- data_3p$norm_WTC[genes_plot, , drop = FALSE]
rownames(wtc_mat) <- gene_converter(rownames(wtc_mat), "ENSEMBL", "SYMBOL")
wtc_df <- as.data.frame(t(wtc_mat)) %>%
    rownames_to_column("sample") %>%
    pivot_longer(cols = -sample, names_to = "gene", values_to = "count") %>%
    left_join(data_3p$meta_WTC, by = "sample") %>%
    mutate(cell_line = "WTC")

combined_df <- bind_rows(lon_df, wtc_df)
plot <- ggplot(combined_df, aes(x = cyclo_dose_qual, y = count, fill = cyclo_dose_qual)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(shape = cell_line), width = 0.2, size = 3, alpha = 0.7) +
    facet_wrap(~gene, scales = "free_y") +
    scale_fill_manual(values = c("#A1A1DE", "#cf6509", "#eaac30", "#80AD3C")) +
    theme_minimal() +
    labs(x = "cyclo_dose_qual", y = "Count (assay)", title = "Merged expression (both cell lines)", shape = "Cell line") +
    guides(fill = "none") +
    theme(legend.position = "right")

print(plot)
ggsave("article/transit_data/3p_combined_expression.png", plot, width = 10, height = 6)

# heatmap (ComplexHeatmap) for the previously plotted genes across both cell lines

genes_plot <- c(consensus_positions_with_stats$Gene, "SHH")
genes_plot
genes_plot <- genes_plot[genes_plot %in% (rownames(data_3p$norm_LON) %>%
    gene_converter("SYMBOL", "ENSEMBL") %>%
    unique())]

library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)

# prepare expression matrix (genes x samples) and sample metadata combined
lon_mat <- data_3p$norm_LON[genes_plot, , drop = FALSE]
rownames(lon_mat) <- gene_converter(rownames(lon_mat), "ENSEMBL", "SYMBOL")
nrow(lon_mat)

# prepare sample ordering and column annotation for LON-only heatmap
meta_lon <- data_3p$meta_LON %>% filter(sample %in% colnames(lon_mat))

# prefer this biological ordering if present, otherwise keep observed order
preferred_order <- c("ventral", "cyclo", "dorsal")
present_pref <- intersect(preferred_order, unique(meta_lon$cyclo_dose_qual))
meta_lon$cyclo_dose_qual <- factor(meta_lon$cyclo_dose_qual,
    levels = c(present_pref, setdiff(unique(meta_lon$cyclo_dose_qual), present_pref))
)

desired_levels <- c("none", "low", "high", "dorsal")
final_levels <- c(
    intersect(desired_levels, unique(meta_lon$cyclo_dose_qual)),
    setdiff(unique(meta_lon$cyclo_dose_qual), desired_levels)
)
meta_lon$cyclo_dose_qual <- factor(meta_lon$cyclo_dose_qual, levels = final_levels)
ordered_samples <- with(meta_lon, sample[order(match(cyclo_dose_qual, final_levels), sample)])

# create a simple column annotation (cyclo_dose_qual) and colors consistent with previous plot
palette_base <- c("#80AD3C", "#eaac30", "#cf6509", "#A1A1DE")
levels_cdq <- levels(meta_lon$cyclo_dose_qual)
ann_colors <- setNames(palette_base[seq_len(length(levels_cdq))], levels_cdq)

# ensure annotation rows follow the same ordering as the heatmap columns
ordered_cdq <- meta_lon$cyclo_dose_qual[match(ordered_samples, meta_lon$sample)]

col_ha <- HeatmapAnnotation(
    cyclo_dose_qual = ordered_cdq,
    col = list(cyclo_dose_qual = ann_colors),
    annotation_legend_param = list(title = "cyclo_dose_qual")
)


# scale rows (z-score); handle constant rows by replacing NA with 0
scaled_mat <- t(scale(t(lon_mat)))
scaled_mat[is.na(scaled_mat)] <- 0

# build heatmap object for later drawing (use same parameters as the heatmap call below)
DEG_ht_plot <- Heatmap(
    scaled_mat[rownames(scaled_mat), ordered_samples],
    name = "Normalized expression",
    row_title = NULL,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = 6),
    show_column_names = TRUE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    width = ncol(scaled_mat) * unit(4, "mm"),
    top_annotation = col_ha,
    col = colorRampPalette(c("black", "purple", "orange", "yellow"))(1000)
)
DEG_ht_plot
# draw and save
png("article/transit_data/3p_genes_heatmap.png", width = 10, height = 6, units = "in", res = 300)
draw(DEG_ht_plot, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

consensus_positions_with_stats$Gene %>%
    unique() %>%
    length()
