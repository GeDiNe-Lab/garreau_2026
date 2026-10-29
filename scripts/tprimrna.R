# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(tibble)

setwd(paste(here::here(), "projects", sep = "/"))

# Loading custom functions
source("article/scripts/custom_fct.R")

# Loading panasenkava data
rawcounts <- readcounts("panasenkava_2024/data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("panasenkava_2024/data/meta.csv", sep = ",", header = TRUE)

# Filtering metadata and counts for the analysis
meta <- filter(rawmeta, CRISPR %in% c("no", "control"), type %in% c("ventral", "cyclo"), sample != "L9C1_2", diff %in% c("diff9", "diff12"))
meta$cyclo_dose_quant_factor <- factor(meta$cyclo_dose_quant, levels = c(0, 0.125, 0.25, 0.5, 1))
rawcounts <- rawcounts[, meta$sample]
counts <- rawcounts[which(rowSums(rawcounts) >= 25), ] # removing lowly expressed genes

# Loading Zeng et al. single cell expressed genes
load("article/transit_data/sc_genes.RData")

# Normalizing counts
norm <- DESeqDataSetFromMatrix(
    countData = counts[which(gene_converter(rownames(counts), "ENSEMBL", "SYMBOL") %in% sc_genes == TRUE), ],
    colData = meta,
    design = ~ line * cyclo_dose_qual
) %>%
    vst(blind = FALSE)

# Removing batch effect
assay(norm) <- limma::removeBatchEffect(assay(norm), batch = meta$line)

# PCA analysis
pca.data <- plotPCA.DESeqTransform(norm, intgroup = c("type", "cyclo_dose_quant_factor", "line"), returnData = TRUE, ntop = 3000)
percentVar <- round(100 * attr(pca.data, "percentVar"))
pca_plot <- ggplot(pca.data, aes(PC1, PC2, color = cyclo_dose_quant_factor, shape = line)) +
    geom_point(size = 3, stroke = 2) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#591b00", "#a63e01", "#b76718", "#e0b45c", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme()
ggsave("article/figures/Figure_1/Figure1supp_PCA.png", plot = pca_plot, width = 10, height = 5)

# Get PCA/covariates ANOVA results
PC_covariate_ANOVA <- pca_anova(
    pca_data = pca.data,
    metadata = meta,
    covariates = c("line", "cyclo_dose_quant_factor")
)
write.csv(PC_covariate_ANOVA, "article/tables/PCA_ANOVA.csv")

# PCs variation percentages :
percentVar_plot <- ggplot(data.frame(perc = percentVar, PC = factor(colnames(pca.data[1:15]), levels = colnames(pca.data[1:15]))), aes(x = PC, y = perc)) +
    geom_bar(stat = "identity") +
    custom_theme(diag_text = TRUE) +
    ylim(0, 100) +
    ylab("% of variance explained") +
    xlab("Principal Components")
ggsave("results/images/Figure_1/Figure1supp_percentVar.png", plot = percentVar_plot, width = 10, height = 5)


# Network construction
# sft <- pickSoftThreshold(t(norm),
#     powerVector = c(1:20),
#     verbose = 5
# )

# filter(sft$fitIndices, SFT.R.sq > 0.8, mean.k. < 100)
coex <- cor(t(assay(norm)))**6 # Soft thresholding with power 6
raw_net <- abs(coex) %>% as.matrix() # Switching to absolute values
raw_net[abs(raw_net) < 0.1] <- 0 # Filtering out weak correlations
keep <- which(rowSums(raw_net != 0) > 1) # Filtering isolated genes
net <- raw_net[keep, keep]

# Get context genes for ORA
rownames(norm) %>%
    gene_converter("ENSEMBL", "SYMBOL") %>%
    write.csv("article/tables/context_genes.txt", row.names = FALSE, quote = FALSE)

# Saving data
data_3p <- list(
    context_genes = norm %>%
        rownames(),
    norm = norm,
    meta = meta,
    coex = net,
    DEGs = DEGs
)
save(data_3p, file = "article/transit_data/data_3p.RData")
