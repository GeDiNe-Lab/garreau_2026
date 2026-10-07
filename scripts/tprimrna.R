# Loading packages and functions
library(Matrix)
library(tidyverse)
library(DESeq2)
library(tibble)

setwd(paste(here::here(), "projects", sep = "/"))

# Loading custom functions
source("article/scripts/custom_fct.R")

rawcounts <- readcounts("panasenkava_2024/data/rawcounts.csv", sep = ",", header = TRUE)
rawmeta <- read.table("panasenkava_2024/data/meta.csv", sep = ",", header = TRUE)

meta_case <- filter(rawmeta, CRISPR %in% c("no", "control"), type %in% c("ventral", "cyclo"), sample != "L9C1_2", diff %in% c("diff9", "diff12"))
meta_case$cyclo_dose_quant_factor <- factor(meta_case$cyclo_dose_quant, levels = c(0, 0.125, 0.25, 0.5, 1))
meta_case$cyclo_dose_quant_factor
meta_case$cyclo_dose_qual
rawcounts_case <- rawcounts[, meta_case$sample]
counts_case <- rawcounts_case[which(rowSums(rawcounts_case) >= 25), ]

load("article/transit_data/sc_genes.RData")
nrow(counts_case)
nrow(counts_case[which(gene_converter(rownames(counts_case), "ENSEMBL", "SYMBOL") %in% sc_genes == TRUE), ])

norm_case <- DESeqDataSetFromMatrix(
    countData = counts_case[which(gene_converter(rownames(counts_case), "ENSEMBL", "SYMBOL") %in% sc_genes == TRUE), ],
    colData = meta_case,
    design = ~ line * cyclo_dose_qual
) %>%
    vst(blind = FALSE) %>%
    assay() %>%
    limma::removeBatchEffect(batch = meta_case$line)

norm_case <- DESeqDataSetFromMatrix(
    countData = counts_case[which(gene_converter(rownames(counts_case), "ENSEMBL", "SYMBOL") %in% sc_genes == TRUE), ],
    colData = meta_case,
    design = ~ line * cyclo_dose_qual
) %>%
    vst(blind = FALSE)
assay(norm_case) <- limma::removeBatchEffect(assay(norm_case), batch = meta_case$line)

pca.data <- plotPCA.DESeqTransform(norm_case, intgroup = c("type", "cyclo_dose_quant_factor", "line"), returnData = TRUE, ntop = 3000)
percentVar <- round(100 * attr(pca.data, "percentVar"))

pca_plot <- ggplot(pca.data, aes(PC1, PC2, color = cyclo_dose_quant_factor, shape = line)) +
    geom_point(size = 3, stroke = 2) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c("#591b00", "#a63e01", "#b76718", "#e0b45c", "#80AD3C")) +
    scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6)) +
    custom_theme()
pca_plot

# Get PCA/covariates ANOVA results
PC_covariate_ANOVA <- pca_anova(
    pca_data = pca.data,
    metadata = meta_case,
    covariates = c("line", "cyclo_dose_quant_factor")
)
#  Saving ANOVA results
write.csv(PC_covariate_ANOVA, "article/tables/PCA_ANOVA.csv")

# PCs variation percentages :
percentVar_plot <- ggplot(data.frame(perc = percentVar, PC = factor(colnames(pca.data[1:15]), levels = colnames(pca.data[1:15]))), aes(x = PC, y = perc)) +
    geom_bar(stat = "identity") +
    custom_theme(diag_text = TRUE) +
    ylim(0, 100) +
    ylab("% of variance explained") +
    xlab("Principal Components")
ggsave("results/images/Figure_1/Figure1supp_percentVar.png", plot = percentVar_plot, width = 10, height = 5)

# sft_case <- pickSoftThreshold(t(norm_case),
#     powerVector = c(1:20),
#     verbose = 5
# )

# filter(sft_case$fitIndices, SFT.R.sq > 0.8, mean.k. < 100)

coex_case <- cor(t(norm_case))**6

raw_case_net <- abs(coex_case) %>% as.matrix()
raw_case_net[abs(raw_case_net) < 0.1] <- 0
case_keep <- which(rowSums(raw_case_net != 0) > 1)
case_net <- raw_case_net[case_keep, case_keep]

nrow(case_net)

rownames(norm_case) %>%
    gene_converter("ENSEMBL", "SYMBOL") %>%
    write.csv("article/tables/context_genes.txt", row.names = FALSE, quote = FALSE)

data_3p <- list(
    context_genes = rawcounts_case[which(rowSums(rawcounts_case) > 0), ] %>%
        rownames(),
    norm_case = norm_case,
    meta_case = meta_case,
    coex_case = case_net,
    DEGs_case = DEGs_case
)
data_3p$context_genes
save(data_3p, file = "article/transit_data/data_3p.RData")
