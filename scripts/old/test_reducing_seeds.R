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
library(plotly)
library(WebGestaltR)
library(scales)

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
source("article/scripts/RECENT_patient_gene_pvalue.R") # process_group() / process_patient()
source("article/scripts/RECENT_rwr_normalization.R") # build_permutation_bank() / normalize_rwr_score()
source("article/scripts/RECENT_webgestalt_ORA.R")

load("article/transit_data/data_3p.RData")

load("article/transit_data/RWR_results_case.RData")
HPE_results <- as.matrix(RWR_results_z)
rm(RWR_results_z)
load("article/transit_data/iop_RWR_results_case.RData")
IOP_results <- as.matrix(RWR_results_z)
rm(RWR_results_z)

load("article/transit_data/iop_patient_variants.RData")
IOP_patient_variants <- patient_variants
rm(patient_variants)
load("article/transit_data/patient_variants.RData")
HPE_patient_variants <- patient_variants
rm(patient_variants)

HPO <- c(
    "B3GALNT2",
    "B9D1",
    "B9D2",
    "BUB1",
    "BUB1B",
    "BUB3",
    "C2CD3",
    "CC2D2A",
    "CDON",
    "CEP290",
    "CEP57",
    "CHD7",
    "CILK1",
    "CNOT1",
    "CRIPTO",
    "CRPPA",
    "CSPP1",
    "DHCR7",
    "DISP1",
    "DLL1",
    "FGF8",
    "FGFR1",
    "FGFR3",
    "FILIP1",
    "FKRP",
    "FKTN",
    "FOXA2",
    "FOXH1",
    "GAS1",
    "GLI2",
    "GLI3",
    "GMPPB",
    "GPKOW",
    "HESX1",
    "HOXA2",
    "KDM6A",
    "KMT2D",
    "L1CAM",
    "LARGE1",
    "LHX4",
    "MKS1",
    "MTHFR",
    "NODAL",
    "OTX2",
    "PLCH1",
    "POMGNT1",
    "POMT1",
    "POMT2",
    "POU1F1",
    "PPP1R12A",
    "PROP1",
    "PRRX1",
    "PTCH1",
    "RAD21",
    "RB1",
    "RPGRIP1",
    "RPGRIP1L",
    "RTTN",
    "SALL1",
    "SEC31A",
    "SEMA3E",
    "SHH",
    "SIX3",
    "SMC1A",
    "SOX2",
    "STAG2",
    "SUFU",
    "TBC1D32",
    "TCTN1",
    "TCTN2",
    "TCTN3",
    "TGIF1",
    "TMEM107",
    "TMEM216",
    "TMEM231",
    "TMEM237",
    "TMEM67",
    "TRIP13",
    "TXNDC15",
    "VANGL2",
    "WDR81",
    "ZIC2",
    "ZRSR2"
)
# Precomputed permutation bank (build_permutation_bank(), cf. RECENT_results_analyses.R)
# -> reused as-is for both cohorts, no RWR recomputation
load("article/transit_data/permutation_bank.RData") # -> `bank`

coex_genes <- rownames(data_3p$coex_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

context_genes <- rownames(data_3p$norm_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

HPE_seeds <- lapply(HPE_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(HPE_seeds) <- colnames(HPE_results)

IOP_seeds <- lapply(IOP_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(IOP_seeds) <- colnames(IOP_results)

data_3p$DEGs_case$gene <- gene_converter(rownames(data_3p$DEGs_case), "ENSEMBL", "SYMBOL")

DEGs <- data_3p$DEGs_case %>% filter(abs(log2FoldChange) >= 1)

HPE_seeds_f <- lapply(HPE_seeds, function(x) {
    return(x[which(x %in% HPO)])
})
HPE_seeds_f


IOP_seeds_f <- lapply(IOP_seeds, function(x) {
    return(x[which(x %in% HPO)] %>% length())
})
IOP_seeds_f %>%
    unlist() %>%
    hist()
