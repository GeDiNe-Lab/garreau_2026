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
source("article/scripts/RECENT_patient_gene_pvalue.R")
source("article/scripts/RECENT_webgestalt_ORA.R")
source("article/scripts/RECENT_rwr_normalization.R")

load("article/transit_data/RWR_results_case.RData")
HPE_results <- as.matrix(RWR_results_z)
rm(RWR_results_z)
load("article/transit_data/RWR_iron_results_case.RData")
FER_results <- as.matrix(RWR_results_z)
rm(RWR_results_z)
load("article/transit_data/iop_RWR_results_case.RData")
IOP_results <- as.matrix(RWR_results_z)
rm(RWR_results_z)
load("article/transit_data/data_3p.RData")
load("article/transit_data/iron_patient_variants.RData")
FER_patient_variants <- patient_variants
rm(patient_variants)
load("article/transit_data/patient_variants.RData")
HPE_patient_variants <- patient_variants
rm(patient_variants)
load("article/transit_data/iop_patient_variants.RData")
IOP_patient_variants <- patient_variants
rm(patient_variants)
load("article/transit_data/permutation_scores.RData")


coex_genes <- rownames(data_3p$coex_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")
HPE_seeds <- lapply(HPE_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(HPE_seeds) <- colnames(HPE_results)

IOP_seeds <- lapply(IOP_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(IOP_seeds) <- colnames(IOP_results)

FER_seeds <- lapply(FER_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(FER_seeds) <- colnames(FER_results)

RWR_net <- data_3p$coex_case
colnames(RWR_net) <- rownames(RWR_net) <- gene_converter(rownames(RWR_net), from = "ENSEMBL", to = "SYMBOL")


context <- read.csv("article/tables/context_genes.txt") %>% pull(x)


# ============================================================
# Exemple d'utilisation
# ============================================================
# 1. Construire la banque UNE FOIS (ex: 2000 permutations, tailles de
#    seeds tirees de la distribution reelle observee chez tes patients)
# observed_seed_counts <- sapply(HPE_seeds, length)
# # remove zeros
# observed_seed_counts <- observed_seed_counts[observed_seed_counts > 0]

# bank <- build_permutation_bank(
#     RWR_net,
#     seed_sizes = observed_seed_counts,
#     n_permutations = 2000, mc.cores = 4
# )
# bank
# save(bank, file = "article/transit_data/permutation_bank.RData")
load("article/transit_data/permutation_bank.RData")

# ============================================================
# Random sampling stability analysis: HPO ORA on signal genes
# from 100 random patient groups of size 20
# ============================================================
#
# Assumes the following objects/functions already exist in the
# session (loaded from your scripts):
#   - RWR_results   : gene x patient score matrix
#   - seeds         : named list of seed genes per patient
#   - process_group : your function (defined above)
#   - run_ora_hpo   : your WebGestaltR wrapper (defined above)
#
# Adjust variable names below if your objects are named differently.

# HPE results without first patient
RWR_results <- HPE_results[, -1]
seeds <- HPE_seeds

length(seeds)

set.seed(123)
n_iter <- 100
sample_size <- 15
all_patients <- colnames(RWR_results)

if (length(all_patients) < sample_size) {
    stop("Moins de patients disponibles que la taille d'échantillon demandée (", sample_size, ")")
}

hpo_enrichment_list <- vector("list", n_iter)

for (i in seq_len(n_iter)) {
    cat("=== Iteration", i, "/", n_iter, "===\n")

    sampled_patients <- sample(all_patients, sample_size)

    # --- Step 1: get signal genes for this random group ---
    group_result <- tryCatch(
        process_group(
            patient_ids = sampled_patients,
            RWR_results = RWR_results,
            seeds = seeds,
            permutation_bank = bank
        ),
        error = function(e) {
            message("Iteration ", i, " : erreur process_group - ", conditionMessage(e))
            return(NULL)
        }
    )

    if (is.null(group_result)) {
        message("Iteration ", i, " : process_group a retourne NULL (ou erreur), ignoree")
        next
    }

    signal_genes <- group_result$gene[group_result$signal_final]

    if (length(signal_genes) < 5) {
        message(
            "Iteration ", i, " : moins de 5 genes signal (",
            length(signal_genes), "), ORA ignoree"
        )
        next
    }

    # --- Step 2: run HPO ORA on the signal genes ---
    hpo_res <- tryCatch(
        run_ora_hpo(
            genelist = signal_genes,
            projectName = paste0("ORA_HPO_iter", i),
            is.output = FALSE,
            referenceGene = context
        ),
        error = function(e) {
            message("Iteration ", i, " : erreur ORA - ", conditionMessage(e))
            return(NULL)
        }
    )

    if (!is.null(hpo_res) && nrow(hpo_res) > 0) {
        hpo_res$iteration <- i
        hpo_res$n_genes <- length(signal_genes)
        hpo_enrichment_list[[i]] <- hpo_res
    }
}

# ============================================================
# Combine all iterations and summarise
# ============================================================
# NOTE: WebGestaltR ORA output typically includes columns such as
# "description" (term name), "geneSet" (term ID), and
# "enrichmentRatio". Adjust the column names below (description /
# enrichmentRatio / geneSet) if your WebGestaltR version differs
# (check colnames(hpo_enrichment_all)).

hpo_enrichment_all <- dplyr::bind_rows(hpo_enrichment_list)
colnames(hpo_enrichment_all)
hpo_enrichment_all
if (nrow(hpo_enrichment_all) == 0) {
    stop("Aucun terme HPO enrichi trouve sur les ", n_iter, " iterations.")
}

hpo_summary <- hpo_enrichment_all %>%
    dplyr::group_by(geneSet, description) %>%
    dplyr::summarise(
        n_iterations = dplyr::n_distinct(iteration),
        pct_iterations = 100 * n_iterations / n_iter,
        mean_enrichmentRatio = mean(enrichmentRatio, na.rm = TRUE),
        sd_enrichmentRatio = sd(enrichmentRatio, na.rm = TRUE),
        mean_FDR = mean(FDR, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_iterations), dplyr::desc(mean_enrichmentRatio))

print(hpo_summary, n = 30)

# IOP_hpo_sum <- hpo_summary
View(IOP_hpo_sum)

# HPE_hpo_sum <- hpo_summary
View(HPE_hpo_sum)

# FER_hpo_sum <- hpo_summary
View(FER_hpo_sum)

# Optional: save results
write.csv(hpo_enrichment_all, "article/enrichment/hpo_random_sampling_all_iterations.csv", row.names = FALSE)
write.csv(hpo_summary, "article/enrichment/hpo_random_sampling_summary.csv", row.names = FALSE)
