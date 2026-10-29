# Loading packages and functions
library(Matrix)
library(tidyverse)
library(parallel)
library(tibble)
library(WGCNA)
library(RANKS)
library(igraph)

# Setting working directory
setwd(paste(here::here(), "projects", sep = "/"))

# Loading custom functions
source("article/scripts/custom_fct.R")
get_patient_seeds <- function(patient_variants, RWR_net,
                              gene_col = "Gene.refGeneWithVer") {
    seeds_list <- lapply(patient_variants, function(patient_df) {
        patient_genes <- intersect(rownames(RWR_net), patient_df[[gene_col]])
        if (length(patient_genes) == 0) {
            return(character(0))
        }
        return(patient_genes)
    })

    if (is.null(names(seeds_list))) {
        names(seeds_list) <- paste0("Patient_", seq_along(seeds_list))
    }

    return(seeds_list)
}

compute_RWR_matrix <- function(RWR_net, seeds_list,
                               gamma = 0.7, tmax = 1000, eps = 1e-10,
                               norm = TRUE, mc.cores = 1, verbose = TRUE) {
    if (is.null(names(seeds_list))) {
        names(seeds_list) <- paste0("Patient_", seq_along(seeds_list))
    }
    patient_ids <- names(seeds_list)

    run_one <- function(pid) {
        seeds <- seeds_list[[pid]]

        if (verbose) {
            message(
                "Processing patient: ", pid
            )
        }

        if (length(seeds) == 0) {
            if (verbose) message("  -> pas de seeds pour ", pid, ", NA")
            return(rep(NA_real_, nrow(RWR_net)))
        }

        score <- RWR(RWR_net, seeds,
            gamma = gamma, tmax = tmax, eps = eps, norm = norm
        )$p %>%
            unname()

        return(score)
    }

    if (mc.cores > 1) {
        score_list <- mclapply(patient_ids, run_one, mc.cores = mc.cores)
    } else {
        score_list <- lapply(patient_ids, run_one)
    }

    RWR_results <- Reduce(cbind, score_list)
    colnames(RWR_results) <- patient_ids
    rownames(RWR_results) <- rownames(RWR_net)

    return(RWR_results)
}

load("article/transit_data/data_3p.RData")

RWR_net <- data_3p$coex_case
colnames(RWR_net) <- rownames(RWR_net) %>% gene_converter("ENSEMBL", "SYMBOL")
rownames(RWR_net) <- rownames(RWR_net) %>% gene_converter("ENSEMBL", "SYMBOL")

load("article/transit_data/hpe_patient_variants.RData")
HPE_seeds <- get_patient_seeds(HPE_patient_variants, RWR_net, gene_col = "Gene.refGeneWithVer")

load("article/transit_data/poi_patient_variants.RData")
IOP_seeds <- get_patient_seeds(POI_patient_variants, RWR_net, gene_col = "Gene.refGeneWithVer")

# Performing RWR for HPE and IOP patient seeds
HPE_results <- compute_RWR_matrix(RWR_net, HPE_seeds, mc.cores = 8)
save(HPE_results, file = "article/transit_data/HPE_RWR_results_cds.RData")
IOP_results <- compute_RWR_matrix(RWR_net, IOP_seeds, mc.cores = 8)
save(IOP_results, file = "article/transit_data/IOP_RWR_results_cds.RData")