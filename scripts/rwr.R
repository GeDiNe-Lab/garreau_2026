# Loading packages and functions
library(Matrix)
library(tidyverse)
library(parallel)
library(tibble)
library(WGCNA)
library(RANKS)
library(igraph)

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

# ============================================================
# get_patient_seeds()
# ------------------------------------------------------------
# Extrait, pour chaque patient, les genes de variants qui
# intersectent avec les genes du reseau (= seeds RWR).
#
# @param patient_variants  liste de data.frames (un par patient),
#                          chacun contenant une colonne de genes
# @param RWR_net           matrice d'adjacence du reseau (genes en
#                          rownames/colnames, deja convertis en SYMBOL)
# @param gene_col           nom de la colonne contenant les genes
#                          dans chaque data.frame patient
#
# @return liste nommee de vecteurs de genes (seeds), un par patient
# ============================================================
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

# ============================================================
# build_seed_weights()
# ------------------------------------------------------------
# Construit un vecteur de poids (NON normalise -- RWR() s'en charge
# en interne via seed.weights/sum(seed.weights)) pour un jeu de
# seeds, a partir d'une table DEG.
#
# Poids d'une seed DE       = 2^|log2FoldChange|  (fold-change reel :
#                              le "n fois plus" demande)
# Poids d'une seed non DE    = 1 (base neutre, 2^0), i.e. seed absente
# (absente de deg_table)       de la table DEG (variant sans signal
#                              d'expression differentielle)
#
# @param seeds       vecteur de genes seeds pour ce patient
# @param deg_table   data.frame DEG deja filtre (ex: FDR <= 0.05), avec
#                    au moins les colonnes gene_col et lfc_col
# @param gene_col    nom de la colonne genes dans deg_table
# @param lfc_col     nom de la colonne log2FoldChange dans deg_table
# @param use_abs     utiliser |log2FC| (TRUE, defaut) plutot que le
#                    log2FC signe -- on pondere ici l'AMPLEUR du
#                    changement d'expression, pas son sens
#
# @return vecteur numerique, DANS LE MEME ORDRE que `seeds`, pret a
#         etre passe a RWR(..., seed.weights = ...)
# ============================================================
build_seed_weights <- function(seeds, deg_table,
                               gene_col = "gene", lfc_col = "log2FoldChange",
                               use_abs = TRUE) {
    lfc <- deg_table[[lfc_col]][match(seeds, deg_table[[gene_col]])]
    lfc[is.na(lfc)] <- 0 # seed non DE -> log2FC = 0 -> poids de base 2^0 = 1

    if (use_abs) lfc <- abs(lfc)

    w <- 2^lfc # pas de normalisation ici : RWR() le fait (seed.weights/sum(seed.weights))

    return(w)
}


# ============================================================
# compute_RWR_matrix()
# ------------------------------------------------------------
# Calcule le RWR pour une liste de patients (chacun defini par son
# propre vecteur de seeds), avec option de ponderation des seeds par
# log2FoldChange (table DEG), via le parametre natif seed.weights de
# RANKS::RWR(). Retourne une matrice genes x patients.
#
# @param RWR_net      matrice d'adjacence du reseau
# @param seeds_list   liste nommee de vecteurs de genes seeds, un par patient
# @param deg_list     optionnel. Soit :
#                       - NULL (defaut) : pas de ponderation, comportement
#                         original (seeds a poids uniforme)
#                       - un data.frame DEG UNIQUE, applique a tous les
#                         patients (meme table DEG globale)
#                       - une liste nommee de data.frames DEG, un par
#                         patient (meme structure que seeds_list)
# @param gene_col, lfc_col, use_abs  transmis a build_seed_weights()
# @param gamma, tmax, eps, norm  parametres transmis a RANKS::RWR()
# @param mc.cores     nombre de coeurs (mclapply si > 1)
# @param verbose      afficher la progression
#
# @return matrice numerique genes x patients
# ============================================================
compute_RWR_matrix <- function(RWR_net, seeds_list,
                               deg_list = NULL,
                               gene_col = "gene", lfc_col = "log2FoldChange", use_abs = TRUE,
                               gamma = 0.9, tmax = 1000, eps = 1e-10,
                               norm = TRUE, mc.cores = 1, verbose = TRUE) {
    if (is.null(names(seeds_list))) {
        names(seeds_list) <- paste0("Patient_", seq_along(seeds_list))
    }
    patient_ids <- names(seeds_list)

    weighted <- !is.null(deg_list)
    global_deg <- if (weighted && is.data.frame(deg_list)) deg_list else NULL

    run_one <- function(pid) {
        seeds <- seeds_list[[pid]]

        if (verbose) {
            message(
                "Processing patient: ", pid,
                if (weighted) " (seeds ponderees log2FC)" else ""
            )
        }

        if (length(seeds) == 0) {
            if (verbose) message("  -> pas de seeds pour ", pid, ", NA")
            return(rep(NA_real_, nrow(RWR_net)))
        }

        seed_weights <- NULL

        if (weighted) {
            deg_table <- if (!is.null(global_deg)) global_deg else deg_list[[pid]]

            if (is.null(deg_table)) {
                if (verbose) {
                    message(
                        "  -> pas de table DEG pour ", pid,
                        ", poids uniformes (equivalent non pondere)"
                    )
                }
            } else {
                seed_weights <- build_seed_weights(seeds, deg_table,
                    gene_col = gene_col, lfc_col = lfc_col,
                    use_abs = use_abs
                )
            }
        }

        score <- RWR(RWR_net, seeds,
            gamma = gamma, tmax = tmax, eps = eps, norm = norm,
            seed.weights = seed_weights
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

load("article/transit_data/patient_variants.RData")
HPE_seeds <- get_patient_seeds(patient_variants, RWR_net, gene_col = "Gene.refGeneWithVer")
rm(patient_variants)
load("article/transit_data/iop_patient_variants.RData")
IOP_seeds <- get_patient_seeds(patient_variants, RWR_net, gene_col = "Gene.refGeneWithVer")
rm(patient_variants)

# load("article/transit_data/iron_patient_variants.RData")
# FER_seeds <- get_patient_seeds(patient_variants, RWR_net, gene_col = "Gene.refGeneWithVer")
# rm(patient_variants)

# Avec ponderation, table DEG globale identique pour tous les patients
# deg_table : data.frame avec colonnes "gene" et "log2FoldChange"
HPE_results <- compute_RWR_matrix(RWR_net, HPE_seeds, mc.cores = 4)
save(HPE_results, file = "article/transit_data/HPE_RWR_results.RData")

IOP_results <- compute_RWR_matrix(RWR_net, IOP_seeds, mc.cores = 4)
save(IOP_results, file = "article/transit_data/IOP_RWR_results.RData")

# FER_results <- compute_RWR_matrix(RWR_net, FER_seeds, mc.cores = 4)
# save(FER_results, file = "article/transit_data/FER_RWR_results.RData")
