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

# Loading custom functions
source("article/scripts/custom_fct.R")
source("article/scripts/RECENT_rwr_normalization.R") # build_permutation_bank() / normalize_rwr_score()

my_locfdr_robust <- function(z) {
    z <- z[is.finite(z) & !is.na(z)]

    # --- Garde-fous : pas assez de points, ou distribution degeneree ---
    if (length(z) < 30) {
        message(
            "my_locfdr_robust : moins de 30 valeurs finies (", length(z),
            "), N(0,1) utilisee sans scaling robuste"
        )
        log_pval <- pnorm(z, mean = 0, sd = 1, lower.tail = FALSE, log.p = TRUE)
        names(log_pval) <- names(z)
        return(log_pval)
    }

    mad_z <- mad(z)
    if (!is.finite(mad_z) || mad_z == 0) {
        message(
            "my_locfdr_robust : MAD nul ou non fini (distribution degeneree, ",
            "trop de valeurs identiques) — fallback N(0,1) sur z brut"
        )
        log_pval <- pnorm(z, mean = 0, sd = 1, lower.tail = FALSE, log.p = TRUE)
        names(log_pval) <- names(z)
        return(log_pval)
    }

    z_rob <- (z - median(z)) / mad_z

    N <- length(z_rob)
    breaks <- seq(min(z_rob), max(z_rob), length.out = max(30, round(sqrt(N) * 10 / 2)) + 1)
    mids <- (breaks[-1] + breaks[-length(breaks)]) / 2
    counts <- as.numeric(table(cut(z_rob, breaks = breaks, include.lowest = TRUE)))
    bin_width <- diff(breaks)[1]

    # --- Garde-fou supplementaire : bin_width nul (min == max apres arrondi) ---
    if (!is.finite(bin_width) || bin_width == 0) {
        message("my_locfdr_robust : bin_width nul — fallback N(0,1) sur z brut")
        log_pval <- pnorm(z, mean = 0, sd = 1, lower.tail = FALSE, log.p = TRUE)
        names(log_pval) <- names(z)
        return(log_pval)
    }

    fit_f <- glm(counts ~ splines::ns(mids, df = 7),
        family  = poisson(),
        offset  = rep(log(N * bin_width), length(mids))
    )
    f_hat <- fitted(fit_f) / (N * bin_width)

    central <- mids >= quantile(z_rob, 0.25) & mids <= quantile(z_rob, 0.75)
    valid <- central & is.finite(f_hat) & f_hat > 0

    if (sum(valid) < 3) {
        message("Pas assez de bins valides pour estimer la null — N(0,1) utilisée")
        mu0 <- 0
        sigma0 <- 1
    } else {
        fit_null <- lm(log(f_hat[valid] + 1e-10) ~ mids[valid] + I(mids[valid]^2))
        cc <- coef(fit_null)[3]
        b <- coef(fit_null)[2]

        if (is.na(cc) || cc >= 0) {
            mu0 <- 0
            sigma0 <- 1
        } else {
            sigma0 <- sqrt(-1 / (2 * cc))
            mu0 <- -b / (2 * cc)
        }
    }

    log_pval <- pnorm(z_rob, mean = mu0, sd = sigma0, lower.tail = FALSE, log.p = TRUE)
    return(log_pval)
}


# ============================================================
# process_patient()
# ------------------------------------------------------------
# CHANGEMENT : HPE_results est desormais attendu en scores RWR
# BRUTS (plus de Z-score en amont). La normalisation se fait ici,
# au niveau du patient (N = 1), via la banque de permutations
# precalculee (build_permutation_bank), reutilisee sans recalcul RWR.
# ============================================================
process_patient <- function(patient_id, HPE_results, HPE_seeds, permutation_bank) {
    print(paste("Processing patient:", patient_id))

    # Extraction des scores bruts sans seeds
    score_wo_seeds <- HPE_results[
        !rownames(HPE_results) %in% HPE_seeds[[patient_id]],
        patient_id
    ]

    # Filtre NA
    if (all(is.na(score_wo_seeds))) {
        message("Patient ", patient_id, " : que des NA, ignoré")
        return(NULL)
    }
    score_wo_seeds <- score_wo_seeds[!is.na(score_wo_seeds)]

    # --- Normalisation Z-score via banque de permutations (niveau patient unique, N = 1) ---
    z_scores <- normalize_rwr_score(
        score_wo_seeds,
        permutation_bank = permutation_bank,
        verbose = FALSE
    )

    # Log p-values + correction, calcules sur les scores normalises
    log_pvalues <- my_locfdr_robust(z_scores)
    pvalues_safe <- exp(pmax(log_pvalues, -745))
    pvalues_corrected <- p.adjust(pvalues_safe, method = "BH")
    threshold_95 <- quantile(z_scores, 0.90, na.rm = TRUE)

    data.frame(
        patient = patient_id,
        gene = names(z_scores),
        score_raw = score_wo_seeds[names(z_scores)],
        z_score = z_scores,
        log_pval = log_pvalues,
        pval_bh = pvalues_corrected,
        signal_locfdr = pvalues_corrected < 0.05,
        signal_95 = z_scores >= threshold_95,
        signal_final = pvalues_corrected < 0.05 & z_scores >= threshold_95
    )
}


# ============================================================
# process_group()
# ------------------------------------------------------------
# CHANGEMENT : RWR_results est desormais attendu en scores RWR
# BRUTS. L'imputation des seeds (median des non-seeds) se fait
# toujours en amont sur les scores bruts. La normalisation Z-score
# est ensuite appliquee APRES le regroupement, directement sur la
# matrice scores_imputed (genes x patients valides) : normalize_rwr_score()
# fait la moyenne ET construit un null matche sur N = nb de patients
# valides (ncol(scores_imputed)), en reechantillonnant la banque de
# permutations precalculee (aucun recalcul RWR).
# ============================================================
process_group <- function(patient_ids = NULL, RWR_results, seeds, permutation_bank) {
    # Par défaut : tous les patients
    if (is.null(patient_ids)) patient_ids <- colnames(RWR_results)

    # Vérification
    invalid <- patient_ids[!patient_ids %in% colnames(RWR_results)]
    if (length(invalid) > 0) warning("Patients non trouvés : ", paste(invalid, collapse = ", "))
    patient_ids <- patient_ids[patient_ids %in% colnames(RWR_results)]

    # Pour chaque patient : remplacer les scores des seeds par la médiane des non-seeds (scores bruts)
    scores_imputed <- sapply(patient_ids, function(pid) {
        scores <- RWR_results[, pid]
        seed_genes <- seeds[[pid]]
        is_seed <- rownames(RWR_results) %in% seed_genes

        # Médiane calculée sans les seeds
        median_nonseed <- median(scores[!is_seed], na.rm = TRUE)

        # Remplacement des scores seeds par la médiane
        scores[is_seed] <- median_nonseed

        return(scores)
    })

    # scores_imputed : matrice gènes x patients
    rownames(scores_imputed) <- rownames(RWR_results)

    # Filtre patients avec trop de NA
    na_prop <- colMeans(is.na(scores_imputed))
    excluded <- patient_ids[na_prop > 0.5]
    if (length(excluded) > 0) message("Patients exclus (>50% NA) : ", paste(excluded, collapse = ", "))
    scores_imputed <- scores_imputed[, na_prop <= 0.5, drop = FALSE]

    if (ncol(scores_imputed) == 0) {
        message("Aucun patient valide dans ce groupe")
        return(NULL)
    }

    # Score moyen brut par gène, conserve pour reporting/comparaison
    mean_scores_raw <- rowMeans(scores_imputed, na.rm = TRUE)

    # --- Normalisation Z-score via banque de permutations ---
    # scores_imputed (genes x N patients valides) est passe directement :
    # normalize_rwr_score fait la moyenne des patients ET construit le null
    # en agregeant N profils tires de la banque (N = ncol(scores_imputed)),
    # donc le null correspond exactement a la statistique testee (moyenne de groupe).
    z_scores <- normalize_rwr_score(
        scores_imputed,
        permutation_bank = permutation_bank,
        verbose = FALSE
    )
    z_scores <- z_scores[is.finite(z_scores)]

    # Pipeline locfdr + correction, sur les z-scores de groupe
    log_pvalues <- my_locfdr_robust(z_scores)
    pvalues_safe <- exp(pmax(log_pvalues, -745))
    pvalues_corrected <- p.adjust(pvalues_safe, method = "BH")
    threshold_90 <- quantile(z_scores, 0.75, na.rm = TRUE)

    df <- data.frame(
        gene           = names(z_scores),
        mean_score_raw = mean_scores_raw[names(z_scores)],
        z_score        = z_scores,
        n_patients     = ncol(scores_imputed),
        log_pval       = log_pvalues,
        pval_bh        = pvalues_corrected,
        signal_locfdr  = pvalues_corrected < 0.05,
        signal_90      = z_scores >= threshold_90,
        signal_final   = pvalues_corrected < 0.05 & z_scores >= threshold_90
    )

    cat(
        "Groupe :", length(patient_ids), "patients demandés |",
        ncol(scores_imputed), "patients valides |",
        sum(df$signal_final), "gènes signal\n"
    )

    return(df)
}


# ============================================================
# Exemple d'utilisation
# ============================================================
# # Construire la banque UNE FOIS (voir rwr_normalization.R), a partir
# # de la distribution reelle du nombre de seeds par patient :
# observed_seed_counts <- sapply(seeds, length)
# permutation_bank <- build_permutation_bank(
#     RWR_net, seed_sizes = observed_seed_counts,
#     n_permutations = 2000, mc.cores = 4
# )
#
# # Patient unique :
# res_patient <- process_patient("Patient_1", HPE_results, HPE_seeds, permutation_bank)
#
# # Groupe :
# res_group <- process_group(sample(colnames(RWR_results), 20), RWR_results, seeds, permutation_bank)
