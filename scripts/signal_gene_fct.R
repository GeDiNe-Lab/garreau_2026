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

    # --- Safeguards: not enough points, or degenerate distribution ---
    if (length(z) < 30) {
        message(
            "my_locfdr_robust: fewer than 30 finite values (", length(z),
            "), using N(0,1) without robust scaling"
        )
        log_pval <- pnorm(z, mean = 0, sd = 1, lower.tail = FALSE, log.p = TRUE)
        names(log_pval) <- names(z)
        return(log_pval)
    }

    mad_z <- mad(z)
    if (!is.finite(mad_z) || mad_z == 0) {
        message(
            "my_locfdr_robust: MAD is zero or non-finite (degenerate distribution, ",
            "too many identical values) - fallback to N(0,1) on raw z"
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

    # --- Additional safeguard: zero bin_width (min == max after rounding) ---
    if (!is.finite(bin_width) || bin_width == 0) {
        message("my_locfdr_robust: bin_width is zero - fallback to N(0,1) on raw z")
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
        message("Not enough valid bins to estimate the null - N(0,1) used")
        mu0 <- 0
        sigma0 <- 1
    } else {
        fit_null <- lm(log(f_hat[valid] + 1e-10) ~ mids[valid] + I(mids[valid]^2))
        cc <- coef(fit_null)[3]
        b <- coef(fit_null)[2]

        if (is.na(cc) || cc >= 0) {
            message(
                "my_locfdr_robust: invalid quadratic coefficient (c = ", round(cc, 4),
                ", not concave or NA) - fallback to N(0,1) on robust z"
            )
            mu0 <- 0
            sigma0 <- 1
        } else {
            sigma0 <- sqrt(-1 / (2 * cc))
            mu0 <- -b / (2 * cc)
        }
    }

    log_pval <- pnorm(z_rob, mean = mu0, sd = sigma0, lower.tail = FALSE, log.p = TRUE)
    message(
        "my_locfdr_robust: empirical null estimated successfully (mu0 = ", round(mu0, 4),
        ", sigma0 = ", round(sigma0, 4), ", n = ", N, ")"
    )
    return(log_pval)
}


# ============================================================
# process_profiles()
# ------------------------------------------------------------
# Replaces process_patient() and process_group(). Same signature
# as the former process_group() so existing call sites do not
# break. Handles a single profile or a group of profiles in a
# unified way:
#   - patient_ids of length 1  -> equivalent to the former
#     process_patient() (but with seed imputation to the median,
#     no exclusion of seeds)
#   - patient_ids of length > 1 -> equivalent to the former
#     process_group()
#
# In both cases: seed scores are imputed to the median of the
# non-seed scores (per patient, on raw scores), then the genes x
# N patients matrix (N >= 1) is passed to normalize_rwr_score(),
# which averages the patients AND builds a null matched on N by
# resampling the precomputed permutation bank (no RWR recomputation,
# N = 1 included).
#
# The FDR computation logic (my_locfdr_robust, BH, quantile
# threshold) is unchanged from the former process_group().
# ============================================================
process_profiles <- function(patient_ids = NULL, RWR_results, seeds, permutation_bank) {
    # Default: all patients
    if (is.null(patient_ids)) patient_ids <- colnames(RWR_results)

    # Check
    invalid <- patient_ids[!patient_ids %in% colnames(RWR_results)]
    if (length(invalid) > 0) warning("Patients not found: ", paste(invalid, collapse = ", "))
    patient_ids <- patient_ids[patient_ids %in% colnames(RWR_results)]

    if (length(patient_ids) == 0) {
        message("No valid patient in the selection")
        return(NULL)
    }

    # For each patient: replace seed scores with the median of non-seed scores (raw scores)
    # sapply(simplify = FALSE) + do.call(cbind, ...) guarantees a matrix
    # even when patient_ids has length 1 (plain sapply would return a vector)
    scores_list <- sapply(patient_ids, function(pid) {
        scores <- RWR_results[, pid]
        seed_genes <- seeds[[pid]]
        is_seed <- rownames(RWR_results) %in% seed_genes

        # Median computed without the seeds
        median_nonseed <- median(scores[!is_seed], na.rm = TRUE)

        # Replace seed scores with the median
        scores[is_seed] <- median_nonseed

        return(scores)
    }, simplify = FALSE)

    scores_imputed <- do.call(cbind, scores_list)
    rownames(scores_imputed) <- rownames(RWR_results)
    colnames(scores_imputed) <- patient_ids

    # Filter patients with too many NAs
    na_prop <- colMeans(is.na(scores_imputed))
    excluded <- patient_ids[na_prop > 0.5]
    if (length(excluded) > 0) message("Patients excluded (>50% NA): ", paste(excluded, collapse = ", "))
    scores_imputed <- scores_imputed[, na_prop <= 0.5, drop = FALSE]

    if (ncol(scores_imputed) == 0) {
        message("No valid patient in this group")
        return(NULL)
    }

    # Mean raw score per gene (= patient's raw score if N = 1), kept for reporting/comparison
    mean_scores_raw <- rowMeans(scores_imputed, na.rm = TRUE)

    # --- Z-score normalization via the permutation bank ---
    # scores_imputed (genes x N valid patients) is passed directly:
    # normalize_rwr_score averages the patients AND builds the null by
    # aggregating N profiles drawn from the bank (N = ncol(scores_imputed),
    # N = 1 included), so the null exactly matches the tested statistic.

    z_scores <- normalize_rwr_score(
        scores_imputed,
        permutation_bank = permutation_bank,
        verbose = FALSE
    )
    
    
    z_scores <- z_scores[is.finite(z_scores)]
    # locfdr + correction pipeline, on the z-scores (single patient or group)
    log_pvalues <- my_locfdr_robust(z_scores)
    pvalues_safe <- exp(pmax(log_pvalues, -745))
    pvalues_corrected <- p.adjust(pvalues_safe, method = "BH")
    threshold_75 <- quantile(z_scores, 0.75, na.rm = TRUE)

    df <- data.frame(
        gene           = names(z_scores),
        mean_score_raw = mean_scores_raw[names(z_scores)],
        z_score        = z_scores,
        n_patients     = ncol(scores_imputed),
        log_pval       = log_pvalues,
        pval_bh        = pvalues_corrected,
        signal_locfdr  = pvalues_corrected < 0.05,
        signal_75      = z_scores >= threshold_75,
        signal_final   = pvalues_corrected < 0.05 & z_scores >= threshold_75
    )

    cat(
        "Profile(s):", length(patient_ids), "patient(s) requested |",
        ncol(scores_imputed), "patient(s) valid |",
        sum(df$signal_final), "signal genes\n"
    )

    return(df)
}


# ============================================================
# Usage example
# ============================================================
# # Build the bank ONCE (see rwr_normalization.R), from the observed
# # distribution of the number of seeds per patient:
# observed_seed_counts <- sapply(seeds, length)
# permutation_bank <- build_permutation_bank(
#     RWR_net, seed_sizes = observed_seed_counts,
#     n_permutations = 2000, mc.cores = 4
# )
#
# # Single patient (N = 1, same logic as the group) :
# res_patient <- process_profiles("Patient_1", HPE_results, HPE_seeds, permutation_bank)
#
# # Group :
# res_group <- process_profiles(sample(colnames(RWR_results), 20), RWR_results, seeds, permutation_bank)
