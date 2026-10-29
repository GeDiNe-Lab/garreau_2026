# ============================================================
# Cohort label permutation test (HPE vs IOP): number of signal genes
#
# Question: is the gap in signal gene counts between HPE and IOP
# specific to the HPE grouping, or would any random split of the
# same patients into two groups of the same sizes give a similar gap?
#
# Procedure:
#   1. Pool all valid HPE and IOP patients. A patient is valid if it
#      keeps at least one seed after network filtering AND passes the
#      NA rule of process_profiles(). Per-patient RWR profiles and
#      seeds are kept as they are, nothing is recomputed.
#   2. Run the real analysis once (true labels) as a sanity check:
#      the counts must match Figure_3.R (run on the same patients).
#   3. B times: shuffle the labels, split the pool into a group of
#      n_HPE patients (pseudo-HPE) and a group of n_IOP patients
#      (pseudo-IOP), run process_profiles() on each group and count
#      the signal genes (signal_final).
#   4. Empirical one-sided p-values for:
#        - n_A    : signal genes in the group of size n_HPE (upper tail)
#        - n_B    : signal genes in the group of size n_IOP (lower tail)
#        - n_diff : n_A - n_B (upper tail)
#
# The group-level pipeline is exactly the one used in Figure_3.R
# (process_profiles(): seed imputation, mean, normalisation against
# the permutation bank with /sqrt(N), my_locfdr_robust, BH, quantile
# threshold). No RWR is recomputed.
#
# Signal gene lists of every permutation are saved, so the HPO
# comparison can be done later without rerunning this script.
# ============================================================

library(Matrix)
library(tidyverse)
library(parallel)
library(patchwork)
library(openxlsx)

setwd(paste(here::here(), "projects", sep = "/"))

source("article/scripts/custom_fct.R")
source("article/scripts/RECENT_rwr_normalization.R") # normalize_rwr_score()
source("article/scripts/RECENT_patient_gene_pvalue.R") # process_profiles() / my_locfdr_robust()

dir.create("article/figures/Figure_3", recursive = TRUE, showWarnings = FALSE)
dir.create("article/tables/Figure_3", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
n_perm <- 1000
n_cores <- max(1, detectCores() - 2)
rng_seed <- 42
make_figure <- TRUE # supplementary figure (not used in the article for now)

col_observed <- "#2E86AB" # HPE color used in Figure_3.R
col_null <- "grey60"

# Same export helper as Figure_3.R
write_single_xlsx <- function(df, path, sheet_name = "Sheet1") {
    wb <- createWorkbook()
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df, headerStyle = createStyle(textDecoration = "bold"))
    freezePane(wb, sheet_name, firstRow = TRUE)
    setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
    saveWorkbook(wb, path, overwrite = TRUE)
}

# Runs fun(i) for i in 1..n in chunks (mclapply inside each chunk) and
# prints progress from the main process after each chunk.
format_duration <- function(secs) {
    secs <- round(secs)
    sprintf("%dh%02dm%02ds", secs %/% 3600, (secs %% 3600) %/% 60, secs %% 60)
}

run_with_progress <- function(n, fun, label, n_cores, chunk_size = max(10, n_cores * 5)) {
    chunks <- split(seq_len(n), ceiling(seq_len(n) / chunk_size))
    out <- vector("list", n)
    t0 <- Sys.time()
    message("[", label, "] start: ", n, " iterations on ", n_cores, " core(s)")
    for (ch in chunks) {
        res <- if (n_cores > 1) mclapply(ch, fun, mc.cores = n_cores) else lapply(ch, fun)
        out[ch] <- res
        done <- max(ch)
        elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        eta <- elapsed / done * (n - done)
        message(sprintf(
            "[%s] %d/%d (%.0f%%) | elapsed %s | ETA %s",
            label, done, n, 100 * done / n, format_duration(elapsed), format_duration(eta)
        ))
    }
    out
}


# ============================================================
# 0. Data loading (same objects as Figure_3.R)
# ============================================================

load("article/transit_data/data_3p.RData")

load("article/transit_data/RWR/HPE_RWR_results_.RData") # -> HPE_results
load("article/transit_data/RWR/IOP_RWR_results_.RData") # -> IOP_results

load("article/transit_data/patient_variants.RData")
HPE_patient_variants <- patient_variants
rm(patient_variants)
load("article/transit_data/iop_patient_variants.RData")
IOP_patient_variants <- patient_variants
rm(patient_variants)

load("article/transit_data/permutation_bank_cds_weighted.RData")
bank <- bank_blend

coex_genes <- rownames(data_3p$coex_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

HPE_seeds <- lapply(HPE_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(HPE_seeds) <- colnames(HPE_results)

IOP_seeds <- lapply(IOP_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(IOP_seeds) <- colnames(IOP_results)


# ============================================================
# 1. Pool the two cohorts
# ============================================================

# Both matrices must share the same genes (same network)
stopifnot(setequal(rownames(HPE_results), rownames(IOP_results)))
IOP_results <- IOP_results[rownames(HPE_results), , drop = FALSE]

# Prefix patient IDs by cohort so they are unique once pooled
colnames(HPE_results) <- names(HPE_seeds) <- paste0("HPE_", colnames(HPE_results))
colnames(IOP_results) <- names(IOP_seeds) <- paste0("IOP_", colnames(IOP_results))

# Rule 1: at least one seed left after network filtering
has_network_seed <- function(RWR_results, seeds) {
    vapply(colnames(RWR_results), function(pid) {
        length(intersect(unique(seeds[[pid]]), rownames(RWR_results))) > 0
    }, logical(1))
}

# Rule 2: same NA rule as process_profiles() (seed imputation, then
# patients with > 50% NA are excluded)
is_valid_patient <- function(RWR_results, seeds) {
    vapply(colnames(RWR_results), function(pid) {
        scores <- RWR_results[, pid]
        is_seed <- rownames(RWR_results) %in% seeds[[pid]]
        scores[is_seed] <- median(scores[!is_seed], na.rm = TRUE)
        mean(is.na(scores)) <= 0.5
    }, logical(1))
}

HPE_has_seed <- has_network_seed(HPE_results, HPE_seeds)
IOP_has_seed <- has_network_seed(IOP_results, IOP_seeds)
HPE_valid_na <- is_valid_patient(HPE_results, HPE_seeds)
IOP_valid_na <- is_valid_patient(IOP_results, IOP_seeds)

HPE_ids <- colnames(HPE_results)[HPE_has_seed & HPE_valid_na]
IOP_ids <- colnames(IOP_results)[IOP_has_seed & IOP_valid_na]

exclusion_df <- tibble(
    cohort = c("HPE", "IOP"),
    n_total = c(ncol(HPE_results), ncol(IOP_results)),
    excluded_no_network_seed = c(sum(!HPE_has_seed), sum(!IOP_has_seed)),
    excluded_na_rule_only = c(sum(HPE_has_seed & !HPE_valid_na), sum(IOP_has_seed & !IOP_valid_na)),
    n_kept = c(length(HPE_ids), length(IOP_ids))
)
print(exclusion_df)

pool_results <- cbind(HPE_results[, HPE_ids, drop = FALSE], IOP_results[, IOP_ids, drop = FALSE])
pool_seeds <- c(HPE_seeds[HPE_ids], IOP_seeds[IOP_ids])
pool_ids <- colnames(pool_results)

n_HPE <- length(HPE_ids)
n_IOP <- length(IOP_ids)
n_pool <- length(pool_ids)

cat("Pool: ", n_pool, " patients (", n_HPE, " HPE + ", n_IOP, " IOP)\n", sep = "")


# ============================================================
# 2. Group-level signal genes (wraps process_profiles, silent)
# ============================================================

get_group_signal <- function(patient_ids) {
    res <- NULL
    invisible(capture.output(
        suppressMessages(
            res <- process_profiles(
                patient_ids = patient_ids,
                RWR_results = pool_results,
                seeds = pool_seeds,
                permutation_bank = bank
            )
        )
    ))
    if (is.null(res)) {
        return(character(0))
    }
    res %>%
        filter(signal_final) %>%
        pull(gene)
}

run_split <- function(ids_A, ids_B) {
    sig_A <- get_group_signal(ids_A)
    sig_B <- get_group_signal(ids_B)
    list(
        n_A = length(sig_A),
        n_B = length(sig_B),
        n_true_HPE_in_A = sum(ids_A %in% HPE_ids),
        signal_A = sig_A,
        signal_B = sig_B
    )
}


# ============================================================
# 3. Observed values (true labels) = sanity check vs Figure_3.R
# ============================================================

observed <- run_split(HPE_ids, IOP_ids)
obs_n_A <- observed$n_A
obs_n_B <- observed$n_B
obs_diff <- obs_n_A - obs_n_B

cat("Observed: HPE = ", obs_n_A, " signal genes | IOP = ", obs_n_B,
    " signal genes | difference = ", obs_diff, "\n",
    sep = ""
)


# ============================================================
# 4. Label permutations
# ------------------------------------------------------------
# All splits are drawn up front in the main process, so results are
# reproducible whatever the number of cores.
# ============================================================

set.seed(rng_seed)
perm_splits <- lapply(seq_len(n_perm), function(i) sample(pool_ids, n_HPE))

run_one_perm <- function(i) {
    ids_A <- perm_splits[[i]]
    ids_B <- setdiff(pool_ids, ids_A)
    out <- run_split(ids_A, ids_B)
    out$perm <- i
    out
}

t0 <- Sys.time()
perm_results <- run_with_progress(n_perm, run_one_perm, label = "Label permutation", n_cores = n_cores)
cat("Permutations done in", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")

# mclapply returns try-error objects on worker failure: report and drop them
failed <- vapply(perm_results, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) {
    warning(sum(failed), " permutation(s) failed and were dropped")
    perm_results <- perm_results[!failed]
}

null_df <- tibble(
    perm            = vapply(perm_results, `[[`, numeric(1), "perm"),
    n_true_HPE_in_A = vapply(perm_results, `[[`, numeric(1), "n_true_HPE_in_A"),
    n_A             = vapply(perm_results, `[[`, numeric(1), "n_A"),
    n_B             = vapply(perm_results, `[[`, numeric(1), "n_B")
) %>%
    mutate(n_diff = n_A - n_B)

B <- nrow(null_df)


# ============================================================
# 5. Empirical p-values (one-sided, +1 correction)
# ============================================================

emp_p_upper <- function(null, obs) (1 + sum(null >= obs)) / (length(null) + 1)
emp_p_lower <- function(null, obs) (1 + sum(null <= obs)) / (length(null) + 1)

summary_df <- tibble(
    statistic = c("n_A (signal genes, group of size n_HPE)", "n_B (signal genes, group of size n_IOP)", "n_A - n_B"),
    observed = c(obs_n_A, obs_n_B, obs_diff),
    null_mean = c(mean(null_df$n_A), mean(null_df$n_B), mean(null_df$n_diff)),
    null_median = c(median(null_df$n_A), median(null_df$n_B), median(null_df$n_diff)),
    null_q95 = c(quantile(null_df$n_A, 0.95), quantile(null_df$n_B, 0.95), quantile(null_df$n_diff, 0.95)),
    null_max = c(max(null_df$n_A), max(null_df$n_B), max(null_df$n_diff)),
    empirical_p = c(
        emp_p_upper(null_df$n_A, obs_n_A),
        emp_p_lower(null_df$n_B, obs_n_B), # IOP: is it LOWER than expected?
        emp_p_upper(null_df$n_diff, obs_diff)
    ),
    tail = c("upper", "lower", "upper"),
    n_permutations = B,
    n_HPE = n_HPE,
    n_IOP = n_IOP
)

print(summary_df)
