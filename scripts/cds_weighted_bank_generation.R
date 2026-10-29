# ============================================================
# Gene length (CDS) bias: weighted null, label permutation and
# seed frequency HPE vs IOP
#
# Review comment: Tier 1/2 candidates are very long genes (LRP2,
# LRP1, MACF1, MYCBP2, FAT3, CELSR1/2, TENM4, FBN2, KMT2C, CREBBP,
# DMD, COL4A1, ROBO2, DCC). Long genes carry more rare variants by
# chance, so they become seeds more often in ANY cohort. The current
# permutation bank draws random seeds uniformly, which ignores this.
#
# Sections:
#   1. CDS length per gene (UCSC refGene hg19, longest coding
#      transcript per gene). Downloaded once and cached.
#   2. Diagnostics: is seed frequency / z-score linked to CDS length?
#   3. Two new permutation banks built with the SAME seed sizes and
#      the SAME random draws of sizes:
#        - uniform  : seeds drawn uniformly (control, isolates the
#                     effect of the weighting from any other change)
#        - cds      : seeds drawn with probability proportional to
#                     CDS length
#   4. Signal genes for HPE and IOP under each bank + status of the
#      candidate genes.
#   5. Label permutation test (same procedure and same random splits
#      as Figure_3_label_permutation.R) with the CDS-weighted bank.
#   6. Seed carrier frequency per gene, HPE vs IOP (Fisher, BH).
# ============================================================

library(Matrix)
library(tidyverse)
library(parallel)
library(patchwork)
library(openxlsx)

setwd(paste(here::here(), "projects", sep = "/"))

source("article/scripts/custom_fct.R") # RWR(), gene_converter(), custom_theme()
source("article/scripts/RECENT_rwr_normalization.R") # normalize_rwr_score()
source("article/scripts/RECENT_patient_gene_pvalue.R") # process_profiles() / my_locfdr_robust()

dir.create("article/figures/Figure_3", recursive = TRUE, showWarnings = FALSE)
dir.create("article/tables/Figure_3", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
# IMPORTANT: must be the gamma used for the patients' RWR
# (compute_RWR_matrix() in rwr.R uses gamma = 0.9 by default, while
# build_permutation_bank() defaults to 0.7: check which value the
# current bank was built with).
rwr_gamma <- 0.9
rwr_tmax <- 1000
rwr_eps <- 1e-10

n_perm <- 1000 # label permutations (section 5)
n_cores <- max(1, detectCores() - 2)
rng_seed <- 42

candidates <- c(
    "LRP2", "LRP1", "MACF1", "MYCBP2", "FAT3", "CELSR1", "CELSR2", "TENM4",
    "FBN2", "KMT2C", "CREBBP", "DMD", "COL4A1", "ROBO2", "DCC"
)

col_observed <- "#2E86AB" # HPE color used in Figure_3.R
col_null <- "grey60"

write_single_xlsx <- function(df, path, sheet_name = "Sheet1") {
    wb <- createWorkbook()
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df, headerStyle = createStyle(textDecoration = "bold"))
    freezePane(wb, sheet_name, firstRow = TRUE)
    setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
    saveWorkbook(wb, path, overwrite = TRUE)
}

export_table <- function(df, name) {
    write_csv(df, paste0("article/tables/Figure_3/", name, ".csv"))
    write_single_xlsx(df, paste0("article/tables/Figure_3/", name, ".xlsx"), substr(name, 1, 31))
}

# ------------------------------------------------------------
# run_with_progress()
# Runs fun(i) for i in 1..n, in chunks. Each chunk is parallelised
# with mclapply, and progress (done / total, elapsed, ETA) is printed
# from the main process after each chunk, so it shows in RStudio too
# (messages from mclapply workers do not).
# ------------------------------------------------------------
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

load("article/transit_data/RWR_results_case.RData")
HPE_results <- as.matrix(RWR_results_z)
rm(RWR_results_z)
load("article/transit_data/iop_RWR_results_case.RData")
IOP_results <- as.matrix(RWR_results_z)
rm(RWR_results_z)

load("article/transit_data/patient_variants.RData")
HPE_patient_variants <- patient_variants
rm(patient_variants)
load("article/transit_data/iop_patient_variants.RData")
IOP_patient_variants <- patient_variants
rm(patient_variants)

load("article/transit_data/permutation_bank.RData") # -> `bank` (current, uniform)
bank_original <- bank
rm(bank)

coex_genes <- rownames(data_3p$coex_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

HPE_seeds <- lapply(HPE_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(HPE_seeds) <- colnames(HPE_results)

IOP_seeds <- lapply(IOP_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(IOP_seeds) <- colnames(IOP_results)

# Network (same construction as rwr.R)
RWR_net <- data_3p$coex_case
rownames(RWR_net) <- colnames(RWR_net) <- coex_genes
network_genes <- rownames(RWR_net)

# Prefix patient IDs by cohort so they are unique once pooled
colnames(HPE_results) <- names(HPE_seeds) <- paste0("HPE_", colnames(HPE_results))
colnames(IOP_results) <- names(IOP_seeds) <- paste0("IOP_", colnames(IOP_results))

stopifnot(setequal(rownames(HPE_results), rownames(IOP_results)))
IOP_results <- IOP_results[rownames(HPE_results), , drop = FALSE]
stopifnot(setequal(rownames(HPE_results), network_genes))
stopifnot(setequal(rownames(bank_original), network_genes))

# Valid patients: same rule as process_profiles()
is_valid_patient <- function(RWR_results, seeds) {
    vapply(colnames(RWR_results), function(pid) {
        scores <- RWR_results[, pid]
        is_seed <- rownames(RWR_results) %in% seeds[[pid]]
        scores[is_seed] <- median(scores[!is_seed], na.rm = TRUE)
        mean(is.na(scores)) <= 0.5
    }, logical(1))
}

HPE_ids <- colnames(HPE_results)[is_valid_patient(HPE_results, HPE_seeds)]
IOP_ids <- colnames(IOP_results)[is_valid_patient(IOP_results, IOP_seeds)]

pool_results <- cbind(HPE_results[, HPE_ids, drop = FALSE], IOP_results[, IOP_ids, drop = FALSE])
pool_seeds <- c(HPE_seeds[HPE_ids], IOP_seeds[IOP_ids])
pool_ids <- colnames(pool_results)

n_HPE <- length(HPE_ids)
n_IOP <- length(IOP_ids)
cat("Pool: ", length(pool_ids), " patients (", n_HPE, " HPE + ", n_IOP, " IOP)\n", sep = "")


# ============================================================
# 1. CDS length per gene, hg19 (UCSC refGene)
# ------------------------------------------------------------
# Same annotation as the variant calls (Gene.refGeneWithVer, hg19).
# CDS length of a transcript = sum of the exon parts lying between
# cdsStart and cdsEnd (stop codon included). Gene CDS length = the
# longest coding (NM_) transcript on a main chromosome.
# ============================================================

refgene_path <- "article/transit_data/refGene_hg19.txt.gz"
cds_path <- "article/transit_data/cds_length_hg19.csv"

if (!file.exists(cds_path)) {
    if (!file.exists(refgene_path)) {
        download.file(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz",
            refgene_path,
            mode = "wb"
        )
    }

    refgene <- read_tsv(
        refgene_path,
        col_names = c(
            "bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
            "exonCount", "exonStarts", "exonEnds", "score", "name2",
            "cdsStartStat", "cdsEndStat", "exonFrames"
        ),
        col_types = cols(.default = col_character())
    ) %>%
        filter(
            str_starts(name, "NM_"),
            str_detect(chrom, "^chr([0-9]+|X|Y)$")
        ) %>%
        mutate(across(c(cdsStart, cdsEnd), as.numeric))

    tx_cds_length <- function(cds_start, cds_end, ex_starts, ex_ends) {
        if (cds_end <= cds_start) {
            return(0)
        }
        s <- as.numeric(strsplit(ex_starts, ",")[[1]])
        e <- as.numeric(strsplit(ex_ends, ",")[[1]])
        sum(pmax(0, pmin(e, cds_end) - pmax(s, cds_start)))
    }

    refgene$cds_length <- mapply(
        tx_cds_length,
        refgene$cdsStart, refgene$cdsEnd, refgene$exonStarts, refgene$exonEnds
    )

    cds_table <- refgene %>%
        filter(cds_length > 0) %>%
        group_by(gene = name2) %>%
        arrange(dplyr::desc(cds_length), .by_group = TRUE) %>%
        summarise(
            cds_length = dplyr::first(cds_length),
            transcript = dplyr::first(name),
            chrom = dplyr::first(chrom),
            .groups = "drop"
        )

    write_csv(cds_table, cds_path)
}

cds_table <- read_csv(cds_path, show_col_types = FALSE)

# Map to network genes; genes missing from refGene get the median length
cds_network <- tibble(gene = network_genes) %>%
    left_join(cds_table, by = "gene") %>%
    mutate(
        cds_missing = is.na(cds_length),
        cds_length = if_else(cds_missing, median(cds_length, na.rm = TRUE), cds_length),
        sex_chrom = chrom %in% c("chrX", "chrY")
    )

cat(
    "CDS length found for ", sum(!cds_network$cds_missing), "/", nrow(cds_network),
    " network genes (", sum(cds_network$cds_missing), " imputed to the median)\n",
    sep = ""
)
export_table(cds_network, "cds_length_network_genes")

cds_weights <- setNames(cds_network$cds_length, cds_network$gene)


# ============================================================
# 3. New permutation banks: uniform (control) and CDS-weighted
# ------------------------------------------------------------
# Same number of columns as the current bank. Seed-set sizes are
# drawn from the sizes of the current bank (attr "seed_sizes") when
# available, otherwise from the observed seed counts of all patients.
# The sizes are drawn ONCE and shared by both banks, so the only
# difference between them is how genes are picked.
# ============================================================

n_bank <- ncol(bank_original)
size_pool <- attr(bank_original, "seed_sizes")
if (is.null(size_pool)) {
    size_pool <- lengths(pool_seeds)
    size_pool <- size_pool[size_pool > 0]
}

# Positional weights (same order as network_genes), safe even if some
# network gene names are NA after ID conversion
cds_prob <- cds_network$cds_length
stopifnot(length(cds_prob) == length(network_genes), all(is.finite(cds_prob)))

set.seed(rng_seed)
bank_sizes <- size_pool[sample.int(length(size_pool), n_bank, replace = TRUE)]
seed_sets_uniform <- lapply(bank_sizes, function(s) sample(network_genes, s))
seed_sets_cds <- lapply(bank_sizes, function(s) network_genes[sample.int(length(network_genes), s, prob = cds_prob)])

build_bank_from_seed_sets <- function(seed_sets, label) {
    run_one <- function(i) {
        RWR(RWR_net, seed_sets[[i]], gamma = rwr_gamma, tmax = rwr_tmax, eps = rwr_eps, norm = TRUE)$p %>% unname()
    }
    scores <- run_with_progress(length(seed_sets), run_one, label = label, n_cores = n_cores)
    failed <- vapply(scores, function(x) inherits(x, "try-error"), logical(1))
    if (any(failed)) stop(sum(failed), " RWR run(s) failed while building the bank '", label, "'")
    b <- do.call(cbind, scores)
    rownames(b) <- network_genes
    colnames(b) <- paste0("perm_", seq_along(seed_sets))
    attr(b, "seed_sizes") <- lengths(seed_sets)
    b
}

source("article/scripts/RECENT_rwr_normalization.R") # normalize_rwr_score()

bank_uniform_path <- "article/transit_data/permutation_bank_uniform_.RData"
bank_cds_path <- "article/transit_data/permutation_bank_cds_weighted_.RData"

if (file.exists(bank_uniform_path)) {
    load(bank_uniform_path)
} else {
    bank_uniform <- build_bank_from_seed_sets(seed_sets_uniform, label = "Bank uniform")
    save(bank_uniform, file = bank_uniform_path)
}

if (file.exists(bank_cds_path)) {
    load(bank_cds_path)
} else {
    bank_cds <- build_bank_from_seed_sets(seed_sets_cds, label = "Bank CDS-weighted")
    save(bank_cds, file = bank_cds_path)
}
