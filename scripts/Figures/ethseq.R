## ============================================================
## EthSEQ ancestry analysis from per-patient variant TSVs
## ============================================================

library(dplyr)
library(tidyverse)
library(EthSEQ)

## ------------------------------------------------------------
## 1. List all patient files (not loaded yet, just paths)
## ------------------------------------------------------------
setwd(paste(here::here(), "projects", sep = "/"))


patient_files_all <- list.files(path = "/home/jules/Documents/phd/Data/vcf_IOP/RESULTS/", full.names = TRUE)
length(patient_files_all)


## ------------------------------------------------------------
## 2. Build a single multi-sample VCF, reading one patient file
##    at a time (two passes) to keep RAM usage flat regardless
##    of cohort size.
##
##    Pass 1 : read each file once, keep only its distinct
##             (Chr, Start, Ref, Alt) site keys + genotype,
##             discard the rest, gc() before moving to the
##             next file.
##    Pass 2 : re-read each file once, build its genotype column
##             against the now-fixed site list, discard, gc().
##
##    NOTE: EthSEQ requires biallelic SNVs at fixed panel positions.
##    Absence of a variant in a patient's filtered table is encoded
##    here as 0/0 (assumed hom-ref) rather than true missing (./.),
##    since these files are a filtered candidate-variant set, not
##    a full joint-genotyped callset. Flag this as an approximation
##    if you don't have access to the pre-filter VCFs.
##
##    Genotype (0/1, 1/1) is read directly from column 16 of each
##    patient file (column name is sample-specific, hence positional
##    selection rather than a fixed name).
## ------------------------------------------------------------

clean_variants <- function(df) {
    df <- df %>% rename(GT_raw = 16)

    df %>%
        filter(!is.na(Chr), !is.na(Start), !is.na(Ref), !is.na(Alt)) %>%
        filter(
            nchar(Ref) == 1, nchar(Alt) == 1,
            Ref %in% c("A", "C", "G", "T"), Alt %in% c("A", "C", "G", "T")
        ) %>%
        mutate(GT_raw = case_when(
            GT_raw == "HET" ~ "0/1",
            GT_raw == "HOM" ~ "1/1",
            TRUE ~ NA_character_ # HOM? (uncertain) and anything else -> dropped
        )) %>%
        filter(!is.na(GT_raw)) %>%
        distinct(Chr, Start, Ref, Alt, .keep_all = TRUE)
}

normalize_chr <- function(chr_vec, chr_prefix = TRUE) {
    chr_vec <- as.character(chr_vec)
    if (chr_prefix) {
        ifelse(grepl("^chr", chr_vec), chr_vec, paste0("chr", chr_vec))
    } else {
        sub("^chr", "", chr_vec)
    }
}

build_ethseq_vcf <- function(patient_files, out_vcf, chr_prefix = TRUE, verbose = TRUE) {
    sample_ids <- tools::file_path_sans_ext(basename(patient_files))

    ## ---- Pass 1: collect the union of variant sites across patients ----
    if (verbose) message("Pass 1/2 : collecting variant sites (", length(patient_files), " files)...")

    site_list <- vector("list", length(patient_files))

    for (i in seq_along(patient_files)) {
        df <- read.table(patient_files[i], header = TRUE, sep = "\t", fill = TRUE, quote = "")
        df <- clean_variants(df)
        df$Chr <- normalize_chr(df$Chr, chr_prefix = chr_prefix)

        site_list[[i]] <- distinct(df, Chr, Start, Ref, Alt)

        rm(df)
        gc(verbose = FALSE)

        if (verbose && i %% 20 == 0) message("  Pass 1 : ", i, "/", length(patient_files))
    }

    all_sites <- bind_rows(site_list) %>%
        distinct(Chr, Start, .keep_all = TRUE) %>% # one variant per position (avoids duplicate-position conflicts)
        arrange(Chr, Start)

    rm(site_list)
    gc(verbose = FALSE)

    ## Unique per-site ID, required by EthSEQ/SNPRelate (VCF ID = "." for every
    ## row would otherwise collide and trigger "snp.id is not unique!")
    site_ids <- paste0(all_sites$Chr, "_", all_sites$Start, "_", all_sites$Ref, "_", all_sites$Alt)

    site_key <- paste(all_sites$Chr, all_sites$Start, all_sites$Ref, all_sites$Alt)

    ## ---- Pass 2: build the genotype matrix one patient column at a time ----
    if (verbose) message("Pass 2/2 : building genotype matrix (", nrow(all_sites), " sites x ", length(patient_files), " samples)...")

    gt_matrix <- matrix("0/0",
        nrow = nrow(all_sites), ncol = length(patient_files),
        dimnames = list(NULL, sample_ids)
    )

    for (i in seq_along(patient_files)) {
        df <- read.table(patient_files[i], header = TRUE, sep = "\t", fill = TRUE, quote = "")
        df <- clean_variants(df)
        df$Chr <- normalize_chr(df$Chr, chr_prefix = chr_prefix)

        pt_key <- paste(df$Chr, df$Start, df$Ref, df$Alt)
        match_idx <- match(pt_key, site_key) # position in all_sites for each variant row
        valid <- !is.na(match_idx)

        gt_matrix[match_idx[valid], i] <- df$GT_raw[valid] # real genotype (0/1 or 1/1)

        rm(df, pt_key, match_idx, valid)
        gc(verbose = FALSE)

        if (verbose && i %% 20 == 0) message("  Pass 2 : ", i, "/", length(patient_files))
    }

    ## ---- Write the VCF ----
    vcf_body <- cbind(
        data.frame(
            `#CHROM` = all_sites$Chr, POS = all_sites$Start, ID = site_ids,
            REF = all_sites$Ref, ALT = all_sites$Alt, QUAL = ".",
            FILTER = "PASS", INFO = ".", FORMAT = "GT", check.names = FALSE
        ),
        gt_matrix
    )

    con <- file(out_vcf, "w")
    writeLines(c(
        "##fileformat=VCFv4.2",
        paste0("#", paste(colnames(vcf_body), collapse = "\t"))
    ), con)
    close(con)
    write.table(vcf_body, out_vcf,
        append = TRUE, sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = FALSE
    )

    rm(vcf_body, gt_matrix)
    gc(verbose = FALSE)

    invisible(out_vcf)
}

## ------------------------------------------------------------
## 3. Build VCF and run EthSEQ
## ------------------------------------------------------------
getwd()

out_dir <- "ethseq_analysis/IOP_cohort/"
vcf_path <- file.path(out_dir, "cohort.vcf")
options(timeout = 3600) # 1 hour instead of the default 60s
# print working directory

build_ethseq_vcf(patient_files_all, vcf_path)

ethseq.Analysis(
    target.vcf = vcf_path,
    model.available = "Gencode.Exome", # ⚠️ to verify: model choice
    model.assembly = "hg19", # ⚠️ to verify: must match Start/Chr coordinate build
    model.pop = "All",
    out.dir = out_dir,
    verbose = TRUE,
    cores = 4,
    composite.model.call.rate = 0.6, # ⚠️ to verify: check actual per-SNP call rate before trusting this
    space = "2D"
)

## ------------------------------------------------------------
## 4. Load final report
## ------------------------------------------------------------

ethseq_report <- read.delim(file.path(out_dir, "Report.txt"), sep = "\t", as.is = TRUE, header = TRUE)
head(ethseq_report)
View(ethseq_report)
ethseq_report$pop %>%
    table() %>%
    as.data.frame()



HPE_report <- read.table("article/ethseq_analysis/HPE_cohort/Report.txt", header = TRUE)
IOP_report <- read.table("article/ethseq_analysis/IOP_cohort/Report.txt", header = TRUE)

table(HPE_report$pop)
14 / 199 * 100
185 / 199 * 100
table(IOP_report$pop)
18 / 298 * 100
280 / 298 * 100


tab <- rbind(
    HPE = table(HPE_report$pop),
    IOP = table(IOP_report$pop)
)
tab
#     AFR EUR
# HPE  14 185
# IOP  18 280

fisher.test(tab)


prop.test(c(14, 18), c(199, 298))$conf.int
