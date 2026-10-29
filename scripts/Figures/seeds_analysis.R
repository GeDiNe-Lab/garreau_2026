# Loading packages and functions
library(Matrix)
library(tidyverse)
library(parallel)
library(tibble)
library(igraph)
library(ggplot2)
library(cluster)
library(stats)
library(circlize)
library(plotly)
library(WebGestaltR)
library(scales)
library(viridis)

# Setting working directory
setwd(paste(here::here(), "projects", sep = "/"))

# Loading custom functions
source("article/scripts/custom_fct.R")
source("article/scripts/RECENT_webgestalt_ORA.R")

load("article/transit_data/data_3p.RData")

load("article/transit_data/patient_variants.RData")
HPE_variants <- patient_variants
rm(patient_variants)
load("article/transit_data/iop_patient_variants.RData")
IOP_variants <- patient_variants
rm(patient_variants)

# Filter patient variants to keep only those falling in genes present in the co-expression network
coex_genes <- rownames(data_3p$coex_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

# ============================================================
# Seed count per patient, before network filtering
# ============================================================
seed_count_before <- bind_rows(
    tibble(
        patient = seq_along(HPE_variants),
        n_seeds = map_int(HPE_variants, ~ n_distinct(.x$Gene.refGeneWithVer)),
        cohort = "HPE"
    ),
    tibble(
        patient = seq_along(IOP_variants),
        n_seeds = map_int(IOP_variants, ~ n_distinct(.x$Gene.refGeneWithVer)),
        cohort = "IOP"
    )
) %>%
    mutate(filter_status = "Before")

HPE_variants <- lapply(HPE_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)
})
IOP_variants <- lapply(IOP_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)
})

# ============================================================
# Seed count per patient, after network filtering
# ============================================================
seed_count_after <- bind_rows(
    tibble(
        patient = seq_along(HPE_variants),
        n_seeds = map_int(HPE_variants, ~ n_distinct(.x$Gene.refGeneWithVer)),
        cohort = "HPE"
    ),
    tibble(
        patient = seq_along(IOP_variants),
        n_seeds = map_int(IOP_variants, ~ n_distinct(.x$Gene.refGeneWithVer)),
        cohort = "IOP"
    )
) %>%
    mutate(filter_status = "After")

seed_count_summary <- bind_rows(seed_count_before, seed_count_after) %>%
    mutate(filter_status = factor(filter_status, levels = c("Before", "After")))

p_seed_distribution <- ggplot(seed_count_summary, aes(x = filter_status, y = n_seeds, fill = filter_status)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    facet_wrap(~cohort) +
    scale_fill_viridis_d() +
    labs(
        x = NULL,
        y = "Number of seeds"
    ) +
    custom_theme() +
    theme(
        legend.position = "none",
        strip.text = element_text(size = 20, face = "bold"), # titres des panels (HPE / IOP)
        axis.text.x = element_text(size = 18), # Before / After
        axis.text.y = element_text(size = 18) # nombres de l'axe Y
    )
ggsave("article/figures/Figure_2/seed_count_distribution.svg", plot = p_seed_distribution, width = 9, height = 6)
# ============================================================
# recode_func_refgene()
# ------------------------------------------------------------
# Groups every Func.refGene value that is NOT in the `keep`
# list under the "Other" label.
# ============================================================
recode_pred_effect <- function(x) {
    keep <- c(
        "frameshift deletion", "frameshift insertion",
        "nonframeshift deletion", "nonframeshift insertion", "nonframeshift substitution",
        "startloss", "stopgain", "stoploss"
    )
    ifelse(x %in% keep, x, "Other")
}

# ============================================================
# HPE
# ============================================================

#  Variant genomic location table (Func.refGene)
HPE_variants_function <- lapply(HPE_variants, function(df) {
    return(df$Func.refGene)
}) %>%
    unlist() %>%
    table() %>%
    as.data.frame()

#  ClinVar significance table
HPE_CLNSIG <- lapply(HPE_variants, function(df) {
    return(df$CLNSIG)
}) %>%
    unlist() %>%
    table() %>%
    as.data.frame()
colnames(HPE_CLNSIG)[1] <- "category"
HPE_CLNSIG <- HPE_CLNSIG %>%
    mutate(category = ifelse(category == ".", "Not reported", as.character(category)))

#  Predicted effect table
HPE_PRED_effect <- lapply(HPE_variants, function(df) {
    recode_pred_effect(df$PRED_effect)
}) %>%
    unlist() %>%
    table() %>%
    as.data.frame()

#  Union of all seeds
HPE_seeds <- HPE_variants %>%
    lapply(function(df) {
        df$Gene.refGeneWithVer
    }) %>%
    unlist() %>%
    unique()

#  Enrichment on seeds
HPE_seeds_GO <- run_ora_gobp(HPE_seeds)
HPE_seeds_HPO <- run_ora_hpo(HPE_seeds)
HPE_seeds_Reactome <- run_ora_reactome(HPE_seeds)

# ============================================================
# IOP
# ============================================================

IOP_variants <- lapply(IOP_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)
})

#  Variant genomic location table (Func.refGene)
IOP_variants_function <- lapply(IOP_variants, function(df) {
    return(df$Func.refGene)
}) %>%
    unlist() %>%
    table() %>%
    as.data.frame()

#  ClinVar significance table
IOP_CLNSIG <- lapply(IOP_variants, function(df) {
    return(df$CLNSIG)
}) %>%
    unlist() %>%
    table() %>%
    as.data.frame()
colnames(IOP_CLNSIG)[1] <- "category"
IOP_CLNSIG <- IOP_CLNSIG %>%
    mutate(category = ifelse(category == ".", "Not reported", as.character(category)))

#  Predicted effect table
IOP_PRED_effect <- lapply(IOP_variants, function(df) {
    recode_pred_effect(df$PRED_effect)
}) %>%
    unlist() %>%
    table() %>%
    as.data.frame()

#  Union of all seeds
IOP_seeds <- IOP_variants %>%
    lapply(function(df) {
        df$Gene.refGeneWithVer
    }) %>%
    unlist() %>%
    unique()

#  Enrichment on seeds
IOP_seeds_GO <- run_ora_gobp(IOP_seeds)
IOP_seeds_HPO <- run_ora_hpo(IOP_seeds)
IOP_seeds_Reactome <- run_ora_reactome(IOP_seeds)


p_HPE_go <- plot_ora_dotplot(HPE_seeds_GO, n = 20, title = "HPE - GO-BP")
p_HPE_hpo <- plot_ora_dotplot(HPE_seeds_HPO, n = 20, title = "HPE - HPO")
p_HPE_reactome <- plot_ora_dotplot(HPE_seeds_Reactome, n = 20, title = "HPE - Reactome")


p_IOP_go <- plot_ora_dotplot(IOP_seeds_GO, n = 20, title = "IOP - GO-BP")
p_IOP_hpo <- plot_ora_dotplot(IOP_seeds_HPO, n = 20, title = "IOP - HPO")
p_IOP_reactome <- plot_ora_dotplot(IOP_seeds_Reactome, n = 20, title = "IOP - Reactome")

# ============================================================
# Saving plots
# ============================================================

ggsave("article/figures/Figure_2/HPE_seeds_GO_dotplot.svg", plot = p_HPE_go, width = 7, height = 4)
ggsave("article/figures/Figure_2/HPE_seeds_HPO_dotplot.svg", plot = p_HPE_hpo, width = 7, height = 4)
ggsave("article/figures/Figure_2/HPE_seeds_Reactome_dotplot.svg", plot = p_HPE_reactome, width = 7, height = 4)

ggsave("article/figures/Figure_2/IOP_seeds_GO_dotplot.svg", plot = p_IOP_go, width = 7, height = 4)
ggsave("article/figures/Figure_2/IOP_seeds_HPO_dotplot.svg", plot = p_IOP_hpo, width = 7, height = 4)
ggsave("article/figures/Figure_2/IOP_seeds_Reactome_dotplot.svg", plot = p_IOP_reactome, width = 7, height = 4)
# ============================================================
# Summary tables (readable by a biologist)
# ============================================================

seed_count_per_patient <- seed_count_summary %>%
    arrange(cohort, filter_status, patient)

HPE_variant_function_summary <- HPE_variants_function
colnames(HPE_variant_function_summary)[1] <- "category"
HPE_variant_function_summary <- HPE_variant_function_summary %>%
    arrange(desc(Freq))

IOP_variant_function_summary <- IOP_variants_function
colnames(IOP_variant_function_summary)[1] <- "category"
IOP_variant_function_summary <- IOP_variant_function_summary %>%
    arrange(desc(Freq))

HPE_CLNSIG_summary <- HPE_CLNSIG %>%
    arrange(desc(Freq))

IOP_CLNSIG_summary <- IOP_CLNSIG %>%
    arrange(desc(Freq))

HPE_PRED_effect_summary <- HPE_PRED_effect
colnames(HPE_PRED_effect_summary)[1] <- "category"
HPE_PRED_effect_summary <- HPE_PRED_effect_summary %>%
    arrange(desc(Freq))

IOP_PRED_effect_summary <- IOP_PRED_effect
colnames(IOP_PRED_effect_summary)[1] <- "category"
IOP_PRED_effect_summary <- IOP_PRED_effect_summary %>%
    arrange(desc(Freq))

HPE_seeds_GO_summary <- HPE_seeds_GO %>%
    dplyr::select(geneSet, description, size, overlap, enrichmentRatio, pValue, FDR, userId) %>%
    arrange(-enrichmentRatio)

IOP_seeds_GO_summary <- IOP_seeds_GO %>%
    dplyr::select(geneSet, description, size, overlap, enrichmentRatio, pValue, FDR, userId) %>%
    arrange(-enrichmentRatio)

HPE_seeds_HPO_summary <- HPE_seeds_HPO %>%
    dplyr::select(geneSet, description, size, overlap, enrichmentRatio, pValue, FDR, userId) %>%
    arrange(-enrichmentRatio)

IOP_seeds_HPO_summary <- IOP_seeds_HPO %>%
    dplyr::select(geneSet, description, size, overlap, enrichmentRatio, pValue, FDR, userId) %>%
    arrange(-enrichmentRatio)

HPE_seeds_Reactome_summary <- HPE_seeds_Reactome %>%
    dplyr::select(geneSet, description, size, overlap, enrichmentRatio, pValue, FDR, userId) %>%
    arrange(-enrichmentRatio)

IOP_seeds_Reactome_summary <- IOP_seeds_Reactome %>%
    dplyr::select(geneSet, description, size, overlap, enrichmentRatio, pValue, FDR, userId) %>%
    arrange(-enrichmentRatio)

write_csv(seed_count_per_patient, "article/tables/Figure_2/seed_count_per_patient.csv")
write_csv(HPE_variant_function_summary, "article/tables/Figure_2/HPE_variant_function_summary.csv")
write_csv(IOP_variant_function_summary, "article/tables/Figure_2/IOP_variant_function_summary.csv")
write_csv(HPE_CLNSIG_summary, "article/tables/Figure_2/HPE_CLNSIG_summary.csv")
write_csv(IOP_CLNSIG_summary, "article/tables/Figure_2/IOP_CLNSIG_summary.csv")
write_csv(HPE_PRED_effect_summary, "article/tables/Figure_2/HPE_PRED_effect_summary.csv")
write_csv(IOP_PRED_effect_summary, "article/tables/Figure_2/IOP_PRED_effect_summary.csv")
write_csv(HPE_seeds_GO_summary, "article/tables/Figure_2/HPE_seeds_GO_summary.csv")
write_csv(IOP_seeds_GO_summary, "article/tables/Figure_2/IOP_seeds_GO_summary.csv")
write_csv(HPE_seeds_HPO_summary, "article/tables/Figure_2/HPE_seeds_HPO_summary.csv")
write_csv(IOP_seeds_HPO_summary, "article/tables/Figure_2/IOP_seeds_HPO_summary.csv")
write_csv(HPE_seeds_Reactome_summary, "article/tables/Figure_2/HPE_seeds_Reactome_summary.csv")
write_csv(IOP_seeds_Reactome_summary, "article/tables/Figure_2/IOP_seeds_Reactome_summary.csv")

library(openxlsx)

write_single_xlsx <- function(df, path, sheet_name = "Sheet1") {
    wb <- createWorkbook()
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df, headerStyle = createStyle(textDecoration = "bold"))
    freezePane(wb, sheet_name, firstRow = TRUE)
    setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
    saveWorkbook(wb, path, overwrite = TRUE)
}

write_single_xlsx(seed_count_per_patient, "article/tables/Figure_2/seed_count_per_patient.xlsx", "seed_count_per_patient")
write_single_xlsx(HPE_variant_function_summary, "article/tables/Figure_2/HPE_variant_function_summary.xlsx", "HPE_variant_function")
write_single_xlsx(IOP_variant_function_summary, "article/tables/Figure_2/IOP_variant_function_summary.xlsx", "IOP_variant_function")
write_single_xlsx(HPE_CLNSIG_summary, "article/tables/Figure_2/HPE_CLNSIG_summary.xlsx", "HPE_CLNSIG")
write_single_xlsx(IOP_CLNSIG_summary, "article/tables/Figure_2/IOP_CLNSIG_summary.xlsx", "IOP_CLNSIG")
write_single_xlsx(HPE_PRED_effect_summary, "article/tables/Figure_2/HPE_PRED_effect_summary.xlsx", "HPE_PRED_effect")
write_single_xlsx(IOP_PRED_effect_summary, "article/tables/Figure_2/IOP_PRED_effect_summary.xlsx", "IOP_PRED_effect")
write_single_xlsx(HPE_seeds_GO_summary, "article/tables/Figure_2/HPE_seeds_GO_summary.xlsx", "HPE_seeds_GO")
write_single_xlsx(IOP_seeds_GO_summary, "article/tables/Figure_2/IOP_seeds_GO_summary.xlsx", "IOP_seeds_GO")
write_single_xlsx(HPE_seeds_HPO_summary, "article/tables/Figure_2/HPE_seeds_HPO_summary.xlsx", "HPE_seeds_HPO")
write_single_xlsx(IOP_seeds_HPO_summary, "article/tables/Figure_2/IOP_seeds_HPO_summary.xlsx", "IOP_seeds_HPO")
write_single_xlsx(HPE_seeds_Reactome_summary, "article/tables/Figure_2/HPE_seeds_Reactome_summary.xlsx", "HPE_seeds_Reactome")
write_single_xlsx(IOP_seeds_Reactome_summary, "article/tables/Figure_2/IOP_seeds_Reactome_summary.xlsx", "IOP_seeds_Reactome")

# ============================================================
# Seed frequency per gene (post network filtering), HPE vs IOP
# ------------------------------------------------------------
#  One row per gene carried as a seed (network-filtered) by at
#  least one patient in either cohort, with the number of
#  patients carrying it as a seed in each cohort.
# ============================================================

HPE_seed_freq <- lapply(HPE_variants, function(df) {
    unique(df$Gene.refGeneWithVer)
}) %>%
    unlist() %>%
    table() %>%
    as.data.frame() %>%
    setNames(c("gene", "n_seed_HPE"))

IOP_seed_freq <- lapply(IOP_variants, function(df) {
    unique(df$Gene.refGeneWithVer)
}) %>%
    unlist() %>%
    table() %>%
    as.data.frame() %>%
    setNames(c("gene", "n_seed_IOP"))

seed_freq <- full_join(HPE_seed_freq, IOP_seed_freq, by = "gene") %>%
    mutate(
        n_seed_HPE = replace_na(n_seed_HPE, 0L),
        n_seed_IOP = replace_na(n_seed_IOP, 0L)
    ) %>%
    arrange(desc(n_seed_HPE + n_seed_IOP))

write_csv(seed_freq, "article/tables/Figure_2/seed_freq.csv")
write_single_xlsx(seed_freq, "article/tables/Figure_2/seed_freq.xlsx", "seed_freq")


# ============================================================
# Statistical tests: seed retention after network filtering, HPE vs IOP
# ------------------------------------------------------------
# Answers review comment 1.2 (31.8% vs 30.0% of seeds kept, no test).
#  (a) Per-patient Wilcoxon rank-sum tests on the number of seeds
#      before filtering, the number after filtering and the proportion
#      of seeds kept. Main test: patients are the independent units.
#  (b) Same tests restricted to the patients kept for the RWR
#      (>= 1 seed in the network), i.e. the 178 HPE / 260 IOP.
#  (c) Pooled proportion of seeds kept (sum after / sum before), with
#      a chi-square test. Reference only: the seeds of one patient are
#      not independent.
#  (d) Fisher's exact test on the number of patients left without any
#      seed in the network.
# ============================================================

seed_retention <- seed_count_summary %>%
    dplyr::select(cohort, patient, filter_status, n_seeds) %>%
    pivot_wider(names_from = filter_status, values_from = n_seeds) %>%
    dplyr::rename(n_before = Before, n_after = After) %>%
    mutate(
        prop_kept = ifelse(n_before > 0, n_after / n_before, NA_real_),
        no_seed_after = n_after == 0
    )

wilcox_hpe_iop <- function(df, var, subset_label) {
    x <- na.omit(df[[var]][df$cohort == "HPE"])
    y <- na.omit(df[[var]][df$cohort == "IOP"])
    wt <- suppressWarnings(wilcox.test(x, y, exact = FALSE, conf.int = TRUE))
    tibble(
        subset = subset_label,
        variable = var,
        test = "Wilcoxon rank-sum",
        n_HPE = length(x),
        n_IOP = length(y),
        value_type = "median (mean)",
        HPE_value = sprintf("%.3g (%.3g)", median(x), mean(x)),
        IOP_value = sprintf("%.3g (%.3g)", median(y), mean(y)),
        effect = unname(wt$estimate), # Hodges-Lehmann location shift, HPE - IOP
        effect_type = "location shift (HPE - IOP)",
        p_value = wt$p.value
    )
}

# (a) all patients
tests_all <- map_dfr(
    c("n_before", "n_after", "prop_kept"),
    ~ wilcox_hpe_iop(seed_retention, .x, "all patients")
)

# (b) patients kept for the RWR (>= 1 seed in the network)
seed_retention_kept <- filter(seed_retention, n_after > 0)
tests_kept <- map_dfr(
    c("n_before", "n_after", "prop_kept"),
    ~ wilcox_hpe_iop(seed_retention_kept, .x, "patients with >= 1 network seed")
)

# (c) pooled proportion of seeds kept
pooled <- seed_retention %>%
    group_by(cohort) %>%
    summarise(kept = sum(n_after), lost = sum(n_before - n_after), .groups = "drop")
pooled_tab <- as.matrix(pooled[, c("kept", "lost")])
rownames(pooled_tab) <- pooled$cohort
chi_pooled <- chisq.test(pooled_tab)
test_pooled <- tibble(
    subset = "all patients",
    variable = "pooled seeds kept",
    test = "Chi-square (reference only)",
    n_HPE = sum(pooled_tab["HPE", ]),
    n_IOP = sum(pooled_tab["IOP", ]),
    value_type = "kept / total",
    HPE_value = sprintf("%d / %d (%.1f%%)", pooled_tab["HPE", "kept"], sum(pooled_tab["HPE", ]),
                        100 * pooled_tab["HPE", "kept"] / sum(pooled_tab["HPE", ])),
    IOP_value = sprintf("%d / %d (%.1f%%)", pooled_tab["IOP", "kept"], sum(pooled_tab["IOP", ]),
                        100 * pooled_tab["IOP", "kept"] / sum(pooled_tab["IOP", ])),
    effect = pooled_tab["HPE", "kept"] / sum(pooled_tab["HPE", ]) -
        pooled_tab["IOP", "kept"] / sum(pooled_tab["IOP", ]),
    effect_type = "difference in proportion (HPE - IOP)",
    p_value = chi_pooled$p.value
)

# (d) patients without any seed in the network
no_seed_tab <- table(
    cohort = seed_retention$cohort,
    no_seed = factor(seed_retention$no_seed_after, levels = c(TRUE, FALSE))
)
fisher_no_seed <- fisher.test(no_seed_tab)
test_no_seed <- tibble(
    subset = "all patients",
    variable = "patients without network seed",
    test = "Fisher's exact",
    n_HPE = sum(no_seed_tab["HPE", ]),
    n_IOP = sum(no_seed_tab["IOP", ]),
    value_type = "without seed / total",
    HPE_value = sprintf("%d / %d (%.1f%%)", no_seed_tab["HPE", "TRUE"], sum(no_seed_tab["HPE", ]),
                        100 * no_seed_tab["HPE", "TRUE"] / sum(no_seed_tab["HPE", ])),
    IOP_value = sprintf("%d / %d (%.1f%%)", no_seed_tab["IOP", "TRUE"], sum(no_seed_tab["IOP", ]),
                        100 * no_seed_tab["IOP", "TRUE"] / sum(no_seed_tab["IOP", ])),
    effect = unname(fisher_no_seed$estimate),
    effect_type = "odds ratio (HPE vs IOP)",
    p_value = fisher_no_seed$p.value
)

seed_retention_tests <- bind_rows(tests_all, tests_kept, test_pooled, test_no_seed)
print(seed_retention_tests, width = Inf)

write_csv(seed_retention, "article/tables/Figure_2/seed_retention_per_patient.csv")
write_csv(seed_retention_tests, "article/tables/Figure_2/seed_retention_tests.csv")
write_single_xlsx(seed_retention_tests, "article/tables/Figure_2/seed_retention_tests.xlsx", "seed_retention_tests")


# ============================================================
# Established HPE genes: seed occurrences (network-filtered), HPE vs IOP
# ============================================================

HPE_genes <- read.csv("article/transit_data/HPE_genes.csv",header = TRUE)
HPE_genes
HPE_genes_seeds <- HPE_genes %>%
    dplyr::select(-any_of(c("n_seed_HPE", "n_seed_IOP"))) %>%
    left_join(mutate(seed_freq, gene = as.character(gene)), by = "gene") %>%
    mutate(
        n_seed_HPE = replace_na(n_seed_HPE, 0L),
        n_seed_IOP = replace_na(n_seed_IOP, 0L)
    )

HPE_genes_seeds %>%
    as_tibble() %>%
    filter(n_seed_HPE > 0 | n_seed_IOP > 0) %>%
    dplyr::select(gene, n_seed_HPE, n_seed_IOP) %>%
    print(n = Inf)

write_csv(HPE_genes_seeds, "article/tables/Figure_2/HPE_genes_seeds.csv")
write_single_xlsx(HPE_genes_seeds, "article/tables/Figure_2/HPE_genes_seeds.xlsx", "HPE_genes_seeds")
