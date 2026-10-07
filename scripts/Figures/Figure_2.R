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
        title = "Number of seeds per patient, before/after network filtering",
        x = NULL,
        y = "Number of seeds"
    ) +
    custom_theme() +
    theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5)
    )

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
# make_pie_chart()
# ------------------------------------------------------------
# ggplot2 pie chart for a frequency table (2 columns: category +
# Freq, regardless of the 1st column's name -- automatically
# renamed to "category").
#
# @param freq_df   2-column data.frame (category, Freq)
# @param title     plot title
# @param top_n     if provided, keeps only the top_n categories by
#                  frequency and groups the rest under "Other"
#                  (useful when there are many categories, e.g. CLNSIG)
# ============================================================
make_pie_chart <- function(freq_df, title = NULL, top_n = NULL, other_label = "Other") {
    df <- freq_df
    colnames(df)[1] <- "category"
    df$category <- as.character(df$category)
    df$category <- gsub("_", " ", df$category)

    if (!is.null(top_n) && nrow(df) > top_n) {
        df <- df %>% arrange(desc(Freq))
        df <- bind_rows(
            df[seq_len(top_n), ],
            data.frame(category = other_label, Freq = sum(df$Freq[(top_n + 1):nrow(df)]))
        )
    }

    df <- df %>%
        arrange(desc(Freq)) %>%
        mutate(
            fraction = Freq / sum(Freq),
            label = sprintf("%s (%.1f%%)", category, fraction * 100),
            label = factor(label, levels = label)
        )

    ggplot(df, aes(x = "", y = Freq, fill = label)) +
        geom_col(width = 1, color = "white") +
        coord_polar(theta = "y") +
        scale_fill_viridis_d() +
        labs(title = title, fill = NULL) +
        theme_void() +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            legend.text = element_text(size = 9),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
        )
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

p_HPE_func <- make_pie_chart(HPE_variants_function, title = "HPE - Variant genomic location")
p_HPE_clnsig <- make_pie_chart(HPE_CLNSIG, title = "HPE - ClinVar significance", top_n = 8)
p_HPE_pred <- make_pie_chart(HPE_PRED_effect, title = "HPE - Predicted effect", top_n = 8)

p_HPE_go <- plot_ora_dotplot(HPE_seeds_GO, n = 20, title = "HPE - Top 20 GO-BP")
p_HPE_hpo <- plot_ora_dotplot(HPE_seeds_HPO, n = 20, title = "HPE - Top 20 HPO")
p_HPE_reactome <- plot_ora_dotplot(HPE_seeds_Reactome, n = 20, title = "HPE - Top 20 Reactome pathways")

p_IOP_func <- make_pie_chart(IOP_variants_function, title = "IOP - Variant genomic location")
p_IOP_clnsig <- make_pie_chart(IOP_CLNSIG, title = "IOP - ClinVar significance", top_n = 8)
p_IOP_pred <- make_pie_chart(IOP_PRED_effect, title = "IOP - Predicted effect", top_n = 8)

p_IOP_go <- plot_ora_dotplot(IOP_seeds_GO, n = 20, title = "IOP - Top 20 GO-BP")
p_IOP_hpo <- plot_ora_dotplot(IOP_seeds_HPO, n = 20, title = "IOP - Top 20 HPO")
p_IOP_reactome <- plot_ora_dotplot(IOP_seeds_Reactome, n = 20, title = "IOP - Top 20 Reactome pathways")

# ============================================================
# Saving plots
# ============================================================
ggsave("article/figures/Figure_2/seed_count_distribution.png", plot = p_seed_distribution, width = 9, height = 6, dpi = 300)

ggsave("article/figures/Figure_2/HPE_Func_refGene_pie.png", plot = p_HPE_func, width = 12, height = 6, dpi = 300)
ggsave("article/figures/Figure_2/HPE_CLNSIG_pie.png", plot = p_HPE_clnsig, width = 12, height = 6, dpi = 300)
ggsave("article/figures/Figure_2/HPE_PRED_effect_pie.png", plot = p_HPE_pred, width = 12, height = 6, dpi = 300)
ggsave("article/figures/Figure_2/HPE_seeds_GO_dotplot.png", plot = p_HPE_go, width = 10, height = 8, dpi = 300)
ggsave("article/figures/Figure_2/HPE_seeds_HPO_dotplot.png", plot = p_HPE_hpo, width = 10, height = 8, dpi = 300)
ggsave("article/figures/Figure_2/HPE_seeds_Reactome_dotplot.png", plot = p_HPE_reactome, width = 10, height = 8, dpi = 300)

ggsave("article/figures/Figure_2/IOP_Func_refGene_pie.png", plot = p_IOP_func, width = 12, height = 6, dpi = 300)
ggsave("article/figures/Figure_2/IOP_CLNSIG_pie.png", plot = p_IOP_clnsig, width = 12, height = 6, dpi = 300)
ggsave("article/figures/Figure_2/IOP_PRED_effect_pie.png", plot = p_IOP_pred, width = 12, height = 6, dpi = 300)
ggsave("article/figures/Figure_2/IOP_seeds_GO_dotplot.png", plot = p_IOP_go, width = 10, height = 8, dpi = 300)
ggsave("article/figures/Figure_2/IOP_seeds_HPO_dotplot.png", plot = p_IOP_hpo, width = 10, height = 8, dpi = 300)
ggsave("article/figures/Figure_2/IOP_seeds_Reactome_dotplot.png", plot = p_IOP_reactome, width = 10, height = 8, dpi = 300)
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
