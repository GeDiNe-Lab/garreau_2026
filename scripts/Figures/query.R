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
library(scales)

# Setting working directory
rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-4) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()

# Loading custom functions
source("article/scripts/custom_fct.R")
source("article/scripts/RECENT_patient_gene_pvalue.R") # process_group() / process_patient()
source("article/scripts/RECENT_rwr_normalization.R") # build_permutation_bank() / normalize_rwr_score()
source("article/scripts/RECENT_webgestalt_ORA.R")

load("article/transit_data/data_3p.RData")

load("article/transit_data/RWR_results_case.RData")
HPE_results <- as.matrix(RWR_results_z)
rm(RWR_results_z)
load("article/transit_data/iop_RWR_results_case.RData")
IOP_results <- as.matrix(RWR_results_z)
rm(RWR_results_z)

load("article/transit_data/iop_patient_variants.RData")
IOP_patient_variants <- patient_variants
rm(patient_variants)
load("article/transit_data/patient_variants.RData")
HPE_patient_variants <- patient_variants
rm(patient_variants)

# Precomputed permutation bank (build_permutation_bank(), cf. RECENT_results_analyses.R)
# -> reused as-is for both cohorts, no RWR recomputation
load("article/transit_data/permutation_bank.RData") # -> `bank`

coex_genes <- rownames(data_3p$coex_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

context_genes <- rownames(data_3p$norm_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

HPE_seeds <- lapply(HPE_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(HPE_seeds) <- colnames(HPE_results)

IOP_seeds <- lapply(IOP_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(IOP_seeds) <- colnames(IOP_results)

source("article/scripts/RECENT_patient_gene_pvalue.R") # process_group() / process_patient()

signal_genes_HPE <- filter(process_group(
    RWR_results = HPE_results,
    seeds = HPE_seeds,
    permutation_bank = bank
), signal_final)$gene

enrich_hpo_HPE <- run_ora_hpo(signal_genes_HPE)
View(enrich_hpo_HPE)

enrich_gopb_HPE <- run_ora_gobp(signal_genes_HPE)
View(enrich_gopb_HPE)

signal_genes_IOP <- filter(process_group(
    RWR_results = IOP_results,
    seeds = IOP_seeds,
    permutation_bank = bank
), signal_final)$gene

enrich_hpo_IOP <- run_ora_hpo(signal_genes_IOP)
View(enrich_hpo_IOP)

enrich_gopb_IOP <- run_ora_gobp(signal_genes_IOP)
View(enrich_gopb_IOP)

# ============================================================
# Etape 1 : pour CHAQUE cohorte, filtrer les genes qui passent le
#           seuil de significativite statistique via process_group()
#           (signal_final : pval_bh < 0.01 & z_score >= 95e pct),
#           puis garder le top 50 en termes de score moyen
#           (mean_score_raw). Les seeds ne sont PAS exclus : ils
#           sont deja neutralises en amont par process_group() qui
#           remplace leur valeur par la mediane des non-seeds avant
#           de calculer mean_score_raw / z_score.
# ============================================================
get_top_significant <- function(RWR_results, seeds, permutation_bank, n_top = 50) {
    group_stats <- process_group(
        RWR_results = RWR_results,
        seeds = seeds,
        permutation_bank = permutation_bank
    )

    top_genes <- group_stats %>%
        filter(signal_final) %>%
        arrange(desc(mean_score_raw)) %>%
        slice_head(n = n_top) %>%
        pull(gene)

    if (length(top_genes) < n_top) {
        warning(
            length(top_genes), " genes seulement passent le seuil de significativite ",
            "(sur ", n_top, " demandes)"
        )
    }

    list(group_stats = group_stats, top_genes = top_genes)
}


HPE_top <- get_top_significant(HPE_results, HPE_seeds, bank, n_top = 50)
IOP_top <- get_top_significant(IOP_results, IOP_seeds, bank, n_top = 50)

cat("HPE : ", length(HPE_top$top_genes), " genes significatifs (top 50 demande)\n", sep = "")
cat("IOP : ", length(IOP_top$top_genes), " genes significatifs (top 50 demande)\n", sep = "")

# ============================================================
# Etape 2 : z-scores individuels calcules pour TOUTE la cohorte en
#           un seul appel vectorise a normalize_rwr_score() (N = 1
#           par patient, applique colonne par colonne via apply()),
#           plutot qu'un lapply patient par patient + process_patient()
#           (qui refait en plus tout le pipeline p-value/locfdr,
#           inutile ici puisqu'on veut seulement la distribution des
#           z-scores). Les patients avec au moins un NA sont exclus
#           avant normalisation. Les seeds ne sont plus exclus des
#           genes consideres : leur score est impute par la mediane
#           des non-seeds, PAR PATIENT, exactement comme dans
#           process_group(), avant le calcul du z-score individuel.
# ============================================================
build_patient_distribution <- function(RWR_results, seeds, permutation_bank, gene_set) {
    # Exclure les patients avec au moins un NA dans leur profil
    valid_patients <- colnames(RWR_results)[colSums(is.na(RWR_results)) == 0]
    RWR_valid <- RWR_results[, valid_patients, drop = FALSE]

    # Imputation des seeds par la mediane des non-seeds, patient par patient
    RWR_imputed <- sapply(valid_patients, function(pid) {
        scores <- RWR_valid[, pid]
        is_seed <- rownames(RWR_valid) %in% seeds[[pid]]
        scores[is_seed] <- median(scores[!is_seed], na.rm = TRUE)
        scores
    })
    rownames(RWR_imputed) <- rownames(RWR_valid)

    # Z-score de chaque patient (N = 1), toute la cohorte en un seul apply()
    z_matrix <- apply(RWR_imputed, 2, function(patient_scores) {
        normalize_rwr_score(
            patient_scores,
            permutation_bank = permutation_bank,
            verbose = FALSE
        )
    })
    rownames(z_matrix) <- rownames(RWR_imputed)

    z_matrix[gene_set, , drop = FALSE] %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        pivot_longer(cols = -gene, names_to = "patient", values_to = "z_score")
}

HPE_dist <- build_patient_distribution(HPE_results, HPE_seeds, bank, HPE_top$top_genes)
IOP_dist <- build_patient_distribution(IOP_results, IOP_seeds, bank, IOP_top$top_genes)

# ============================================================
# Etape 3 : boxplot par gene (dispersion des z-scores across
#           patients), jitter uniquement sur les outliers (regle
#           de Tukey), genes ordonnes par z-score moyen, meme style
#           que RECENT_score_distrib.R.
# ============================================================
plot_gene_score_distribution <- function(patient_stats, seeds, title_suffix, fill_color = "#2E86AB") {
    df <- patient_stats %>%
        group_by(gene) %>%
        mutate(mean_z = mean(z_score, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(gene = fct_reorder(gene, mean_z))

    df <- df %>%
        group_by(gene) %>%
        mutate(
            q1 = quantile(z_score, 0.25, na.rm = TRUE),
            q3 = quantile(z_score, 0.75, na.rm = TRUE),
            iqr = q3 - q1,
            lower_whisker = q1 - 1.5 * iqr,
            upper_whisker = q3 + 1.5 * iqr,
            is_outlier = z_score < lower_whisker | z_score > upper_whisker
        ) %>%
        ungroup()

    outliers_df <- df %>% filter(is_outlier)

    # Mediane des outliers superieurs a la moyenne, par gene (outliers du cote
    # haut de la distribution uniquement, i.e. z_score > mean_z du gene)
    upper_outlier_median_df <- outliers_df %>%
        filter(z_score > mean_z) %>%
        group_by(gene) %>%
        summarise(median_upper_outlier = median(z_score, na.rm = TRUE), .groups = "drop")

    # Nombre de patients pour lesquels chaque gene est un seed (donc impute a
    # la mediane des non-seeds pour ce patient, cf Etape 2) : rattache
    # directement au nom du gene sur l'axe des y, plutot qu'affiche a part
    seed_count_lookup <- sapply(levels(df$gene), function(g) {
        sum(sapply(seeds, function(s) g %in% s))
    })
    gene_labels <- setNames(
        paste0(names(seed_count_lookup), " (", seed_count_lookup, ")"),
        names(seed_count_lookup)
    )

    ggplot(df, aes(y = gene, x = z_score)) +
        geom_boxplot(fill = fill_color, alpha = 0.4, outlier.shape = NA, width = 0.6) +
        geom_jitter(
            data = outliers_df,
            width = 0.15, alpha = 0.6, size = 1.4, color = "#A23B72"
        ) +
        geom_point(
            data = upper_outlier_median_df,
            aes(y = gene, x = median_upper_outlier),
            inherit.aes = FALSE,
            shape = 23, size = 2.8, fill = "gold", color = "black", stroke = 0.6
        ) +
        scale_y_discrete(labels = gene_labels) +
        labs(
            title = paste0("Distribution of RWR z-scores across patients \u2014 top 50 significant genes (", title_suffix, ")"),
            subtitle = "Genes filtered via process_group() (signal_final), seeds imputed to per-patient non-seed median \u00b7 ranked by mean raw score \u00b7 jitter shows outliers only \u00b7 (n) next to gene = seed count \u00b7 \u25c6 gold = median of upper outliers (> gene mean)",
            y = "Gene",
            x = "RWR z-score"
        ) +
        xlim(-275, 1600) +
        custom_theme() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
            plot.subtitle = element_text(hjust = 0.5, size = 9),
            axis.title = element_text(size = 12),
            axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 9),
            panel.grid.minor = element_blank()
        )
}

HPE_plot <- plot_gene_score_distribution(HPE_dist, HPE_seeds, "HPE / case", fill_color = "#2E86AB")
HPE_plot
ggsave("article/figures/Figure_3A_HPE_top50_zscore_distribution.png", plot = HPE_plot, width = 12, height = 9, dpi = 300)

IOP_plot <- plot_gene_score_distribution(IOP_dist, IOP_seeds, "IOP / control", fill_color = "#A23B72")
IOP_plot
ggsave("article/figures/Figure_3A_IOP_top50_zscore_distribution.png", plot = IOP_plot, width = 12, height = 9, dpi = 300)

# ============================================================
# Etape 4 : z-score des signal genes (signal_final, seuil donne par
#           process_group()) vs. degree centrality dans le reseau.
#           On reutilise directement group_stats calcule a l'Etape 1
#           (HPE_top$group_stats / IOP_top$group_stats) : le z-score
#           de cohorte y est deja calcule sur les scores IMPUTES
#           (seed -> mediane des non-seeds par patient), donc cohabile
#           avec le traitement des seeds fait aux Etapes 1-3.
# ============================================================

# --- Degree centrality dans le reseau de co-expression ---
g <- graph_from_adjacency_matrix(
    data_3p$coex_case,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
)

centrality_df <- data.frame(
    gene = V(g)$name %>% gene_converter(from = "ENSEMBL", to = "SYMBOL"),
    degree = igraph::degree(g)
)

# --- Genes signal (signal_final == TRUE, seuil de process_group), avec
#     leur z-score de cohorte (deja impute) directement issu de l'Etape 1 ---
signal_zscore_degree <- bind_rows(
    HPE_top$group_stats %>%
        filter(signal_final) %>%
        dplyr::select(gene, z_score) %>%
        mutate(cohort = "HPE (case)"),
    IOP_top$group_stats %>%
        filter(signal_final) %>%
        dplyr::select(gene, z_score) %>%
        mutate(cohort = "IOP (control)")
) %>%
    left_join(centrality_df, by = "gene")

save_zscore_degree_table <- function(df, output_path) {
    table_df <- df %>%
        dplyr::select(gene, score = z_score, degree_centrality = degree) %>%
        dplyr::arrange(dplyr::desc(score))

    write.csv(table_df, output_path, row.names = FALSE)
    invisible(table_df)
}

plot_zscore_vs_degree <- function(df, title_suffix, point_color = "#2E86AB") {
    cor_test <- cor.test(df$degree, df$z_score, method = "pearson")
    cor_label <- paste0(
        "Pearson r = ", round(cor_test$estimate, 2),
        ", p = ", signif(cor_test$p.value, 3)
    )

    ggplot(df, aes(x = degree, y = z_score)) +
        geom_point(alpha = 0.7, size = 2, color = point_color) +
        geom_smooth(method = "lm", se = TRUE, linewidth = 0.6, color = point_color) +
        ggrepel::geom_text_repel(
            aes(label = gene),
            size = 2.8, max.overlaps = Inf, seed = 42
        ) +
        labs(
            title = paste0("Signal genes: cohort-level RWR z-score vs. degree centrality (", title_suffix, ")"),
            subtitle = paste0(
                "Signal genes = signal_final (process_group) \u00b7 z-score = group-level z-score from process_group() (seeds imputed to per-patient non-seed median)\n",
                cor_label
            ),
            x = "Degree centrality",
            y = "RWR z-score (cohort-level)"
        ) +
        custom_theme() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
            plot.subtitle = element_text(hjust = 0.5, size = 9),
            axis.title = element_text(size = 12),
            panel.grid.minor = element_blank()
        )
}

HPE_zscore_degree_plot <- plot_zscore_vs_degree(
    signal_zscore_degree %>% filter(cohort == "HPE (case)"),
    "HPE / case",
    point_color = "#2E86AB"
)
HPE_zscore_degree_plot
ggsave("article/figures/Figure_3X_HPE_signal_genes_zscore_vs_degree.png", plot = HPE_zscore_degree_plot, width = 8, height = 6, dpi = 300)
save_zscore_degree_table(
    signal_zscore_degree %>% filter(cohort == "HPE (case)"),
    "article/figures/Figure_3X_HPE_signal_genes_zscore_vs_degree_table.csv"
)

IOP_zscore_degree_plot <- plot_zscore_vs_degree(
    signal_zscore_degree %>% filter(cohort == "IOP (control)"),
    "IOP / control",
    point_color = "#A23B72"
)
IOP_zscore_degree_plot
ggsave("article/figures/Figure_3X_IOP_signal_genes_zscore_vs_degree.png", plot = IOP_zscore_degree_plot, width = 8, height = 6, dpi = 300)
save_zscore_degree_table(
    signal_zscore_degree %>% filter(cohort == "IOP (control)"),
    "article/figures/Figure_3X_IOP_signal_genes_zscore_vs_degree_table.csv"
)


enrich_reactome <- run_ora_reactome(filter(process_group(
    RWR_results = HPE_results,
    seeds = HPE_seeds,
    permutation_bank = bank
), signal_final)$gene)

enrich_gobp <- run_ora_gobp(filter(process_group(
    RWR_results = HPE_results,
    seeds = HPE_seeds,
    permutation_bank = bank
), signal_final)$gene)

source("article/scripts/RECENT_patient_gene_pvalue.R") # process_group() / process_patient()
source("article/scripts/RECENT_webgestalt_ORA.R")

signal_genes_HPE <- filter(process_group(
    RWR_results = HPE_results,
    seeds = HPE_seeds,
    permutation_bank = bank
), signal_final)$gene
signal_genes_HPE %>% sort()

signal_genes_IOP <- filter(process_group(
    RWR_results = IOP_results,
    seeds = IOP_seeds,
    permutation_bank = bank
), signal_final)$gene
signal_genes_IOP %>% sort()

enrich_hpo_HPE <- run_ora_hpo(signal_genes_HPE)
enrich_gopb_HPE <- run_ora_gopb(signal_genes_HPE)

enrich_hpo_IOP <- run_ora_hpo(signal_genes_IOP)
enrich_gopb_IOP <- run_ora_gopb(signal_genes_IOP)
