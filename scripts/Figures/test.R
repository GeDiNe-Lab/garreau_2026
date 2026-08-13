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

load("article/transit_data/HPE_GENES_RWR_results.RData")
HPE_results <- as.matrix(HPE_genes_results)
rm(HPE_genes_results)


# Precomputed permutation bank (build_permutation_bank(), cf. RECENT_results_analyses.R)
# -> reused as-is for both cohorts, no RWR recomputation
load("article/transit_data/permutation_bank.RData") # -> `bank`

coex_genes <- rownames(data_3p$coex_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

context_genes <- rownames(data_3p$norm_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")


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

cat("HPE : ", length(HPE_top$top_genes), " genes significatifs (top 50 demande)\n", sep = "")


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

# ============================================================
# Etape 5 : enrichissement fonctionnel par sous-echantillonnage
#           iteratif de patients, cohortes HPE et IOP, meme logique
#           que RECENT_results_analyses.R (process_group() sur un
#           echantillon aleatoire de patients -> genes signal -> ORA).
#           HPO ET GO-BP sont calcules DANS LA MEME BOUCLE, sur le
#           MEME echantillon de patients / memes genes signal a
#           chaque iteration (pas de re-echantillonnage separe par
#           ontologie) -> tables et graphs restent separes par
#           ontologie, mais directement comparables entre eux.
# ------------------------------------------------------------
# Pour n_iter iterations : tirer un echantillon aleatoire de
# sample_size patients, calculer les genes signal via process_group()
# (seuil de significativite deja applique en interne : pval_bh < 0.01
# & z_score >= 95e pct -> signal_final), lancer l'ORA (une fois par
# ontologie fournie dans ora_funs) sur ces genes. Un terme n'apparait
# dans les resultats d'une iteration QUE s'il a deja passe le seuil
# de FDR de WebGestaltR a la detection ; aucun filtre de
# significativite supplementaire n'est donc applique ici sur les
# resultats d'ORA eux-memes.
# Les patients avec au moins un NA dans leur profil RWR sont exclus
# AVANT tout calcul (donc avant meme l'echantillonnage).
# ============================================================

# --- Parametres d'iteration : a ajuster librement pour tester
#     difIOPentes configurations (n_iter, sample_size) ---
n_iter <- 100 # nombre d'iterations de sous-echantillonnage
sample_size <- 15 # nombre de patients tires aleatoirement par iteration
ora_min_signal_genes <- 5 # nb minimum de genes signal pour lancer l'ORA sur une iteration
sampling_seed <- 123 # graine pour la reproductibilite du tirage aleatoire

run_iterative_ora <- function(RWR_results, seeds, permutation_bank, context_genes,
                              n_iter, sample_size,
                              ora_funs = list(HPO = run_ora_hpo, GOBP = run_ora_gobp),
                              min_signal_genes = 5, sampling_seed = 123,
                              cohort_label = "") {
    # Exclure les patients avec au moins un NA dans leur profil RWR, AVANT tout
    valid_patients <- colnames(RWR_results)[colSums(is.na(RWR_results)) == 0]
    if (length(valid_patients) < sample_size) {
        stop(
            cohort_label, " : seulement ", length(valid_patients),
            " patients valides (sans NA) pour un echantillon de taille ", sample_size
        )
    }

    ontology_names <- names(ora_funs)
    enrichment_lists <- setNames(
        lapply(ontology_names, function(x) vector("list", n_iter)),
        ontology_names
    )

    set.seed(sampling_seed)

    for (i in seq_len(n_iter)) {
        cat("===", cohort_label, "- iteration", i, "/", n_iter, "===\n")

        sampled_patients <- sample(valid_patients, sample_size)

        # --- Genes signal, UNE SEULE FOIS par iteration, partages entre ontologies ---
        group_result <- tryCatch(
            process_group(
                patient_ids = sampled_patients,
                RWR_results = RWR_results,
                seeds = seeds,
                permutation_bank = permutation_bank
            ),
            error = function(e) {
                message("Iteration ", i, " : erreur process_group - ", conditionMessage(e))
                NULL
            }
        )

        if (is.null(group_result)) next

        signal_genes <- group_result$gene[group_result$signal_final]

        if (length(signal_genes) < min_signal_genes) {
            message(
                "Iteration ", i, " : moins de ", min_signal_genes,
                " genes signal (", length(signal_genes), "), ORA ignoree"
            )
            next
        }

        # --- ORA pour chaque ontologie, sur les MEMES genes signal ---
        for (ont in ontology_names) {
            ora_res <- tryCatch(
                ora_funs[[ont]](
                    genelist = signal_genes,
                    referenceGene = context_genes,
                    projectName = paste0("ORA_", cohort_label, "_", ont, "_iter", i),
                    is.output = FALSE
                ),
                error = function(e) {
                    message("Iteration ", i, " (", ont, ") : erreur ORA - ", conditionMessage(e))
                    NULL
                }
            )

            if (!is.null(ora_res) && nrow(ora_res) > 0) {
                ora_res$iteration <- i
                ora_res$n_genes <- length(signal_genes)
                enrichment_lists[[ont]][[i]] <- ora_res
            }
        }
    }

    # --- Resume par ontologie (n_iterations, pct_iterations, mean/sd enrichmentRatio, mean FDR) ---
    summarise_ontology <- function(enrichment_list, ont_label) {
        enrichment_all <- bind_rows(enrichment_list)
        if (nrow(enrichment_all) == 0) {
            warning(cohort_label, " (", ont_label, ") : aucun terme enrichi trouve sur les ", n_iter, " iterations")
            return(list(all_iterations = enrichment_all, summary = NULL))
        }

        summary_df <- enrichment_all %>%
            group_by(geneSet, description) %>%
            summarise(
                n_iterations = n_distinct(iteration),
                pct_iterations = 100 * n_iterations / n_iter,
                mean_enrichmentRatio = mean(enrichmentRatio, na.rm = TRUE),
                sd_enrichmentRatio = sd(enrichmentRatio, na.rm = TRUE),
                mean_FDR = mean(FDR, na.rm = TRUE),
                .groups = "drop"
            ) %>%
            arrange(desc(mean_enrichmentRatio))

        list(all_iterations = enrichment_all, summary = summary_df)
    }

    setNames(
        lapply(ontology_names, function(ont) summarise_ontology(enrichment_lists[[ont]], ont)),
        ontology_names
    )
}

HPE_ora <- run_iterative_ora(
    HPE_results, HPE_seeds, bank, context_genes,
    n_iter = n_iter, sample_size = sample_size,
    ora_funs = list(HPO = run_ora_hpo, GOBP = run_ora_gobp),
    min_signal_genes = ora_min_signal_genes,
    sampling_seed = sampling_seed, cohort_label = "HPE"
)

IOP_ora <- run_iterative_ora(
    IOP_results, IOP_seeds, bank, context_genes,
    n_iter = n_iter, sample_size = sample_size,
    ora_funs = list(HPO = run_ora_hpo, GOBP = run_ora_gobp),
    min_signal_genes = ora_min_signal_genes,
    sampling_seed = sampling_seed, cohort_label = "IOP"
)

# --- Tables separees par cohorte x ontologie ---
write.csv(HPE_ora$HPO$all_iterations, "article/enrichment/HPE_hpo_iterative_all.csv", row.names = FALSE)
write.csv(HPE_ora$HPO$summary, "article/enrichment/HPE_hpo_iterative_summary.csv", row.names = FALSE)
write.csv(HPE_ora$GOBP$all_iterations, "article/enrichment/HPE_gobp_iterative_all.csv", row.names = FALSE)
write.csv(HPE_ora$GOBP$summary, "article/enrichment/HPE_gobp_iterative_summary.csv", row.names = FALSE)

write.csv(IOP_ora$HPO$all_iterations, "article/enrichment/IOP_hpo_iterative_all.csv", row.names = FALSE)
write.csv(IOP_ora$HPO$summary, "article/enrichment/IOP_hpo_iterative_summary.csv", row.names = FALSE)
write.csv(IOP_ora$GOBP$all_iterations, "article/enrichment/IOP_gobp_iterative_all.csv", row.names = FALSE)
write.csv(IOP_ora$GOBP$summary, "article/enrichment/IOP_gobp_iterative_summary.csv", row.names = FALSE)


# ============================================================
# Dot plot : top 20 termes classes par mean_enrichmentRatio.
#   - position x        : enrichment ratio moyen sur les iterations
#   - barre d'erreur     : +/- sd de l'enrichment ratio (seule barre
#                          d'erreur du plot, sur la mesure qui compte)
#   - taille du point    : % d'iterations ou le terme est ressorti
#                          (indicatif, sans marge d'erreur)
#   - couleur du point   : FDR moyenne (indicatif, sans marge d'erreur)
# ============================================================
plot_iterative_enrichment <- function(summary_df, n_top = 20, title_suffix = "") {
    top_terms <- summary_df %>%
        arrange(desc(mean_enrichmentRatio)) %>%
        slice_head(n = n_top) %>%
        mutate(description = fct_reorder(description, mean_enrichmentRatio))

    ggplot(top_terms, aes(y = description, x = mean_enrichmentRatio)) +
        geom_segment(
            aes(
                x = pmax(mean_enrichmentRatio - sd_enrichmentRatio, 0),
                xend = mean_enrichmentRatio + sd_enrichmentRatio,
                yend = description
            ),
            color = "grey50", linewidth = 0.6
        ) +
        geom_point(aes(size = pct_iterations, color = mean_FDR)) +
        scale_color_gradient(low = "#D7263D", high = "#2E86AB", name = "Mean FDR") +
        scale_size_continuous(name = "% iterations", range = c(2, 9)) +
        labs(
            title = paste0("Iterative functional enrichment \u2014 top 20 terms (", title_suffix, ")"),
            subtitle = paste0(
                n_iter, " iterations \u00b7 ", sample_size, " patients/iteration \u00b7 ",
                "error bars = SD of enrichment ratio across iterations where the term passed FDR"
            ),
            y = "Term",
            x = "Mean enrichment ratio"
        ) +
        custom_theme() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
            plot.subtitle = element_text(hjust = 0.5, size = 8.5),
            axis.title = element_text(size = 12),
            axis.text.y = element_text(size = 9),
            panel.grid.minor = element_blank()
        )
}

# --- HPO : HPE et IOP ---
HPE_hpo_plot <- plot_iterative_enrichment(HPE_ora$HPO$summary, title_suffix = "HPE / case \u2014 HPO")
HPE_hpo_plot
ggsave("article/figures/Figure_3C_HPE_iterative_HPO_enrichment.png", plot = HPE_hpo_plot, width = 11, height = 8, dpi = 300)

IOP_hpo_plot <- plot_iterative_enrichment(IOP_ora$HPO$summary, title_suffix = "IOP / control \u2014 HPO")
IOP_hpo_plot
ggsave("article/figures/Figure_3C_IOP_iterative_HPO_enrichment.png", plot = IOP_hpo_plot, width = 11, height = 8, dpi = 300)

# --- GO-BP : HPE et IOP ---
HPE_gobp_plot <- plot_iterative_enrichment(HPE_ora$GOBP$summary, title_suffix = "HPE / case \u2014 GO-BP")
HPE_gobp_plot
ggsave("article/figures/Figure_3D_HPE_iterative_GOBP_enrichment.png", plot = HPE_gobp_plot, width = 11, height = 8, dpi = 300)

IOP_gobp_plot <- plot_iterative_enrichment(IOP_ora$GOBP$summary, title_suffix = "IOP / control \u2014 GO-BP")
IOP_gobp_plot
ggsave("article/figures/Figure_3D_IOP_iterative_GOBP_enrichment.png", plot = IOP_gobp_plot, width = 11, height = 8, dpi = 300)
