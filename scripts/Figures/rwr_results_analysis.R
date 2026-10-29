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
library(openxlsx)

setwd(paste(here::here(), "projects", sep = "/"))

source("article/scripts/custom_fct.R")
source("article/scripts/RECENT_patient_gene_pvalue.R") # process_group() / process_patient()
source("article/scripts/RECENT_rwr_normalization.R") # build_permutation_bank() / normalize_rwr_score()
source("article/scripts/RECENT_webgestalt_ORA.R") # run_ora_reactome() / run_ora_gobp() / run_ora_hpo()

dir.create("article/figures/Figure_3", recursive = TRUE, showWarnings = FALSE)
dir.create("article/tables/Figure_3", recursive = TRUE, showWarnings = FALSE)

seed_color_palette <- c("#1dea65", "#2bbe8b", "#24a5a3", "#376d9c", "#1a2a47")

# ============================================================
# write_single_xlsx()
# ------------------------------------------------------------
# Exports a single data frame to a single-sheet .xlsx file, with
# a bold header row, a frozen header row, and auto-fit column
# widths -- used throughout this script as the .xlsx counterpart
# to write_csv()/write.csv().
# ============================================================
write_single_xlsx <- function(df, path, sheet_name = "Sheet1") {
    wb <- createWorkbook()
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, df, headerStyle = createStyle(textDecoration = "bold"))
    freezePane(wb, sheet_name, firstRow = TRUE)
    setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
    saveWorkbook(wb, path, overwrite = TRUE)
}

# ============================================================
# 0. Data loading (once)
# ============================================================
load("article/transit_data/data_3p.RData")

load("article/transit_data/RWR/HPE_RWR_results_.RData")
load("article/transit_data/RWR/IOP_RWR_results_.RData")

load("article/transit_data/patient_variants.RData")
HPE_patient_variants <- patient_variants
rm(patient_variants)
load("article/transit_data/iop_patient_variants.RData")
IOP_patient_variants <- patient_variants
rm(patient_variants)

load("article/transit_data/permutation_bank_cds_weighted.RData")

coex_genes <- rownames(data_3p$coex_case) %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

HPE_seeds <- lapply(HPE_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(HPE_seeds) <- colnames(HPE_results)

IOP_seeds <- lapply(IOP_patient_variants, function(df) {
    filter(df, Gene.refGeneWithVer %in% coex_genes)$Gene.refGeneWithVer
})
names(IOP_seeds) <- colnames(IOP_results)

HPE_seeds %>% unlist() %>% unique() %>% length() # 1,042 unique genes across all HPE patients
IOP_seeds %>% unlist() %>% unique() %>% length() # 1,042 unique genes across all IOP patients
HPE_patient_variants %>% Reduce(rbind, .) %>% nrow() # 1,042 unique genes across all HPE patients


# ============================================================
# 1. Co-expression network: degree + mean correlation per gene
#    (reused in section 4)
# ============================================================

g <- graph_from_adjacency_matrix(
    data_3p$coex_case,
    mode = "undirected", weighted = TRUE, diag = FALSE
)
V(g)$name <- V(g)$name %>% gene_converter(from = "ENSEMBL", to = "SYMBOL")

centrality_df <- data.frame(
    gene     = V(g)$name,
    degree   = igraph::degree(g),
    mean_cor = igraph::strength(g, weights = E(g)$weight) / pmax(igraph::degree(g), 1)
)


# ============================================================
# 2. Cohorts: common setup (group stats, all signal genes)
# ============================================================

get_signal_genes <- function(RWR_results, seeds, permutation_bank) {
    group_stats <- process_profiles(RWR_results = RWR_results, seeds = seeds, permutation_bank = permutation_bank)
    signal_genes <- group_stats %>%
        filter(signal_final) %>%
        arrange(desc(mean_score_raw)) %>%
        pull(gene)
    list(group_stats = group_stats, signal_genes = signal_genes)
}

cohorts <- list(
    HPE = list(results = HPE_results, seeds = HPE_seeds, color = "#2E86AB", label = "HPE / case"),
    IOP = list(results = IOP_results, seeds = IOP_seeds, color = "#A23B72", label = "IOP / control")
)


cohorts <- imap(cohorts, function(co, cname) {
    res <- get_signal_genes(co$results, co$seeds, bank_cds)
    co$group_stats <- res$group_stats
    co$signal_genes <- res$signal_genes
    cat(cname, ": ", length(co$signal_genes), " signal genes\n", sep = "")
    co
})




# ============================================================
# 3. Distribution of individual (per patient) z-scores, used by
#    section 4
# ============================================================
build_patient_distribution <- function(RWR_results, seeds, permutation_bank, gene_set) {
    valid_patients <- colnames(RWR_results)[colSums(is.na(RWR_results)) == 0]
    RWR_valid <- RWR_results[, valid_patients, drop = FALSE]

    RWR_imputed <- sapply(valid_patients, function(pid) {
        scores <- RWR_valid[, pid]
        is_seed <- rownames(RWR_valid) %in% seeds[[pid]]
        scores[is_seed] <- median(scores[!is_seed], na.rm = TRUE)
        scores
    })
    rownames(RWR_imputed) <- rownames(RWR_valid)

    z_matrix <- apply(RWR_imputed, 2, function(patient_scores) {
        normalize_rwr_score(patient_scores, permutation_bank = permutation_bank, verbose = FALSE)
    })
    rownames(z_matrix) <- rownames(RWR_imputed)

    z_matrix[gene_set, , drop = FALSE] %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        pivot_longer(cols = -gene, names_to = "patient", values_to = "z_score")
}


# ============================================================
# 4. Significant genes figure: y = gene, x = z_score, shared
#    x-scale across cohorts (min/max computed on the combined data).
#    - Genes WITHOUT any seed (background genes): minimal
#      representation, no label, no jitter, just a thick segment
#      for the IQR ("mustache box" as a line) and a thin segment
#      extending to the true min/max. Grey, no color scale.
#    - Genes WITH at least one seed: full boxplot, jitter shown
#      only for outliers, diamond = median of upper outliers.
#      Boxplot fill, jitter, diamond, and gene label are all
#      colored on a shared continuous, colorblind-friendly scale
#      (viridis "mako") by seed count.
#    - Genes ordered top to bottom by median of upper outliers.
#    - DISPLAY IS TRUNCATED for readability: only the top n_edge,
#      bottom n_edge, and every seed-carrying gene are shown.
#      Exactly two separators are drawn (never one per internal
#      gap): a thick dark line + count right after the top block,
#      and a matching thick dark line (no text) right before the
#      bottom block. Both separators create a slightly larger
#      visual gap than the normal row spacing. If nothing is
#      omitted (typically the smaller IOP cohort), no separator
#      is drawn at all.
#    A full, untruncated recap table (all genes, every stat, plus
#    a shown_in_plot flag) is exported separately so the biologist
#    can access what the truncated HPE plot leaves out.
# ============================================================

compute_gene_stats <- function(z) {
    z <- z[is.finite(z)]
    q1 <- quantile(z, 0.25, na.rm = TRUE)
    q3 <- quantile(z, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    tibble(
        q1 = q1,
        mean_z = mean(z, na.rm = TRUE),
        median_z = median(z, na.rm = TRUE),
        q3 = q3,
        whisker_low = q1 - 1.5 * iqr,
        whisker_high = q3 + 1.5 * iqr,
        range_min = min(z),
        range_max = max(z)
    )
}

build_full_gene_stats <- function(patient_stats, seeds) {
    seed_count_lookup <- sapply(unique(patient_stats$gene), function(gn) sum(sapply(seeds, function(s) gn %in% s)))
    seed_count_df <- tibble(gene = names(seed_count_lookup), n_seeds = as.numeric(seed_count_lookup))

    gene_stats <- patient_stats %>%
        group_by(gene) %>%
        summarise(compute_gene_stats(z_score), .groups = "drop") %>%
        left_join(seed_count_df, by = "gene")

    df <- patient_stats %>%
        left_join(gene_stats, by = "gene") %>%
        mutate(is_outlier = z_score < whisker_low | z_score > whisker_high)

    outlier_counts <- df %>%
        group_by(gene) %>%
        summarise(
            n_outliers = sum(is_outlier, na.rm = TRUE),
            n_upper_outliers = sum(is_outlier & z_score > mean_z, na.rm = TRUE),
            .groups = "drop"
        )

    upper_outlier_median_df <- df %>%
        filter(is_outlier, z_score > mean_z) %>%
        group_by(gene) %>%
        summarise(median_upper_outlier = median(z_score, na.rm = TRUE), .groups = "drop")

    gene_stats <- gene_stats %>%
        left_join(outlier_counts, by = "gene") %>%
        left_join(upper_outlier_median_df, by = "gene") %>%
        mutate(has_seed = n_seeds > 0) %>%
        arrange(desc(median_upper_outlier)) %>%
        mutate(rank = row_number())

    list(df = df, gene_stats = gene_stats)
}

plot_gene_score_distribution <- function(built, title_suffix, xlim, n_edge = 25) {
    gene_stats <- built$gene_stats
    df <- built$df
    n_total <- nrow(gene_stats)

    top_block <- gene_stats %>% filter(rank <= n_edge)
    bottom_block <- gene_stats %>% filter(rank > n_total - n_edge)
    middle_seed_block <- gene_stats %>% filter(rank > n_edge, rank <= n_total - n_edge, n_seeds > 0)

    kept_genes <- bind_rows(top_block, middle_seed_block, bottom_block) %>%
        distinct(gene, .keep_all = TRUE) %>%
        arrange(rank) %>%
        mutate(
            row_number = row_number(),
            block = case_when(
                rank <= n_edge ~ "top",
                rank > n_total - n_edge ~ "bottom",
                TRUE ~ "middle"
            )
        )

    n_kept <- nrow(kept_genes)
    n_omitted_total <- n_total - n_kept
    show_separators <- n_omitted_total > 0

    if (show_separators) {
        kept_genes <- kept_genes %>%
            mutate(slot = row_number + case_when(
                block == "top" ~ 0,
                block == "middle" ~ 1,
                block == "bottom" ~ 2
            ))
    } else {
        kept_genes <- kept_genes %>% mutate(slot = row_number)
    }

    n_slots <- max(kept_genes$slot)
    kept_genes <- kept_genes %>% mutate(y_pos = n_slots - slot + 1)

    top_sep_y <- NA_real_
    bottom_sep_y <- NA_real_
    if (show_separators) {
        top_last_slot <- kept_genes %>%
            filter(block == "top") %>%
            pull(slot) %>%
            max()
        next_slot_after_top <- kept_genes %>%
            filter(slot > top_last_slot) %>%
            pull(slot) %>%
            min()
        top_sep_slot <- (top_last_slot + next_slot_after_top) / 2
        top_sep_y <- n_slots - top_sep_slot + 1

        bottom_first_slot <- kept_genes %>%
            filter(block == "bottom") %>%
            pull(slot) %>%
            min()
        prev_slot_before_bottom <- kept_genes %>%
            filter(slot < bottom_first_slot) %>%
            pull(slot) %>%
            max()
        bottom_sep_slot <- (bottom_first_slot + prev_slot_before_bottom) / 2
        bottom_sep_y <- n_slots - bottom_sep_slot + 1
    }

    y_pos_lookup <- kept_genes %>% dplyr::select(gene, y_pos)

    gene_stats_disp <- gene_stats %>% inner_join(y_pos_lookup, by = "gene")
    df_disp <- df %>% inner_join(y_pos_lookup, by = "gene")

    outliers_df <- df_disp %>% filter(is_outlier)
    upper_outlier_median_df <- gene_stats_disp %>%
        filter(!is.na(median_upper_outlier)) %>%
        dplyr::select(gene, n_seeds, median_upper_outlier, y_pos)

    stats_with_seed <- gene_stats_disp %>% filter(n_seeds > 0)
    stats_without_seed <- gene_stats_disp %>% filter(n_seeds == 0)
    df_with_seed <- df_disp %>% filter(n_seeds > 0)
    outliers_with_seed <- outliers_df %>% filter(n_seeds > 0)
    upper_outlier_with_seed <- upper_outlier_median_df %>% filter(n_seeds > 0)

    # gene name in bold italic, seed count in bold upright (plotmath expression)
    label_df_seed <- stats_with_seed %>%
        mutate(label = paste0("bolditalic('", gene, "')~bold('(", n_seeds, ")')"))

    p <- ggplot() +
        # Background genes (no seed): thin line to true min/max, thick line for IQR
        geom_segment(
            data = stats_without_seed,
            aes(y = y_pos, yend = y_pos, x = range_min, xend = range_max),
            linewidth = 0.3, color = "grey60"
        ) +
        geom_segment(
            data = stats_without_seed,
            aes(y = y_pos, yend = y_pos, x = q1, xend = q3),
            linewidth = 1.3, color = "grey60", lineend = "butt"
        ) +
        # Seed-carrying genes: full boxplot, colored by seed count
        geom_boxplot(
            data = df_with_seed, aes(y = y_pos, x = z_score, group = y_pos, fill = n_seeds),
            color = "grey30", alpha = 0.85, outlier.shape = NA, width = 0.6
        ) +
        geom_jitter(
            data = outliers_with_seed, aes(y = y_pos, x = z_score, fill = n_seeds),
            shape = 21, color = "black", width = 0.15, height = 0, alpha = 0.85, size = 1.6, stroke = 0.3
        ) +
        geom_point(
            data = upper_outlier_with_seed, aes(y = y_pos, x = median_upper_outlier, fill = n_seeds),
            shape = 23, size = 2.8, color = "black", stroke = 0.6
        ) +
        geom_text(
            data = label_df_seed,
            aes(y = y_pos, x = xlim[1], label = label, color = n_seeds),
            hjust = 1, size = 3, parse = TRUE
        )

    if (show_separators) {
        p <- p +
            geom_hline(yintercept = top_sep_y, linetype = "solid", color = "grey15", linewidth = 1) +
            geom_hline(yintercept = bottom_sep_y, linetype = "solid", color = "grey15", linewidth = 1) +
            geom_text(
                aes(x = mean(xlim), y = top_sep_y),
                label = paste0(n_omitted_total, " genes without seed omitted (middle of distribution)"),
                color = "grey15", fontface = "bold", size = 3.5, vjust = -0.6
            )
    }

    p +
        scale_fill_gradientn(name = "Seed count", colours = seed_color_palette) +
        scale_color_gradientn(colours = seed_color_palette, guide = "none") +
        scale_y_continuous(limits = c(0.4, n_slots + 0.6), expand = c(0, 0)) +
        labs(
            y = NULL, x = "RWR z-score"
        ) +
        coord_cartesian(xlim = xlim, clip = "off") +
        custom_theme() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
            plot.subtitle = element_text(hjust = 0.5, size = 9),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.background = element_rect(fill = "white", color = "black", linewidth = 1),
            axis.title.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            plot.margin = margin(5.5, 5.5, 5.5, 70),
            legend.title = element_text(size = 16, face = "bold"), # titre "Seed count"
            legend.text = element_text(size = 14),                 # valeurs de la barre
            legend.key.height = unit(1, "cm"),                     # hauteur de la barre
            legend.key.width = unit(0.8, "cm")                     # largeur de la barre
        )
}

build_gene_score_recap <- function(gene_stats, n_edge = 25) {
    n_total <- nrow(gene_stats)
    gene_stats %>%
        mutate(shown_in_plot = rank <= n_edge | rank > (n_total - n_edge) | n_seeds > 0) %>%
        arrange(rank) %>%
        dplyr::select(
            gene, has_seed, n_seeds, rank, shown_in_plot,
            median_z, mean_z, q1, q3, whisker_low, whisker_high, range_min, range_max,
            n_outliers, n_upper_outliers, median_upper_outlier
        )
}

plot_gene_score_distribution_all_cohorts <- function(cohorts, permutation_bank, n_edge = 25) {
    dist_list <- map(cohorts, ~ build_patient_distribution(.x$results, .x$seeds, permutation_bank, .x$signal_genes))

    global_xlim <- range(
        map_dbl(dist_list, ~ min(.x$z_score, na.rm = TRUE)),
        map_dbl(dist_list, ~ max(.x$z_score, na.rm = TRUE))
    )

    built_list <- imap(dist_list, function(d, cname) {
        build_full_gene_stats(d, cohorts[[cname]]$seeds)
    })

    plots <- imap(built_list, function(built, cname) {
        plot_gene_score_distribution(built, cohorts[[cname]]$label, global_xlim, n_edge = n_edge)
    })

    recap <- imap(built_list, function(built, cname) {
        build_gene_score_recap(built$gene_stats, n_edge = n_edge) %>% mutate(cohort = cname)
    }) %>%
        bind_rows() %>%
        dplyr::select(cohort, everything())

    list(plots = plots, recap = recap)
}


scale_fill_viridis_c_manual <- scale_fill_gradientn(
    name = "Seed count",
    colours = seed_color_palette
)
scale_color_viridis_c_manual <- scale_color_gradientn(
    colours = seed_color_palette,
    guide = "none"
)

result_gene_score <- plot_gene_score_distribution_all_cohorts(cohorts, bank)
gene_score_plots <- result_gene_score$plots
gene_score_recap <- result_gene_score$recap

ggsave("article/figures/Figure_3/Figure_3A_HPE_zscore_distribution.svg", plot = gene_score_plots[["HPE"]], width = 8, height = 12)
ggsave("article/figures/Figure_3/Figure_3A_IOP_zscore_distribution.svg", plot = gene_score_plots[["IOP"]], width = 8, height = 6)

write_csv(gene_score_recap, "article/tables/Figure_3/Figure_3A_gene_score_recap.csv")
write_single_xlsx(gene_score_recap, "article/tables/Figure_3/Figure_3A_gene_score_recap.xlsx", "gene_score_recap")


# ============================================================
# 5. Z-score vs degree centrality, colored by mean correlation
# ============================================================
plot_rank_vs_degree <- function(df, title_suffix, seed_limits) {
    cor_test <- cor.test(df$degree, df$rank, method = "pearson")
    cor_label <- paste0("Pearson r = ", round(cor_test$estimate, 2), ", p = ", signif(cor_test$p.value, 3))

    df_no_seed <- dplyr::filter(df, n_seeds == 0)
    df_seed <- dplyr::filter(df, n_seeds > 0)

    ggplot(df, aes(x = degree, y = rank)) +
        geom_point(data = df_no_seed, color = "grey70", alpha = 0.85, size = 2.5) +
        geom_point(data = df_seed, aes(color = n_seeds), alpha = 0.85, size = 2.5) +
        geom_smooth(method = "lm", se = TRUE, linewidth = 0.6, color = "black", alpha = 0.15) +
        ggrepel::geom_text_repel(
            data = df_seed,
            aes(label = gene), size = 2.8, max.overlaps = Inf, seed = 42, color = "black",
            fontface = "italic"
        ) +
        scale_color_gradientn(name = "Seed count", colours = seed_color_palette, limits = seed_limits) +
        scale_y_reverse() +
        labs(
            x = "Degree centrality", y = "Rank (by RWR z-score)"
        ) +
        custom_theme() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
            plot.subtitle = element_text(hjust = 0.5, size = 9),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "white", color = "black", linewidth = 1)
        )
}

all_dfs <- imap(cohorts, function(co, cname) {
    df <- co$group_stats %>%
        filter(signal_final) %>%
        dplyr::select(gene, z_score) %>%
        left_join(centrality_df, by = "gene")
    df$n_seeds <- sapply(df$gene, function(gn) sum(sapply(co$seeds, function(s) gn %in% s)))
    df$rank <- rank(-df$z_score, ties.method = "min")
    df
})

seed_limits <- range(unlist(map(all_dfs, ~ .x$n_seeds)), na.rm = TRUE)

rank_vs_degree_recap <- imap(all_dfs, function(df, cname) {
    df %>%
        mutate(cohort = cname) %>%
        dplyr::select(cohort, gene, rank, z_score, degree, mean_cor, n_seeds) %>%
        arrange(rank)
}) %>%
    bind_rows()

write_csv(rank_vs_degree_recap, "article/tables/Figure_3/Figure_3X_rank_vs_degree_recap.csv")
write_single_xlsx(rank_vs_degree_recap, "article/tables/Figure_3/Figure_3X_rank_vs_degree_recap.xlsx", "rank_vs_degree_recap")

for (cname in names(cohorts)) {
    co <- cohorts[[cname]]
    df <- all_dfs[[cname]]

    if (nrow(df) == 0) {
        message(cname, ": no signal gene, section 5 skipped")
        next
    }

    p <- plot_rank_vs_degree(df, co$label, seed_limits)
    ggsave(paste0("article/figures/Figure_3/Figure_3X_", cname, "_signal_genes_zscore_vs_degree.svg"), plot = p, width = 8, height = 6, dpi = 600)
}

# ============================================================
# 6. ORA (Reactome / GO:BP / HPO), top 25 terms per cohort
# ------------------------------------------------------------
# Uses plot_ora_dotplot() directly (RECENT_webgestalt_ORA.R):
# x axis = enrichment ratio, size = enrichment ratio, color = FDR.
# The PLOT stays limited to the top 25 (readability), but the
# exported TABLE contains ALL enriched terms (no truncation), so
# the biologist can explore the full result set beyond what is
# visually displayed. Each table is exported as both .csv and
# .xlsx.
# If no term is enriched (common for IOP, few signal genes) or if
# the WebGestalt call fails, the plot AND both table formats are
# skipped for that cohort/database combination -> no file
# produced, just a message.
# ============================================================

ora_dbs <- list(
    Reactome = run_ora_reactome,
    `GO:BP`  = run_ora_gobp,
    HPO      = run_ora_hpo
)


for (cname in names(cohorts)) {
    co <- cohorts[[cname]]
    if (length(co$signal_genes) == 0) {
        message(cname, ": no signal gene, section 6 skipped")
        next
    }

    for (db_label in names(ora_dbs)) {
        ora_result <- tryCatch(ora_dbs[[db_label]](co$signal_genes), error = function(e) NULL)

        if (is.null(ora_result) || nrow(ora_result) == 0) {
            message("No enriched term (", db_label, ", ", co$label, ") -> plot/table skipped")
            next
        }

        db_slug <- gsub("[: ]", "_", db_label)

        p <- plot_ora_dotplot(
            ora_result,
            n     = 20,
            title = db_label
        )
        ggsave(paste0("article/figures/Figure_3/Figure_ORA_", db_slug, "_", cname, "_top20.svg"), plot = p, width = 7, height = 4)

        # Full, untruncated table, sorted by decreasing enrichment ratio
        table_full <- ora_result %>%
            dplyr::arrange(dplyr::desc(enrichmentRatio))
        write.csv(table_full, paste0("article/tables/Figure_3/ORA_", db_slug, "_", cname, "_full.csv"), row.names = FALSE)
        write_single_xlsx(
            table_full,
            paste0("article/tables/Figure_3/ORA_", db_slug, "_", cname, "_full.xlsx"),
            paste0(db_slug, "_", cname)
        )
    }
}
