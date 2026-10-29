# Loading packages and functions
library(Matrix)
library(tidyverse)
library(tibble)
library(ggplot2)
library(cluster)
library(WebGestaltR)

# Loading custom functions
source("article/scripts/custom_fct.R")
run_ora_hpo <- function(genelist, referenceGene = NULL, referenceSet = "genome_protein-coding",
                        interestGeneType = "genesymbol", fdrThr = 0.05, minNum = 9, maxNum = 500) {
    #' ORA enrichment on Human Phenotype Ontology (HPO)
    #'
    #' @return A data frame containing the HPO enrichment results.
    #' @export

    ora_args <- list(
        enrichMethod = "ORA",
        organism = "hsapiens",
        enrichDatabase = "phenotype_Human_Phenotype_Ontology",
        enrichDatabaseType = "genesymbol",
        interestGene = genelist,
        interestGeneType = interestGeneType,
        minNum = minNum,
        maxNum = maxNum,
        sigMethod = "fdr",
        fdrMethod = "BH",
        fdrThr = fdrThr,
        isOutput = FALSE,
        useWeightedSetCover = FALSE,
        usekMedoid = FALSE
    )

    if (!is.null(referenceGene)) {
        ora_args$referenceGene <- referenceGene
        ora_args$referenceGeneType <- interestGeneType
    } else {
        ora_args$referenceSet <- referenceSet
    }

    hpo_results <- do.call(WebGestaltR, ora_args)

    return(hpo_results)
}


run_ora_gobp <- function(genelist, referenceGene = NULL, referenceSet = "genome_protein-coding",
                         interestGeneType = "genesymbol", fdrThr = 0.05, minNum = 9, maxNum = 500,
                         noRedundant = TRUE) {
    #' ORA enrichment on Gene Ontology - Biological Process (GO-BP)
    #'
    #' @return A data frame containing the GO-BP enrichment results.
    #' @export

    enrich_db <- if (noRedundant) {
        "geneontology_Biological_Process_noRedundant"
    } else {
        "geneontology_Biological_Process"
    }

    ora_args <- list(
        enrichMethod = "ORA",
        organism = "hsapiens",
        enrichDatabase = enrich_db,
        enrichDatabaseType = "genesymbol",
        interestGene = genelist,
        interestGeneType = interestGeneType,
        minNum = minNum,
        maxNum = maxNum,
        sigMethod = "fdr",
        fdrThr = fdrThr,
        isOutput = FALSE
    )

    if (!is.null(referenceGene)) {
        ora_args$referenceGene <- referenceGene
        ora_args$referenceGeneType <- interestGeneType
    } else {
        ora_args$referenceSet <- referenceSet
    }

    gobp_results <- do.call(WebGestaltR, ora_args)

    return(gobp_results)
}


run_ora_reactome <- function(genelist, referenceGene = NULL, referenceSet = "genome_protein-coding",
                             interestGeneType = "genesymbol", fdrThr = 0.05, minNum = 9, maxNum = 500) {
    #' ORA enrichment on REACTOME pathways
    #'
    #' @return A data frame containing the REACTOME enrichment results.
    #' @export

    ora_args <- list(
        enrichMethod = "ORA",
        organism = "hsapiens",
        enrichDatabase = "pathway_Reactome",
        enrichDatabaseType = "genesymbol",
        interestGene = genelist,
        interestGeneType = interestGeneType,
        minNum = minNum,
        maxNum = maxNum,
        sigMethod = "fdr",
        fdrThr = fdrThr,
        isOutput = FALSE
    )

    if (!is.null(referenceGene)) {
        ora_args$referenceGene <- referenceGene
        ora_args$referenceGeneType <- interestGeneType
    } else {
        ora_args$referenceSet <- referenceSet
    }

    reactome_results <- do.call(WebGestaltR, ora_args)

    return(reactome_results)
}


# ============================================================
# plot_ora_dotplot()
# ------------------------------------------------------------
# Dot plot of the top N most enriched terms (by enrichmentRatio) from
# an ORA result (run_ora_hpo / run_ora_gobp / run_ora_reactome).
#
#   - y axis    : term (description)
#   - x axis    : enrichment ratio
#   - size      : enrichment ratio
#   - color     : FDR
#
# @param ora_results   data.frame returned by run_ora_hpo/gobp/reactome
#                      (standard WebGestaltR columns: description,
#                      enrichmentRatio, FDR)
# @param n              number of terms to display (top n by
#                      enrichmentRatio, default 20)
# @param term_col       name of the term column (default "description")
# @param ratio_col      name of the enrichment ratio column (default "enrichmentRatio")
# @param fdr_col        name of the FDR column (default "FDR")
# @param title           plot title (optional)
#
# @return a ggplot object
# ============================================================
# plot_ora_dotplot <- function(ora_results, n = 20,
#                              term_col = "description",
#                              ratio_col = "enrichmentRatio",
#                              fdr_col = "FDR",
#                              title = NULL,
#                              wrap_width = 30) {
#     if (is.null(ora_results) || nrow(ora_results) == 0) {
#         stop("ora_results is empty or NULL, nothing to plot")
#     }

#     df <- ora_results %>%
#         dplyr::arrange(dplyr::desc(.data[[ratio_col]])) %>%
#         dplyr::slice_head(n = n) %>%
#         dplyr::mutate(
#             !!term_col := factor(.data[[term_col]],
#                 levels = rev(unique(.data[[term_col]]))
#             )
#         )

#     p <- ggplot(df, aes(
#         x = .data[[ratio_col]],
#         y = .data[[term_col]],
#         size = .data[[ratio_col]],
#         color = .data[[fdr_col]]
#     )) +
#         geom_point() +
#         scale_color_gradient(low = "#D7263D", high = "#3B8BD4", name = "FDR") +
#         scale_size_continuous(name = "Enrichment\nratio") +
#         scale_y_discrete(labels = scales::label_wrap(wrap_width)) +
#         labs(
#             x = "Enrichment ratio",
#             y = NULL,
#             title = title
#         ) +
#         custom_theme() +
#         theme(
#             axis.text.y = element_text(size = 9),
#             panel.grid.minor = element_blank(),
#             plot.title = element_text(face = "bold")
#         )

#     return(p)
# }
plot_ora_dotplot <- function(ora_results, n = 20,
                             term_col = "description",
                             ratio_col = "enrichmentRatio",
                             fdr_col = "FDR",
                             title = NULL,
                             wrap_width = 30,
                             base_size = 9,
                             size_step = 1,
                             min_size = 6,
                             legend_vpos = 0.35,
                             text_scale = 0.75) {
    if (is.null(ora_results) || nrow(ora_results) == 0) {
        stop("ora_results is empty or NULL, nothing to plot")
    }

    df <- ora_results %>%
        dplyr::arrange(dplyr::desc(.data[[ratio_col]])) %>%
        dplyr::slice_head(n = n) %>%
        dplyr::mutate(
            !!term_col := factor(.data[[term_col]],
                levels = unique(.data[[term_col]])
            )
        )

    term_levels <- levels(df[[term_col]])

    # wrap a label on 2 lines at most:
    # - fits in wrap_width: 1 line
    # - fits in 2 lines of wrap_width: normal wrap
    # - longer: balanced split in 2 lines at a space (lines exceed wrap_width)
    wrap_two_lines <- function(label, width) {
        if (nchar(label) <= width) return(label)
        lines <- strwrap(label, width = width)
        if (length(lines) <= 2) return(paste(lines, collapse = "\n"))
        words <- strsplit(label, " ")[[1]]
        if (length(words) < 2) return(label)
        splits <- seq_len(length(words) - 1)
        longest <- vapply(splits, function(k) {
            max(nchar(paste(words[1:k], collapse = " ")),
                nchar(paste(words[(k + 1):length(words)], collapse = " ")))
        }, numeric(1))
        k <- splits[which.min(longest)]
        paste(paste(words[1:k], collapse = " "),
              paste(words[(k + 1):length(words)], collapse = " "),
              sep = "\n")
    }
    wrapped_labels <- vapply(term_levels, wrap_two_lines, character(1),
                             width = wrap_width, USE.NAMES = FALSE)

    # font size: base_size for 1 line, minus size_step for 2 lines,
    # then reduced in proportion when a line exceeds wrap_width, floored at min_size
    line_list <- strsplit(wrapped_labels, "\n")
    n_lines <- lengths(line_list)
    max_line_len <- vapply(line_list, function(l) max(nchar(l)), numeric(1))
    label_sizes <- base_size - (n_lines - 1) * size_step
    too_long <- max_line_len > wrap_width
    label_sizes[too_long] <- label_sizes[too_long] * wrap_width / max_line_len[too_long]
    label_sizes <- pmax(label_sizes, min_size)

    # alternate black / grey labels
    label_cols <- rep_len(c("black", "grey55"), length(term_levels))

    # reference sizes: custom_theme() sets title and axis titles to 20,
    # legend sizes come from the ggplot2 default theme (base size 11)
    title_size <- 20 * text_scale
    axis_title_size <- 20 * text_scale
    legend_title_size <- 11 * text_scale
    legend_text_size <- 11 * 0.8 * text_scale
    legend_key_size <- 1.2 * text_scale

    p <- ggplot(df, aes(
        x = .data[[term_col]],
        y = .data[[ratio_col]],
        size = .data[[ratio_col]],
        color = .data[[fdr_col]]
    )) +
        geom_point() +
        scale_color_gradient(low = "#D7263D", high = "#3B8BD4", name = "FDR") +
        scale_size_continuous(name = "Enrichment\nratio") +
        guides(
            color = guide_colorbar(order = 1),
            size = guide_legend(order = 2)
        ) +
        scale_x_discrete(labels = setNames(wrapped_labels, term_levels)) +
        labs(
            x = NULL,
            y = "Enrichment ratio",
            title = title
        ) +
        custom_theme() +
        theme(
            axis.text.x = element_text(
                angle = -55, hjust = 0, vjust = 1,
                lineheight = 0.85,
                size = label_sizes,
                colour = label_cols
            ),
            axis.title.y = element_text(size = axis_title_size),
            plot.title = element_text(face = "bold", size = title_size),
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_text_size),
            legend.key.size = unit(legend_key_size, "lines"),
            panel.grid.minor = element_blank(),
            legend.justification = c(0, legend_vpos),
            plot.background = element_rect(fill = "white", color = "black", linewidth = 1)
        )

    return(p)
}