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
plot_ora_dotplot <- function(ora_results, n = 20,
                             term_col = "description",
                             ratio_col = "enrichmentRatio",
                             fdr_col = "FDR",
                             title = NULL) {
    if (is.null(ora_results) || nrow(ora_results) == 0) {
        stop("ora_results is empty or NULL, nothing to plot")
    }

    df <- ora_results %>%
        dplyr::arrange(dplyr::desc(.data[[ratio_col]])) %>%
        dplyr::slice_head(n = n) %>%
        dplyr::mutate(
            !!term_col := factor(.data[[term_col]],
                levels = rev(unique(.data[[term_col]]))
            )
        )

    p <- ggplot(df, aes(
        x = .data[[ratio_col]],
        y = .data[[term_col]],
        size = .data[[ratio_col]],
        color = .data[[fdr_col]]
    )) +
        geom_point() +
        scale_color_gradient(low = "#D7263D", high = "#3B8BD4", name = "FDR") +
        scale_size_continuous(name = "Enrichment\nratio") +
        labs(
            x = "Enrichment ratio",
            y = NULL,
            title = title %||% paste0("Top ", nrow(df), " enriched terms")
        ) +
        custom_theme() +
        theme(
            axis.text.y = element_text(size = 9),
            panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold")
        )

    return(p)
}
