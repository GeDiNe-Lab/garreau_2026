# Here are all necessary packages (without dependancies)
library(ggrepel)
library(DEGreport)
library(ggsignif)
library(patchwork)
library(grid)
library(ggplot2)
library(WGCNA)
library(tibble)
library(paletteer)
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(DESeq2)
library(Seurat)
library(colorblindr)
library(ggpubr)
library(RANKS)
library(NetPreProc)

RWR <- function(W, ind.positives, gamma = 0.6, tmax = 1000, eps = 1e-10, norm = TRUE, seed.weights = NULL) {
    #' Random Walk with Restart (RWR)
    #'
    #' Performs a random walk with restart algorithm on a weighted network to compute
    #' node proximity scores starting from a set of positive seed nodes.
    #'
    #' @param W A numeric matrix representing the adjacency matrix of the network.
    #' @param ind.positives An integer vector of indices representing the positive seed nodes.
    #' @param gamma A numeric value (0 to 1) controlling the restart probability.
    #'   Default is 0.6. Higher values increase the weight of returning to seed nodes.
    #' @param tmax An integer specifying the maximum number of iterations. Default is 1000.
    #' @param eps A numeric value specifying the convergence threshold. Default is 1e-10.
    #'   The algorithm stops when the L1 norm between consecutive iterations falls below this value.
    #' @param norm A logical value indicating whether to normalize the adjacency matrix
    #'   using the degree-based probability normalization. Default is TRUE.
    #' @param seed.weights A numeric vector of optional weights for the seed nodes.
    #'   If NULL (default), uniform weights are assigned to all positive nodes.
    #'   Must have the same length as ind.positives and contain non-negative values.
    #'
    #' @return A list containing:
    #'   \item{p}{A numeric vector of proximity scores for all nodes, named by node names.}
    #'   \item{ind.positives}{The input indices of positive seed nodes.}
    #'   \item{n.iter}{The number of iterations performed before convergence.}
    #'
    #' @details
    #' The algorithm iteratively computes node proximity scores using the formula:
    #' p(t+1) = (1-gamma) * M * p(t) + gamma * p0, where M is the transition matrix,
    #' p0 is the initial seed distribution, and gamma is the restart probability.
    #'
    #' @examples
    #' \dontrun{
    #' # Example: Run RWR on a sample network
    #' result <- RWR(W = adjacency_matrix, ind.positives = c(1, 5, 10))
    #' }
    if (norm) {
        M <- Prob.norm(W)
    } # M = D^-1 * W
    else {
        M <- W
    }
    n <- nrow(M)
    p0 <- p <- numeric(n)
    names(p) <- names(p0) <- rownames(W)
    rm(W)
    n.positives <- length(ind.positives)
    if (n.positives == 0) {
        stop("RWR: number of core positives is equal to 0!")
    }

    #  original uniform version: set p0 for the positives to 1/n.positives
    # p0[ind.positives] <- 1 / n.positives

    #  weighted version: if seed.weights is provided, use it to set p0 for the positives
    if (is.null(seed.weights)) {
        p0[ind.positives] <- 1 / n.positives # original uniform
    } else {
        if (length(seed.weights) != n.positives) {
            stop("seed.weights must have the same length as ind.positives")
        }
        if (any(seed.weights < 0)) {
            stop("seed.weights must be non-negative")
        }
        p0[ind.positives] <- seed.weights / sum(seed.weights) # normalize to sum to 1
    }
    # end weighted version

    # M <- t(M);
    p <- p0
    for (t in 1:tmax) {
        pold <- p
        p <- ((1 - gamma) * as.vector(pold %*% M)) + gamma * p0
        if (norm1(p - pold) < eps) break()
    }
    return(list(p = p, ind.positives = ind.positives, n.iter = t))
}


kinetic_lineplots <- function(data) {
    #' Generate Kinetic Line Plots
    #'
    #' This function creates a kinetic line plot using ggplot2 for the given data.
    #' The plot displays the mean expression levels of genes over time with error bars
    #' representing the standard deviation. Gene labels are added to the plot for the
    #' data points corresponding to "day12".
    #'
    #' @param data A data frame containing the following columns:
    #'   \itemize{
    #'     \item \code{day}: A factor or character vector representing the time points.
    #'     \item \code{mean_expression}: A numeric vector representing the mean expression levels.
    #'     \item \code{sd_expression}: A numeric vector representing the standard deviation of the expression levels.
    #'     \item \code{gene}: A factor or character vector representing the gene names.
    #'   }
    #'
    #' @return A ggplot object representing the kinetic line plot.
    #'
    #' @import ggplot2
    #' @import ggrepel
    #' @import dplyr
    #'
    #' @examples
    #' \dontrun{
    #' data <- data.frame(
    #'   day = rep(c("day0", "day02", "day04", "day06", "day08", "day10", "day12"), each = 3),
    #'   mean_expression = runif(21, 0, 5),
    #'   sd_expression = runif(21, 0, 1),
    #'   gene = rep(c("gene1", "gene2", "gene3"), 7)
    #' )
    #' plot <- kinetic_lineplots(data)
    #' print(plot)
    #' }
    #'
    plot <- ggplot(data, aes(x = day, y = mean_expression, color = gene, group = gene, label = gene)) +
        geom_point(size = 2) +
        geom_line(size = 2) +
        geom_errorbar(
            aes(
                ymin = mean_expression - sd_expression,
                ymax = mean_expression + sd_expression
            ),
            width = 0, # Width of the horizontal bar on the error bar
            size = 1 # Thickness of the error bars
        ) +
        geom_text_repel(
            data = filter(data, day == "day12"),
            aes(x = 7.75),
            nudge_x = 0.25,
            hjust = 0,
            direction = "y",
            size = 10,
            fontface = "italic",
            segment.color = NA,
            max.overlaps = Inf
        ) +
        scale_x_discrete(breaks = c("day0", "day02", "day04", "day06", "day08", "day10", "day12"), expand = c(0, 0)) +
        ylim(0, 5) +
        ylab("Scaled normalized expression") +
        custom_theme() +
        theme(
            legend.position = "none",
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 30),
            axis.text.y = element_text(size = 20),
            plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
        )
    return(plot)
}

png_save <- function(plot, path, width, height, gg = FALSE, res = 250) {
    #' Save a plot as a PNG file
    #'
    #' This function saves a given plot to a specified file path in PNG format.
    #'
    #' @param plot The plot object to be saved.
    #' @param path A character string specifying the file path where the PNG file will be saved.
    #' @param width The width of the PNG file in pixels.
    #' @param height The height of the PNG file in pixels.
    #' @param res The resolution of the PNG file in dots per inch (DPI). Default is 250.
    #'
    #' @return None. The function is called for its side effect of saving the plot to a file.
    #'
    #' @examples
    #' \dontrun{
    #'   library(ggplot2)
    #'   p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()
    #'   png_save(p, "plot.png", width = 800, height = 600)
    #' }
    if (gg == TRUE) {
        png(filename = path, width = width, height = height, res = res, bg = "transparent")
        print(plot)
        dev.off()
    } else {
        png(filename = path, width = width, height = height, res = res, bg = "transparent")
        draw(plot)
        dev.off()
    }
}

pca_anova <- function(pca_data, metadata, covariates) {
    #' Perform ANOVA on Principal Component Analysis (PCA) Data
    #'
    #' This function computes the correlation and ANOVA between the first five principal components (PCs)
    #' of the PCA data and specified covariates from the metadata.
    #'
    #' @param pca_data A data frame or matrix containing the PCA results, with samples in rows and PCs in columns.
    #' @param metadata A data frame containing the metadata for the samples, with samples in rows and covariates in columns.
    #' @param covariates A character vector specifying the names of the covariates in the metadata to be included in the analysis.
    #'
    #' @return A matrix containing the p-values from the ANOVA tests between each of the first five PCs and the specified covariates.
    #'
    #' @details The function first constructs a combined data frame with the first five PCs and the specified covariates.
    #' The covariates are converted to numeric factors. It then computes the absolute correlation between the PCs and the covariates.
    #' Finally, it performs ANOVA tests to assess the relationship between each PC and each covariate, returning the p-values of these tests.
    #'
    #' @import dplyr
    #' @importFrom stats cor aov summary
    #' @export
    # Building dataframe with first 5 PC and covariates
    PC_covariate <- cbind(pca_data[, 1:5], metadata %>%
        dplyr::select(all_of(covariates)) %>%
        apply(2, function(x) {
            return(as.numeric(factor(x)) - 1)
        }) %>%
        as.matrix())

    # Computing PC-covariate correlation and ANOVA
    PC_covariate_cor <- cor(PC_covariate[, 1:5], PC_covariate[, 6:ncol(PC_covariate)]) %>% abs()
    PC_covariate_ANOVA <- c(6:ncol(PC_covariate)) %>% lapply(function(i) {
        apply(PC_covariate[, 1:5], 2, function(x) {
            aov(x ~ PC_covariate[, i])
        }) %>% sapply(function(x) {
            summary(x)[[1]]$`Pr(>F)`[1]
        })
    })
    PC_covariate_ANOVA <- Reduce(cbind, PC_covariate_ANOVA)
    colnames(PC_covariate_ANOVA) <- colnames(PC_covariate)[6:ncol(PC_covariate)]

    return(PC_covariate_ANOVA)
}


plot_go_term <- function(genelist, path, range = c(1:20), cut = 40, textsize = 20, imgw = 17, imgh = 10) {
    #' Plot GO Term Enrichment
    #'
    #' This function generates a bar plot of Gene Ontology (GO) term enrichment results.
    #'
    #' @param genelist A vector of gene identifiers.
    #' @param path A string specifying the file path to save the plot.
    #' @param range A numeric vector specifying the range of GO terms to display in the plot. Default is c(1:20).
    #' @param cut An integer specifying the width to wrap the GO term descriptions. Default is 40.
    #' @param textsize An integer specifying the size of the text in the plot. Default is 20.
    #' @param imgw A numeric value specifying the width of the saved plot image. Default is 17.
    #' @param imgh A numeric value specifying the height of the saved plot image. Default is 10.
    #'
    #' @return A data frame containing the GO enrichment results.
    #' @export
    #'
    #' @examples
    #' \dontrun{
    #' genelist <- c("ENSG00000141510", "ENSG00000171862", "ENSG00000139618")
    #' plot_go_term(genelist, "output/go_plot")
    #' }

    # performing GO enrichment
    GO_enrichment <- clusterProfiler::enrichGO(genelist,
        OrgDb = "org.Hs.eg.db",
        keyType = "ENSEMBL",
        ont = "BP"
    )

    #  formating GO results
    GO_results <- GO_enrichment@result
    GO_results$GeneRatio <- sapply(GO_enrichment@result$GeneRatio, function(x) {
        eval(parse(text = x))
    }) %>% unname()
    GO_results$rank <- rank(-GO_results$GeneRatio, ties.method = "first")
    GO_results_f <- GO_results[order(GO_results$GeneRatio, decreasing = TRUE)[range], ]
    GO_results_f$Description <- str_wrap(GO_results_f$Description, width = cut) %>% str_to_upper()
    GO_results_f$Description <- factor(GO_results_f$Description, levels = rev(GO_results_f$Description))

    # plottinh
    goplot <- ggplot(GO_results_f, aes(x = GeneRatio, y = reorder(Description, GeneRatio), fill = p.adjust)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = Description),
            hjust = 1.01, # Move text inside the bar, adjust as needed
            color = "#000000", # Make the text white for better visibility
            size = 13
        ) + # Adjust size to fit the text inside the bar
        custom_theme() +
        scale_fill_gradient(name = "p-value\n", low = "#e06663", high = "#327eba") +
        theme(
            axis.title.x = element_text(size = 30), # Adjusts the x-axis title size
            axis.text.x = element_text(size = textsize),
            axis.text.y = element_blank(), # Remove y-axis text
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            legend.text = element_text(size = 30), # Adjusts the legend text size
            legend.title = element_text(size = 40), # Adjusts the legend title size
            legend.key.size = unit(3, "lines")
        )
    ggsave(paste0(path, ".png"), goplot, width = imgw, height = imgh)
    return(GO_results)
}


MyDegPlotCluster <- function(table, time, sign_comp, cluster_i, color = NULL,
                             min_genes = 10,
                             process = FALSE,
                             cluster_column = "cluster",
                             prefix_title = "Group: ") {
    #' Plot Differential Gene Expression for a Specific Cluster
    #' adapted from DegPlotCluster function from DEGReport package
    #' This function generates a plot to visualize the differential gene expression for a specified cluster.
    #'
    #' @param table A data frame containing gene expression data.
    #' @param time A string representing the column name for time points in the data frame.
    #' @param sign_comp A list of comparisons for significance testing.
    #' @param cluster_i An integer specifying the cluster index to be plotted.
    #' @param color An optional string representing the column name for color grouping. Default is NULL.
    #' @param min_genes An integer specifying the minimum number of genes required in a cluster. Default is 10.
    #' @param process A logical value indicating whether to process the data before plotting. Default is FALSE.
    #' @param cluster_column A string representing the column name for clusters in the data frame. Default is "cluster".
    #' @param prefix_title A string to prefix the title of each cluster plot. Default is "Group: ".
    #'
    #' @return A ggplot object representing the differential gene expression plot for the specified cluster.
    #'
    #' @examples
    #' \dontrun{
    #' MyDegPlotCluster(table = my_data, time = "timepoint", sign_comp = list(c("A", "B")), cluster_i = 1)
    #' }
    #'
    #' @import ggplot2
    #' @import dplyr
    #' @importFrom stringi stri_extract_first
    #' @importFrom ggpubr geom_signif
    #' @importFrom stats poly
    #' @importFrom magrittr %>%
    #' @export
    stopifnot(class(table)[1] == "data.frame")

    if (cluster_column %in% colnames(table)) {
        table[["cluster"]] <- table[[cluster_column]]
    }
    if (process) {
        table <- .process(table, time, color)
        print("huu ?")
    }

    if ("cluster" %in% colnames(table)) {
        counts <- table(distinct(table, genes, cluster)[["cluster"]])
        counts <- counts[counts >= min_genes]
        if (length(counts) == 0) {
            stop("No clusters with min_genes > ", min_genes)
        }
        table <- inner_join(table,
            data.frame(
                cluster = as.integer(names(counts)),
                title = paste(
                    prefix_title,
                    names(counts),
                    "- genes:",
                    counts
                ),
                stringsAsFactors = FALSE
            ),
            by = "cluster"
        )
    }

    if (is.null(color)) {
        color <- "dummy"
        table[[color]] <- ""
        lines <- FALSE
    }
    table[["line_group"]] <- paste(
        table[["genes"]],
        table[[color]]
    )
    splan <- length(unique(table[[time]])) - 1L
    table$title <- table$title %>% as.factor()
    old <- table$title %>%
        levels() %>%
        stringi::stri_extract_first(., regex = "\\d+") %>%
        as.integer()
    index <- sort(old, index.return = TRUE)[[2]]
    table$title <- factor(table$title, levels = levels(table$title)[index])

    plot <- ggplot(table, aes_string(
        x = time, y = "value"
    )) +
        geom_boxplot(
            alpha = 0,
            outlier.size = 0,
            outlier.shape = NA,
            color = "#ff9718", # Fixed color for boxplot outline
            fill = "#ff9718"
        ) +
        geom_jitter(
            alpha = 0.4, size = 1,
            color = "#ff9718",
            fill = "#ff9718"
        ) +
        stat_smooth(
            aes_string(
                x = time, y = "value",
                group = color, color = color
            ),
            se = FALSE,
            method = "lm", formula = y ~ poly(x, splan),
            color = "black"
        ) +
        geom_line(aes_string(group = "line_group"), alpha = 0.1, color = "#ff9718") +
        geom_signif(
            comparisons = sign_comp, # Groups being compared
            map_signif_level = TRUE, # Automatically converts p-values to significance stars
            y_position = rep(2, length(sign_comp)),
            color = "black",
            size = 1,
            step_increase = 0.05,
            textsize = 10
        ) +
        ggtitle(paste0("Cluster ", cluster_i, ", ", length(unique(table$genes)), " genes")) +
        ylab("Z-score of gene abundance") +
        xlab("") +
        ylim(-2.1, 3.2) +
        custom_theme(hide_legend = TRUE) +
        theme(
            plot.title = element_text(size = 40, face = "bold"),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 30)
        )
    plot
}

plotPCA.DESeqTransform <- function(object, intgroup = "condition",
                                   ntop = 500, returnData = FALSE, pcsToUse = 1:15, batch = NULL) {
    #' Plot PCA for DESeqTransform Object
    #' adapted from DESeq2 pca function
    #' This function performs Principal Component Analysis (PCA) on a DESeqTransform object and returns the PCA results.
    #'
    #' @param object A DESeqTransform object.
    #' @param intgroup A character vector specifying the columns of colData(object) to use for grouping.
    #' @param ntop An integer specifying the number of top features by variance to use for PCA. Default is 500.
    #' @param returnData A logical indicating whether to return the PCA data. Default is FALSE.
    #' @param pcsToUse A numeric vector specifying which principal components to use. Default is 1:15.
    #' @param batch A parameter for batch effect correction (currently not used).
    #'
    #' @return A data frame containing the PCA results, with additional attributes:
    #' \itemize{
    #'   \item \code{percentVar}: The contribution to the total variance for each component.
    #'   \item \code{factoextra}: PCA variable results from the factoextra package.
    #'   \item \code{pca_var}: PCA results from the prcomp function.
    #' }
    #'
    #' @examples
    #' \dontrun{
    #'   dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ condition)
    #'   vsd <- vst(dds)
    #'   pca_data <- plotPCA.DESeqTransform(vsd, intgroup = "condition")
    #' }
    #'
    #' @importFrom MatrixGenerics rowVars
    #' @importFrom factoextra get_pca_var
    #' @importFrom stats prcomp
    #' @importFrom utils head
    #' @export
    message(paste0("using ntop=", ntop, " top features by variance"))

    # calculate the variance for each gene

    rv <- MatrixGenerics::rowVars(assay(object))

    # select the ntop genes by variance
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

    # perform a PCA on the data in assay(x) for the selected genes
    pca <- prcomp(t(assay(object)[select, ]))

    pca_data_fe <- factoextra::get_pca_var(pca)

    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum(pca$sdev^2)

    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }

    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])

    # add the intgroup factors together to create a new grouping factor
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = ":"))
    } else {
        colData(object)[[intgroup]]
    }

    # assembly the data
    pcs <- paste0("PC", pcsToUse)
    print(pcsToUse)

    d <- pcsToUse %>% lapply(function(PC) {
        return(pca$x[, PC])
    })
    d <- do.call(cbind, d)
    d <- cbind(d, data.frame(group = group, intgroup.df, name = colnames(object)))
    colnames(d)[pcsToUse] <- pcs

    attr(d, "percentVar") <- percentVar[pcsToUse]
    attr(d, "factoextra") <- pca_data_fe
    attr(d, "pca_var") <- prcomp(assay(object)[select, ])
    return(d)
}

custom_theme <- function(diag_text = FALSE, hide_legend = FALSE, hide_x_lab = FALSE) {
    #' Custom ggplot2 Theme
    #'
    #' This function creates a custom ggplot2 theme with options to hide the legend,
    #' rotate diagonal text, and hide the x-axis label.
    #'
    #' @param diag_text Logical. If TRUE, the x-axis text will be rotated by 45 degrees. Default is FALSE.
    #' @param hide_legend Logical. If TRUE, the legend will be hidden. Default is FALSE.
    #' @param hide_x_lab Logical. If TRUE, the x-axis label will be hidden. Default is FALSE.
    #'
    #' @return A ggplot2 theme object with the specified customizations.
    #'
    #' @examples
    #' library(ggplot2)
    #' ggplot(mtcars, aes(x = wt, y = mpg)) +
    #'   geom_point() +
    #'   custom_theme(diag_text = TRUE, hide_legend = TRUE, hide_x_lab = TRUE)
    if (hide_legend == TRUE) {
        hide_legend <- "none"
    } else {
        hide_legend <- "right"
    }
    if (diag_text == TRUE) {
        text_angle <- 45
    } else {
        text_angle <- 0
    }
    if (hide_x_lab == TRUE) {
        axis_text_x <- element_blank()
    } else {
        axis_text_x <- element_text(size = 12, angle = text_angle, hjust = 1)
    }
    return(theme(
        plot.title = element_text(size = 20),
        axis.text.x = axis_text_x,
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "#dcdcdc"),
        axis.line = element_line(size = 0.2, colour = "Black"),
        axis.title = element_text(size = 20),
        legend.position = hide_legend
    ))
}

gene_converter <- function(gene_vec, from, to, species = "human") {
    #' Convert Gene Identifiers
    #'
    #' This function converts a vector of gene identifiers from one type to another for a specified species.
    #'
    #' @param gene_vec A character vector of gene identifiers to be converted.
    #' @param from A character string specifying the type of the input gene identifiers (e.g., "ENSEMBL", "SYMBOL").
    #' @param to A character string specifying the type of the output gene identifiers (e.g., "SYMBOL", "ENTREZID").
    #' @param species A character string specifying the species of the genes. Default is "human". Other option is "mouse".
    #'
    #' @return A character vector of converted gene identifiers.
    #'
    #' @details
    #' The function uses the `org.Hs.eg.db` database for human genes and the `org.Mm.eg.db` database for mouse genes.
    #' If the input gene vector contains the gene "ENSG00000229415", it will be manually converted to "SFTA3".
    #'
    #' @examples
    #' \dontrun{
    #' gene_vec <- c("ENSG00000229415", "ENSG00000139618")
    #' converted_genes <- gene_converter(gene_vec, from = "ENSEMBL", to = "SYMBOL", species = "human")
    #' }
    #'
    #' @importFrom AnnotationDbi mapIds
    #' @importFrom magrittr %>%
    #' @import org.Hs.eg.db
    #' @import org.Mm.eg.db
    #' @export
    if (species == "human") {
        db <- org.Hs.eg.db
    } else if (species == "mouse") {
        db <- org.Mm.eg.db
    }
    SFTA3_i <- which(gene_vec == "ENSG00000257520")
    sym <- mapIds(db,
        key = gene_vec,
        keytype = from, column = to
    ) %>% unname()
    if (!is.null(SFTA3_i)) {
        sym[SFTA3_i] <- "SFTA3"
    }
    return(sym)
}

readcounts <- function(file, sep = ",", header = T, row.names = 1) {
    #' Read Counts from a File
    #'
    #' This function reads a file and returns its contents as a matrix.
    #'
    #' @param file A character string specifying the path to the file to be read.
    #' @param sep A character string specifying the field separator character. Default is ",".
    #' @param header A logical value indicating whether the file contains the names of the variables as its first line. Default is TRUE.
    #' @param row.names A specification of the column to be used as row names. Default is 1.
    #'
    #' @return A matrix containing the data from the file.
    #' @export
    #'
    #' @examples
    #' \dontrun{
    #'   data_matrix <- readcounts("data.csv")
    #' }
    return(as.matrix(read.table(file, sep = sep, header = header, row.names = row.names)))
}

cleaning_pLI <- function(pLI_vec, agg_fun = max) {
    vapply(pLI_vec, function(x) {
        if (is.na(x) || !nzchar(x)) {
            return(0)
        }
        vals <- str_split(x, ";")[[1]]
        vals <- vals[nzchar(vals)]
        vals <- suppressWarnings(as.numeric(vals))
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) {
            return(0)
        }
        agg_fun(vals)
    }, FUN.VALUE = numeric(1), USE.NAMES = FALSE)
}
