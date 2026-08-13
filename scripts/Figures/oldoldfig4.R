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
HPE_seeds %>%
    unlist() %>%
    unique()
# groupes fonctionnels
#  Base voie de signalisation
BMP_TGFbeta <- c(
    "FBN2", "NR3C1", "COL2A1", "BMP2", "SMAD3",
    "BMP4", "RERE"
)

SHH <- c(
    "NDST1", "SHH", "GLI2", "LRP2", "VCAN", "GLI3", "ZIC2",
    "CREBBP", "PTCH1", "SIX3", "FOXA2",
    "PAX6", "EMX2"
)

NOTCH <- c(
    "CIT", "NOTCH1", "NOTCH2",
    "HES1"
)

FGF <- c(
    "APLP2", "FGF8",
    "FGFR1"
)

WNT_PCP <- c(
    "MACF1", "ITGA5", "CELSR1", "CELSR2", "FAT3", "NCOA6",
    "APC2", "AXIN2", "PTK7", "TCF7L2",
    "WNT5B", "FZD2", "FZD8", "DAAM1"
)
#  bonus groups (to test later maybe)
# SHH_BMP <- c(SHH, "GLIS3", "NDST1", "VCAN")
# WNT <- c("MACF1", "CACHD1", "ITGA5", "LGR5")


# ============================================================
# 1. Groupes de patients par voie de signalisation
# ------------------------------------------------------------
# Un patient appartient au groupe d'une voie s'il possede au moins
# un seed (variant) dans un des genes de cette voie.
# ============================================================

pathway_groups <- list(
    BMP_TGFbeta = BMP_TGFbeta,
    SHH         = SHH,
    NOTCH       = NOTCH,
    FGF         = FGF,
    WNT_PCP     = WNT_PCP
)

get_patients_with_seed_in <- function(gene_set, seeds_list) {
    names(seeds_list)[sapply(seeds_list, function(seeds) any(seeds %in% gene_set))]
}

patient_groups <- lapply(pathway_groups, get_patients_with_seed_in, seeds_list = HPE_seeds)

# Recap taille des groupes
sapply(patient_groups, length)


# ============================================================
# 2. Extraction des genes signal par groupe (process_group)
# ============================================================

pathway_group_results <- lapply(patient_groups, function(patient_ids) {
    if (length(patient_ids) == 0) {
        message("Groupe vide, ignore")
        return(NULL)
    }
    process_group(
        patient_ids      = patient_ids,
        RWR_results      = HPE_results,
        seeds            = HPE_seeds,
        permutation_bank = bank
    )
})

# Genes signal (signal_final == TRUE) pour chaque groupe
pathway_signal_genes <- lapply(pathway_group_results, function(df) {
    if (is.null(df)) {
        return(character(0))
    }
    df$gene[df$signal_final]
})

sapply(pathway_signal_genes, length)

# ============================================================
# 3. ORA (Reactome, HPO, GO-BP) pour chaque groupe
# ============================================================

pathway_ORA_reactome <- lapply(names(pathway_signal_genes), function(pw) {
    genelist <- pathway_signal_genes[[pw]]
    if (length(genelist) < 5) {
        message("Groupe ", pw, " : moins de 5 genes signal, ORA ignoree")
        return(NULL)
    }
    run_ora_reactome(
        genelist         = genelist,
        # referenceGene    = context_genes,
        interestGeneType = "genesymbol",
        outputDirectory  = "article/enrichment",
        projectName      = paste0(pw, "_signal_Reactome"),
        fdrThr           = 0.05,
        minNum           = 5,
        maxNum           = 2000,
        is.output        = TRUE
    )
})
names(pathway_ORA_reactome) <- names(pathway_signal_genes)

pathway_ORA_HPO <- lapply(names(pathway_signal_genes), function(pw) {
    genelist <- pathway_signal_genes[[pw]]
    if (length(genelist) < 5) {
        message("Groupe ", pw, " : moins de 5 genes signal, ORA ignoree")
        return(NULL)
    }
    run_ora_hpo(
        genelist         = genelist,
        # referenceGene    = context_genes,
        interestGeneType = "genesymbol",
        outputDirectory  = "article/enrichment",
        projectName      = paste0(pw, "_signal_HPO"),
        fdrThr           = 0.05,
        minNum           = 5,
        maxNum           = 2000,
        is.output        = TRUE
    )
})
names(pathway_ORA_HPO) <- names(pathway_signal_genes)

pathway_ORA_GOBP <- lapply(names(pathway_signal_genes), function(pw) {
    genelist <- pathway_signal_genes[[pw]]
    if (length(genelist) < 5) {
        message("Groupe ", pw, " : moins de 5 genes signal, ORA ignoree")
        return(NULL)
    }
    run_ora_gobp(
        genelist         = genelist,
        # referenceGene    = context_genes,
        interestGeneType = "genesymbol",
        outputDirectory  = "article/enrichment",
        projectName      = paste0(pw, "_signal_GOBP"),
        fdrThr           = 0.05,
        minNum           = 5,
        maxNum           = 2000,
        noRedundant      = TRUE,
        is.output        = TRUE
    )
})
names(pathway_ORA_GOBP) <- names(pathway_signal_genes)


# ============================================================
# Diagramme de Venn - overlap des patients entre groupes de voies
# ============================================================
# install.packages("ggVennDiagram") si pas deja installe
library(ggVennDiagram)

# On ne garde que les groupes non vides
venn_groups <- patient_groups[sapply(patient_groups, length) > 0]

p_venn_pathways <- ggVennDiagram(
    venn_groups,
    label = "count",
    edge_size = 0.5
) +
    scale_fill_gradient(low = "#F4F4F4", high = "#3B8BD4", name = "N patients") +
    labs(title = "Overlap des patients par voie de signalisation") +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right"
    )

print(p_venn_pathways)

ggsave(
    "article/figures/pathway_patient_overlap_venn.png",
    plot = p_venn_pathways, width = 8, height = 7, dpi = 300
)


# ============================================================
# make_ora_heatmap()
# ------------------------------------------------------------
# Heatmap "termes enrichis x groupes", inspiree de la Figure 4.
#   - x : groupes, dans l'ordre fixe SHH > NOTCH > BMP_TGFbeta >
#         WNT_PCP > FGF (tous les groupes demandes apparaissent,
#         meme ceux sans terme significatif -> colonne blanche)
#   - y : termes enrichis (union des groupes, top n par groupe),
#         labels retournes a la ligne (wrap) a 45 caracteres sans
#         couper un mot
#   - couleur des cases : enrichmentRatio, blanc si absent/NA
#   - cases carrees (meme taille en x et y)
#
# @param ora_list       liste NOMMEE de data.frames (ou NULL) issus
#                       de run_ora_gobp/run_ora_hpo/run_ora_reactome
# @param group_order    ordre voulu des groupes en x
# @param n_per_group    top n termes par groupe (NULL = tous)
# @param wrap_width     largeur max (caracteres) avant retour a la
#                       ligne des labels de termes (defaut 45)
# @param term_col       colonne du terme (defaut "description")
# @param ratio_col      colonne enrichment ratio (defaut "enrichmentRatio")
# @param title          titre du graphique
#
# @return liste(plot = <ggplot>, n_terms = <int>, n_groups = <int>)
#         n_terms/n_groups servent a calculer une taille de sauvegarde
#         (ggsave) qui conserve des cases carrees.
# ============================================================
make_ora_heatmap <- function(ora_list,
                             group_order = c("SHH", "NOTCH", "BMP_TGFbeta", "WNT_PCP", "FGF"),
                             n_per_group = 15,
                             wrap_width = 45,
                             term_col = "description",
                             ratio_col = "enrichmentRatio",
                             title = NULL) {
    if (length(ora_list) == 0 || is.null(names(ora_list))) {
        stop("ora_list doit etre une liste NOMMEE (un data.frame ou NULL par groupe)")
    }

    # Ordre des groupes : uniquement ceux presents dans ora_list, dans l'ordre demande
    all_groups <- group_order[group_order %in% names(ora_list)]
    missing_from_order <- setdiff(names(ora_list), group_order)
    if (length(missing_from_order) > 0) {
        message(
            "Groupe(s) absent(s) de group_order, ajoute(s) a la fin : ",
            paste(missing_from_order, collapse = ", ")
        )
        all_groups <- c(all_groups, missing_from_order)
    }

    # Normalise : NULL -> data.frame vide
    ora_list_clean <- lapply(ora_list, function(df) {
        if (is.null(df) || nrow(df) == 0) {
            return(data.frame())
        }
        df
    })

    empty_groups <- all_groups[sapply(ora_list_clean[all_groups], nrow) == 0]
    non_empty_groups <- setdiff(all_groups, empty_groups)

    if (length(empty_groups) > 0) {
        message(
            "Aucun terme significatif (FDR) pour le(s) groupe(s) : ",
            paste(empty_groups, collapse = ", "),
            " -- colonne(s) laissee(s) blanche(s)"
        )
    }

    if (length(non_empty_groups) == 0) {
        warning("Aucun groupe n'a de terme significatif : heatmap entierement blanche")
        df_full <- tidyr::expand_grid(
            term = "aucun terme significatif",
            group = factor(all_groups, levels = all_groups),
            enrichmentRatio = NA_real_
        )
    } else {
        df_long <- purrr::imap_dfr(
            ora_list_clean[non_empty_groups],
            function(df, group_name) {
                df <- df %>% dplyr::arrange(dplyr::desc(.data[[ratio_col]]))
                if (!is.null(n_per_group)) df <- df %>% dplyr::slice_head(n = n_per_group)
                df %>%
                    dplyr::transmute(
                        term            = .data[[term_col]],
                        enrichmentRatio = .data[[ratio_col]],
                        group           = group_name
                    )
            }
        )

        term_order <- df_long %>%
            dplyr::group_by(term) %>%
            dplyr::summarise(max_ratio = max(enrichmentRatio, na.rm = TRUE)) %>%
            dplyr::arrange(max_ratio) %>%
            dplyr::pull(term)

        df_long <- df_long %>%
            dplyr::mutate(term = factor(term, levels = term_order))

        df_full <- df_long %>%
            tidyr::complete(term, group = all_groups) %>%
            dplyr::mutate(group = factor(group, levels = all_groups))
    }

    n_terms <- nlevels(df_full$term)
    n_groups <- length(all_groups)

    p <- ggplot(df_full, aes(x = group, y = term, fill = enrichmentRatio)) +
        geom_tile(color = "grey85") +
        scale_fill_gradient(
            low = "#FDE9E9", high = "#7A1F2B",
            na.value = "white",
            name = "Enrichment\nratio"
        ) +
        scale_y_discrete(labels = scales::label_wrap(wrap_width)) +
        coord_fixed(ratio = 1) +
        labs(x = NULL, y = NULL, title = title) +
        custom_theme() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
            axis.text.y = element_text(size = 8, lineheight = 0.85),
            panel.grid = element_blank(),
            plot.title = element_text(face = "bold", hjust = 0.5)
        )

    list(plot = p, n_terms = n_terms, n_groups = n_groups)
}


# ============================================================
# save_square_heatmap()
# ------------------------------------------------------------
# Sauvegarde un plot produit par make_ora_heatmap() en calculant
# une taille (largeur/hauteur) proportionnelle a n_groups/n_terms,
# pour que les cases restent visuellement carrees malgre
# coord_fixed(ratio = 1). unit_size = taille (pouces) d'une case ;
# margin_w/margin_h = marge reservee aux axes/legende/titre (labels
# retournes a la ligne -> marge gauche un peu plus large).
# ============================================================
save_square_heatmap <- function(heatmap_obj, filename,
                                unit_size = 0.35,
                                margin_w = 5, margin_h = 1.5,
                                dpi = 300) {
    w <- heatmap_obj$n_groups * unit_size + margin_w
    h <- heatmap_obj$n_terms * unit_size + margin_h

    ggsave(filename, plot = heatmap_obj$plot, width = w, height = h, dpi = dpi, limitsize = FALSE)
}


# ============================================================
# Application aux 3 bases (GO-BP, HPO, Reactome)
# ============================================================

hm_GOBP <- make_ora_heatmap(
    pathway_ORA_GOBP,
    title = "GO Biological Process"
)
print(hm_GOBP$plot)

hm_HPO <- make_ora_heatmap(
    pathway_ORA_HPO,
    title = "HPO Phenotypes"
)
print(hm_HPO$plot)

hm_Reactome <- make_ora_heatmap(
    pathway_ORA_reactome,
    title = "Reactome Pathways"
)
print(hm_Reactome$plot)


# ============================================================
# Sauvegarde
# ============================================================
dir.create("article/figures", recursive = TRUE, showWarnings = FALSE)

save_square_heatmap(hm_GOBP, "article/figures/heatmap_GOBP_by_group.png")
save_square_heatmap(hm_HPO, "article/figures/heatmap_HPO_by_group.png")
save_square_heatmap(hm_Reactome, "article/figures/heatmap_Reactome_by_group.png")


# ============================================================
# Figure 4 - Diagramme de Venn et table de l'overlap des
# genes signal (signal_final == TRUE) entre groupes de voies
# ============================================================

dir.create("article/tables", recursive = TRUE, showWarnings = FALSE)

# On ne garde que les groupes non vides
venn_signal_genes <- pathway_signal_genes[sapply(pathway_signal_genes, length) > 0]

p_venn_signal_genes <- ggVennDiagram(
    venn_signal_genes,
    label = "count",
    edge_size = 0.5
) +
    scale_fill_gradient(low = "#F4F4F4", high = "#3B8BD4", name = "N genes") +
    labs(title = "Overlap des genes signal par voie de signalisation") +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right"
    )

print(p_venn_signal_genes)

ggsave(
    "article/figures/pathway_signal_genes_overlap_venn.png",
    plot = p_venn_signal_genes, width = 8, height = 7, dpi = 300
)


# ------------------------------------------------------------
# Table associee : matrice d'appartenance gene x groupe
# ------------------------------------------------------------

all_signal_genes <- sort(unique(unlist(venn_signal_genes)))

signal_gene_membership <- sapply(venn_signal_genes, function(genes) {
    all_signal_genes %in% genes
})
rownames(signal_gene_membership) <- all_signal_genes
signal_gene_membership <- as.data.frame(signal_gene_membership)

write.csv(
    signal_gene_membership %>% rownames_to_column("gene"),
    "article/tables/signal_genes_membership_by_group.csv",
    row.names = FALSE
)


# ------------------------------------------------------------
# Table recap par combinaison exacte de groupes (type UpSet) :
# pour chaque combinaison observee de voies, nombre de genes
# et liste des genes concernes
# ------------------------------------------------------------

signal_gene_combinations <- signal_gene_membership %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "group", values_to = "present") %>%
    filter(present) %>%
    group_by(gene) %>%
    summarise(combination = paste(sort(group), collapse = " & "), .groups = "drop") %>%
    group_by(combination) %>%
    summarise(
        n_genes = n(),
        genes   = paste(sort(gene), collapse = ", "),
        .groups = "drop"
    ) %>%
    arrange(desc(n_genes))

write.csv(
    signal_gene_combinations,
    "article/tables/signal_genes_overlap_combinations.csv",
    row.names = FALSE
)

print(signal_gene_combinations)

HPE_diag_genes <- c(
    "SHH", "GLI2", "FGF8", "ZIC2", "SIX3", "TGIF1", "PTCH1", "GAS1", "FGFR1", "BOC", "SCUBE2", "LRP2"
)

c(
    "Pituitary dwarfism",
    "Ectopic anterior pituitary gland",
    "Anterior pituitary agenesis",
    "Median cleft lip and palate",
    "Hypoplastic iris stroma",
    "Basal encephalocele",
    "Ethmoidal encephalocele",
    "Abnormal hard palate morphology",
    "Midnasal stenosis",
    "Abnormality of the midnasal cavity",
    "Abnormality of the nasal cavity",
    "Median cleft palate",
    "Cyclopia",
    "Median cleft lip",
    "Abnormal hypothalamus physiology",
    "Expressive language delay",
    "Abnormality of the diencephalon",
    "Cleft lip"
)

# >>> A DEFINIR SOI-MEME : vecteur de termes HPO par groupe
group_hpo_terms <- list(
    BMP_TGFbeta = c(
        "Pituitary dwarfism",
        "Ectopic anterior pituitary gland",
        "Anterior pituitary agenesis",
        "Median cleft lip and palate",
        "Hypoplastic iris stroma",
        "Basal encephalocele",
        "Ethmoidal encephalocele",
        "Abnormal hard palate morphology",
        "Midnasal stenosis",
        "Abnormality of the midnasal cavity",
        "Abnormality of the nasal cavity",
        "Median cleft palate"
    ),
    SHH = c(
        "Pituitary dwarfism",
        "Ectopic anterior pituitary gland",
        "Median cleft lip and palate",
        "Median cleft palate",
        "Median cleft lip",
        "Cleft lip"
    ),
    NOTCH = c(
        "Pituitary dwarfism",
        "Ectopic anterior pituitary gland",
        "Anterior pituitary agenesis",
        "Median cleft lip and palate",
        "Basal encephalocele",
        "Ethmoidal encephalocele",
        "Abnormal hard palate morphology",
        "Midnasal stenosis",
        "Abnormality of the midnasal cavity",
        "Abnormality of the nasal cavity",
        "Median cleft palate",
        "Cyclopia"
    ),
    FGF = c(
        "Pituitary dwarfism",
        "Ectopic anterior pituitary gland",
        "Anterior pituitary agenesis",
        "Median cleft lip and palate",
        "Hypoplastic iris stroma",
        "Basal encephalocele",
        "Ethmoidal encephalocele",
        "Midnasal stenosis",
        "Abnormality of the midnasal cavity",
        "Abnormality of the nasal cavity",
        "Median cleft palate",
        "Cyclopia",
        "Median cleft lip"
    ),
    WNT_PCP = c(
        "Abnormal hard palate morphology",
        "Median cleft lip",
        "Abnormal hypothalamus physiology",
        "Expressive language delay",
        "Abnormality of the diencephalon"
    )
)
