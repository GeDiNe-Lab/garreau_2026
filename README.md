# A context-specific co-expression network and patient-wise network propagation to prioritize candidate genes in holoprosencephaly

Code and supplementary tables for:

> Garreau J. et al. (2026). *A context-specific co-expression network and patient-wise network propagation to prioritize candidate genes in holoprosencephaly.* [Journal, DOI to be added]

## Overview

Holoprosencephaly (HPE) remains genetically unsolved in most patients, and part of the missing heritability is thought to be oligogenic. This repository contains the analysis pipeline used to:

1. build a co-expression network from an iPSC model of anterior neuroectoderm exposed to a cyclopamine (SHH inhibitor) gradient, restricted to genes expressed in the early human forebrain;
2. select rare, potentially damaging variants in each patient's exome as seed genes;
3. run a Random Walk with Restart (RWR) from each patient's seeds, aggregate the profiles per cohort and detect signal genes against a permutation null;
4. compare the HPE cohort with a control cohort (premature ovarian insufficiency, POI);
5. rate the genes that are both signal genes and seeds on three lines of evidence (mechanistic link to HPE, mouse phenotype, loss-of-function constraint).

## Repository structure

```
.
├── scripts/                  # data processing, network propagation and function files
│   └── Figures/              # analyses producing the figures and tables of the article
├── supplementary_tables/     # S1 to S13 Tables of the article
└── README.md
```

## Scripts

### `scripts/`

| Script | Role |
|---|---|
| `scrna.R` | Selects the genes expressed in the early forebrain from the single-cell atlas of Zeng et al. 2023. |
| `tprimrna.R` | Normalizes the iPSC RNA-seq data and builds the co-expression network. |
| `patient_data_HPE_POI.R` | Filters the exome variants of the HPE and POI patients and selects the seed variants. |
| `rwr.R` | Runs one Random Walk with Restart per patient from its seed genes. |
| `rwr_normalization_fct.R` | Functions to build the permutation bank and Z-score normalize RWR profiles. |
| `signal_gene_fct.R` | Functions to aggregate patient profiles, estimate the empirical null and call signal genes. |
| `ora_fct.R` | Functions for over-representation analyses (GO-BP, HPO, Reactome) and dot plots. |
| `custom_fct.R` | Shared helper functions (RWR, gene ID conversion, plotting theme). |
| `PUBMED_query.R` | Runs the PubMed queries used to rate the mechanistic link to HPE. |
| `MGI_query.R` | Retrieves and rates the mouse phenotypes of the candidate genes from MGI. |
| `loeuf_query.R` | Retrieves the gnomAD v4.1.2 constraint metrics (LOEUF) of the candidate genes. |

### `scripts/Figures/`

| Script | Role |
|---|---|
| `ethseq_analysis.R` | Infers the ancestry of the HPE and POI patients with EthSEQ (S1 Table). |
| `network_analysis.R` | Network topology and enrichment of the 500 most connected genes (Fig 1–2, S6–S7 Tables). |
| `seeds_analysis.R` | Seed composition, seed counts, retention tests and seed enrichment (Fig 3–4, S4, S8–S10 Tables). |
| `cds_weighted_bank_generation.R` | Builds the permutation bank with random seeds drawn in proportion to CDS length. |
| `rwr_results_analysis.R` | Signal genes, Z-score distributions, rank vs degree and signal gene enrichment (Fig 5–7, S11 Table). |
| `rwr_cohort_permutations.R` | Cohort label permutation test on the number of signal genes. |
| `rwr_subnetwork_analysis.R` | Builds the candidate subnetwork and the gene recap table (Fig 8, S12 Table). |

## Supplementary tables

| Table | Content |
|---|---|
| S1 Table | Ancestry of the HPE and POI patients (EthSEQ) |
| S2 Table | RNA-seq samples used to build the co-expression network |
| S3 Table | In silico predictors and thresholds used in the pathogenicity metascore |
| S4 Table | Established HPE genes (Tekendo-Ngongang et al. 2026; Dubourg et al. 2018), network position and seed counts |
| S5 Table | PubMed queries used to rate the mechanistic link to HPE |
| S6 Table | Network gene summary (degree, betweenness, closeness, eigenvector centrality) |
| S7 Table | Enrichment of the 500 most connected network genes |
| S8 Table | Seed variant composition and seed counts per patient |
| S9 Table | Seed frequency per gene in the HPE and POI cohorts |
| S10 Table | Enrichment of the HPE and POI seed genes |
| S11 Table | Enrichment of the HPE signal genes |
| S12 Table | Candidate subnetwork (node and edge lists) |
| S13 Table | Evidence ratings and tiers of the candidate genes |

## Requirements

R with the following packages:

- **Bioconductor:** `DESeq2`, `limma`, `org.Hs.eg.db`, `AnnotationDbi`, `ComplexHeatmap`, `DEGreport`
- **Networks and propagation:** `igraph`, `WGCNA`, `RANKS`, `NetPreProc`, `EGAD`
- **Enrichment and queries:** `WebGestaltR`, `rentrez`, `ontologyIndex`
- **Data handling:** `tidyverse`, `data.table`, `Matrix`, `here`
- **Single-cell:** `Seurat`
- **Plotting:** `ggplot2`, `ggrepel`, `ggpubr`, `ggsignif`, `patchwork`, `circlize`, `paletteer`, `colorblindr`, `factoextra`, `svglite`
