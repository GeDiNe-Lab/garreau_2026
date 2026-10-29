# A context-specific co-expression network and patient-wise network propagation to prioritize candidate genes in holoprosencephaly

Code and supplementary tables for:

> Garreau J. et al. (2026). *A context-specific co-expression network and patient-wise network propagation to prioritize candidate genes in holoprosencephaly.* [Journal, DOI to be added]

## Overview

Holoprosencephaly (HPE) remains genetically unsolved in most patients, and part of the missing heritability is thought to be oligogenic. This repository contains the analysis pipeline used to:

1. build a co-expression network from an iPSC model of anterior neuroectoderm exposed to a cyclopamine (SHH inhibitor) gradient (3,086 genes);
2. select rare, potentially damaging variants in each patient's exome as seed genes;
3. run a Random Walk with Restart (RWR) from each patient's seeds, aggregate the profiles per cohort and detect signal genes against a CDS length-weighted permutation null;
4. compare the HPE cohort with a control cohort (premature ovarian insufficiency, POI), including a cohort label permutation test;
5. rate the genes that are both signal genes and seeds on three lines of evidence (mechanistic link to HPE, mouse phenotype, loss-of-function constraint) and sort them into tiers.

## Repository structure

```
.
├── scripts/                  # analysis scripts (R)
├── supplementary_tables/     # S1 to S13 Tables of the article
└── README.md
```

⚠️ Adjust the tree to the actual folder names before publishing.

## Scripts

| Script | Role | Article |
|---|---|---|
| `custom_fct.R` | Shared helper functions (gene ID conversion, plotting themes, dot plots) | All |
| `RECENT_webgestalt_ORA.R` | Over-representation analyses with WebGestaltR (GO-BP, HPO, Reactome) | Fig 2, 4, 7 |
| `patient_data_alternate.R` | Seed variant filtering for the HPE cohort | Methods, Fig 3–4 |
| `iop_patient_data.R` | Seed variant filtering for the POI cohort | Methods, Fig 3–4 |
| `rwr.R` | Patient-wise Random Walk with Restart | Methods |
| `RECENT_rwr_normalization.R` | Seed score imputation, cohort aggregation and Z-score normalization | Methods |
| `RECENT_patient_gene_pvalue.R` | Empirical null, p-values and signal gene detection | Methods, Fig 5 |
| `RECENT_score_distrib.R` | Score distributions of the signal genes | Fig 5 |
| `RECENT_results_analyses.R` | Downstream analyses of the RWR results | Fig 5–7 |
| `Figure_1.R` | Network topology and enrichment of the 500 most connected genes | Fig 1–2, S6–S7 Tables |
| `Figure_2.R` | Seed composition, seed counts, seed enrichment and retention tests | Fig 3–4, S8–S10 Tables |
| `Figure_3.R` | Signal genes, rank vs degree, enrichment of the signal genes | Fig 5–7, S11 Table |
| `Figure_3_label_permutation.R` | Cohort label permutation test | Results part 3 |
| `Figure_4.R` | Candidate subnetwork | Fig 8, S12 Table |
| `PUBMED_query.R` | PubMed queries for the mechanistic link to HPE | Methods, S5 Table |
| `MGI_query.R` | Mouse phenotypes from MGI | Methods, S13 Table |

⚠️ Check that every script listed here is in the repository and that its role is right, then remove the scripts that are not used in the final version.

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

## Data availability

- **Patient exomes.** HPE and POI exome data cannot be shared publicly because of patient consent. The scripts are provided so that the pipeline can be reproduced on other cohorts.
- **iPSC RNA-seq data.** Raw data are under submission to the European Genome-phenome Archive (EGA) and can be provided by the authors upon reasonable request.
- **Single-cell RNA-seq reference.** Zeng et al. 2023, Gene Expression Omnibus [GSE155121](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155121).

## Requirements

R (≥ 4.x ⚠️ version to fill in) with the following main packages: `tidyverse`, `Matrix`, `igraph`, `WGCNA`, `WebGestaltR`, `rentrez`, `ggplot2`, `ggrepel`, `svglite`, `openxlsx`.

⚠️ Add the exact versions (for example with `sessionInfo()` or the `renv.lock` file if the project uses renv).

## Contact

Jules Garreau, GeDiNe team, Univ Rennes, CNRS, INSERM, IGDR (Institut de Génétique et Développement de Rennes), UMR 6290, ERL U1305, F-35000 Rennes, France.

## Citation

If you use this code, please cite the article above.
