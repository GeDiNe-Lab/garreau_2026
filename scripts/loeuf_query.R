library(readr)
library(dplyr)

setwd(paste(here::here(), "projects", sep = "/"))
library(data.table)
recap <- fread("article/tables/Figure_4/HPE_full_cohort_signal_seed_node_recap_refined.csv")
tiers <- fread("article/transit_data/Candidate_genes_tiers_EN.tsv")

new_candidates <- recap[condition == "seed+signal" & is.na(HPE_genes), gene] # Orphanet genes excluded
genes_to_rate <- union(setdiff(new_candidates, tiers$Gene), "GLI3") # new genes + GLI3 (no longer a control)
genes_dropped <- setdiff(tiers[!grepl("Known HPE gene", Status), Gene], new_candidates)

cat("To rate:", paste(genes_to_rate, collapse = ", "), "\n")
cat("No longer candidates:", paste(genes_dropped, collapse = ", "), "\n")

# --- Candidate genes from your file ------------------------------------
# cand <- read_tsv("article/transit_data/Genes candidats Garreau 2026.csv",
#     locale = locale(encoding = "UTF-8"),
#     show_col_types = FALSE
# )
# colnames(cand)[1] <- "gene" # rename the first column
# cand <- cand %>%
#     filter(!is.na(gene), gene != "") %>%
#     mutate(gene = trimws(gene))
cand <- tibble(gene = genes_to_rate)

# --- gnomAD constraint ---------------------------------------------------
constraint <- read_tsv("article/transit_data/gnomad.v4.1.2.constraint_metrics.tsv.bgz",
    show_col_types = FALSE
)

loeuf <- constraint %>%
    dplyr::filter(mane_select == TRUE) %>% # one transcript per gene
    dplyr::select(gene,
        LOEUF = lof.oe_ci.upper,
        pLI   = lof.pLI,
        mis_z = mis.z_score
    ) %>%
    dplyr::distinct(gene, .keep_all = TRUE)
    
# --- Join ----------------------------------------------------------------
cand_constraint <- cand %>%
    left_join(loeuf, by = "gene") %>%
    relocate(LOEUF, pLI, mis_z, .after = gene) %>%
    arrange(LOEUF)

# Genes with no MANE match (symbol mismatch or missing in gnomAD)
setdiff(cand$gene, loeuf$gene)

# write_tsv(cand_constraint, "article/transit_data/Genes_candidats_constraint.tsv")
write_tsv(cand_constraint, "article/transit_data/Genes_candidats_constraint_update.tsv")
