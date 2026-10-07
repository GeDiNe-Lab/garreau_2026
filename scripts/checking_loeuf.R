library(readr)
library(dplyr)

setwd(paste(here::here(), "projects", sep = "/"))

# --- Candidate genes from your file ------------------------------------
cand <- read_tsv("article/transit_data/Genes candidats Garreau 2026.csv",
    locale = locale(encoding = "UTF-8"),
    show_col_types = FALSE
)
colnames(cand)[1] <- "gene" # rename the first column
cand <- cand %>%
    filter(!is.na(gene), gene != "") %>%
    mutate(gene = trimws(gene))

# --- gnomAD constraint ---------------------------------------------------
constraint <- read_tsv("article/transit_data/gnomad.v4.1.2.constraint_metrics.tsv.bgz",
    show_col_types = FALSE
)

loeuf <- constraint %>%
    filter(mane_select == TRUE) %>% # one transcript per gene
    select(gene,
        LOEUF = lof.oe_ci.upper,
        pLI   = lof.pLI,
        mis_z = mis.z_score
    ) %>%
    distinct(gene, .keep_all = TRUE)

# --- Join ----------------------------------------------------------------
cand_constraint <- cand %>%
    left_join(loeuf, by = "gene") %>%
    relocate(LOEUF, pLI, mis_z, .after = gene) %>%
    arrange(LOEUF)

# Genes with no MANE match (symbol mismatch or missing in gnomAD)
setdiff(cand$gene, loeuf$gene)

write_tsv(cand_constraint, "article/transit_data/Genes_candidats_constraint.tsv")
