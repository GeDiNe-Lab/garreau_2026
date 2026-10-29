library(data.table)
library(ontologyIndex) # install.packages("ontologyIndex")

setwd(paste(here::here(), "projects", sep = "/"))

out_dir <- "article/transit_data/MGI_reports"
dir.create(out_dir, showWarnings = FALSE)
base_url <- "https://www.informatics.jax.org/downloads/reports/"

reports <- c(
    pheno  = "MGI_PhenoGenoMP.rpt", # genotype -> MP term -> PMID -> marker(s)
    homol  = "HOM_MouseHumanSequence.rpt", # human-mouse orthology
    marker = "MRK_List2.rpt", # marker types (to drop transgenes)
    mp     = "MPheno_OBO.ontology" # Mammalian Phenotype ontology
)
for (f in reports) {
    dest <- file.path(out_dir, f)
    print(dest)
    if (!file.exists(dest)) download.file(paste0(base_url, f), dest, mode = "wb")
}

date_file <- file.path(out_dir, "MGI_download_date.txt")
if (!file.exists(date_file)) writeLines(format(Sys.Date()), date_file)
cat("MGI reports downloaded on:", readLines(date_file), "\n")

# --- Genes to rate after re-analysis (Oct 2026) ---------------------------
recap <- fread("article/tables/Figure_4/HPE_full_cohort_signal_seed_node_recap_refined.csv")
tiers <- fread("article/transit_data/Candidate_genes_tiers_EN.tsv")

new_candidates <- recap[condition == "seed+signal" & is.na(HPE_genes), gene] # Orphanet genes excluded
genes_to_rate <- union(setdiff(new_candidates, tiers$Gene), "GLI3") # new genes + GLI3 (no longer a control)
genes_dropped <- setdiff(tiers[!grepl("Known HPE gene", Status), Gene], new_candidates)

cat("To rate:", paste(genes_to_rate, collapse = ", "), "\n")
cat("No longer candidates:", paste(genes_dropped, collapse = ", "), "\n")

# no previous curated rating for these genes: "agreement" will read
# "MGI higher - check" for every gene with a phenotype, which can be ignored
cand <- data.table(human_symbol = genes_to_rate, Tier = NA_character_, curated_in_vivo = "None")

# -----------------------------------------------------------------------------
# Human -> mouse orthologs
# -----------------------------------------------------------------------------
hom <- fread(file.path(out_dir, reports["homol"]))
setnames(hom, make.names(names(hom)))
hs <- hom[Common.Organism.Name == "human", .(DB.Class.Key, human_symbol = Symbol)]
mm <- hom[
    Common.Organism.Name == "mouse, laboratory",
    .(DB.Class.Key, mouse_symbol = Symbol, mouse_mgi = Mouse.MGI.ID)
]
orth <- merge(hs, mm, by = "DB.Class.Key")[, DB.Class.Key := NULL]
cand <- merge(cand, orth, by = "human_symbol", all.x = TRUE)
if (any(is.na(cand$mouse_mgi))) {
    warning(
        "No mouse ortholog found for: ",
        paste(cand[is.na(mouse_mgi), human_symbol], collapse = ", ")
    )
}

# -----------------------------------------------------------------------------
# Marker types: identify transgenes (Cre drivers etc.)
# -----------------------------------------------------------------------------
mrk <- fread(file.path(out_dir, reports["marker"]), quote = "")
setnames(mrk, make.names(names(mrk)))
transgene_ids <- mrk[Marker.Type == "Transgene", MGI.Accession.ID]

# -----------------------------------------------------------------------------
# Genotype-phenotype annotations (file has no header)
# Columns: allelic composition, allele symbol(s), genetic background,
#          MP ID, PubMed ID, MGI marker ID(s) (comma-separated)
# -----------------------------------------------------------------------------
ph <- fread(file.path(out_dir, reports["pheno"]),
    header = FALSE, sep = "\t",
    quote = "", fill = TRUE
)
setnames(ph, 1:6, c(
    "allelic_composition", "alleles", "background",
    "mp_id", "pmid", "marker_ids"
))
ph <- ph[, 1:6]

# one row per (annotation, marker)
ph[, row_id := .I]
ph_long <- ph[, .(mouse_mgi = trimws(unlist(strsplit(marker_ids, ",")))), by = row_id]
# number of non-transgene markers per genotype
ph_long[, n_genes := sum(!mouse_mgi %in% transgene_ids), by = row_id]
ph <- merge(ph, unique(ph_long[, .(row_id, n_genes)]), by = "row_id")

ann <- merge(ph_long[mouse_mgi %in% cand$mouse_mgi, .(row_id, mouse_mgi)],
    ph,
    by = "row_id"
)
ann[, genotype_type := fifelse(n_genes == 1, "single-gene", "multi-gene")]
ann[, conditional := grepl("Tg\\(|cre|Cre", alleles)]

# -----------------------------------------------------------------------------
# MP ontology: names and the branches used for the rating
# -----------------------------------------------------------------------------
mp <- get_OBO(file.path(out_dir, reports["mp"]),
    propagate_relationships = "is_a", extract_tags = "minimal"
)
id_of <- function(term_name) {
    id <- names(mp$name)[tolower(mp$name) == tolower(term_name)]
    if (length(id) != 1) {
        warning("MP term not found, skipped (check the name in the MP browser): ", term_name)
        return(NULL)
    }
    id
}
desc_of <- function(term_names) {
    ids <- unlist(lapply(term_names, id_of))
    unique(unlist(lapply(ids, function(i) get_descendants(mp, i))))
}

# Moderate: structural brain, neural tube or commissural defect
moderate_branches <- desc_of(c(
    "abnormal brain morphology",
    "abnormal neural tube morphology",
    "abnormal axon guidance",
    "abnormal brain commissure morphology"
))
# Limited: any other nervous system phenotype, or embryonic/perinatal lethality
limited_branches <- desc_of(c(
    "nervous system phenotype",
    "embryonic lethality",
    "perinatal lethality"
))

# Strong: HPE, midline or ventral forebrain defect (matched on term names)
strong_regex <- paste(
    "holoprosenceph", "cyclopia", "synophthalm", "proboscis",
    "hypotelorism", "single median", "midline", "ventral forebrain",
    "ventral telencephalon", "absent forebrain", "absent telencephalon",
    "fused cerebral hemispheres", "medial ganglionic eminence",
    sep = "|"
)
moderate_regex <- paste(
    "neural tube", "exenceph", "anenceph", "spina bifida", "craniorachischisis",
    "commissur", "corpus callosum", "microceph", "hydroceph",
    "cortical layer", "neuron migration", "axon guidance",
    sep = "|"
)

ann[, mp_name := mp$name[mp_id]]
ann[, level := fcase(
    grepl(strong_regex, tolower(mp_name)), "Strong",
    mp_id %in% moderate_branches | grepl(moderate_regex, tolower(mp_name)), "Moderate",
    mp_id %in% limited_branches, "Limited",
    default = "Other"
)]
ann[, level := factor(level,
    levels = c("Strong", "Moderate", "Limited", "Other"),
    ordered = TRUE
)]

# best (highest) level reached by a gene, "None" if only "Other" terms
best_level <- function(l) {
    l <- l[l != "Other"]
    if (length(l) == 0) "None" else as.character(min(l))
}

ann <- merge(cand, ann, by = "mouse_mgi", all.x = FALSE)
fwrite(
    ann[order(human_symbol, level)][
        ,
        .(
            human_symbol, mouse_symbol, Tier, genotype_type, conditional,
            allelic_composition, background, mp_id, mp_name, level, pmid
        )
    ],
    # file.path(out_dir, "MGI_annotations_candidates_long.tsv"),  # Sept 2026 run
    file.path(out_dir, "MGI_annotations_candidates_long_update.tsv"),
    sep = "\t"
)

# -----------------------------------------------------------------------------
# Per-gene summary (single-gene genotypes drive the rating)
# -----------------------------------------------------------------------------
collapse_top <- function(x, n = 5) paste(head(unique(x), n), collapse = "; ")

summ <- ann[genotype_type == "single-gene", .(
    n_genotypes = uniqueN(allelic_composition),
    n_annotations = .N,
    mgi_in_vivo = best_level(level),
    strong_terms = collapse_top(mp_name[level == "Strong"]),
    moderate_terms = collapse_top(mp_name[level == "Moderate"]),
    limited_terms = collapse_top(mp_name[level == "Limited"]),
    pmids_top_level = collapse_top(pmid[as.character(level) == best_level(level)], 10)
), by = .(human_symbol, mouse_symbol)]

summ <- merge(cand[, .(human_symbol, mouse_symbol, Tier, curated_in_vivo)],
    summ,
    by = c("human_symbol", "mouse_symbol"), all.x = TRUE
)
summ[is.na(mgi_in_vivo), `:=`(mgi_in_vivo = "None", n_genotypes = 0L, n_annotations = 0L)]

# multi-gene genotypes only: listed separately, not used in the rating
multi <- ann[genotype_type == "multi-gene", .(
    multi_gene_terms = collapse_top(mp_name[level %in% c("Strong", "Moderate")])
), by = human_symbol]
summ <- merge(summ, multi, by = "human_symbol", all.x = TRUE)

rank_lvl <- c(Strong = 3, Moderate = 2, Limited = 1, None = 0)
summ[, curated_in_vivo := fifelse(
    curated_in_vivo == "" | is.na(curated_in_vivo),
    "None", curated_in_vivo
)]
summ[, agreement := fcase(
    rank_lvl[mgi_in_vivo] == rank_lvl[curated_in_vivo], "same",
    rank_lvl[mgi_in_vivo] > rank_lvl[curated_in_vivo], "MGI higher - check",
    default = "MGI lower - check"
)]

setorder(summ, Tier, human_symbol)
# fwrite(summ, file.path(out_dir, "MGI_summary_candidates.tsv"), sep = "\t")  # Sept 2026 run
fwrite(summ, file.path(out_dir, "MGI_summary_candidates_update.tsv"), sep = "\t")

print(summ[, .N, by = agreement])
print(summ[
    agreement != "same",
    .(
        human_symbol, Tier, curated_in_vivo, mgi_in_vivo,
        strong_terms, moderate_terms
    )
])
