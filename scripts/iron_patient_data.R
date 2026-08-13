# Loading packages and functions
library(Matrix)
library(tidyverse)
library(parallel)
library(tibble)

# Setting working directory
rstudioapi::getSourceEditorContext()$path %>%
    str_split("/") %>%
    unlist() %>%
    head(-3) %>%
    str_c(collapse = "/") %>%
    str_c("/") %>%
    setwd()

# Loading custom functions
source("article/scripts/custom_fct.R")
load("article/transit_data/data_3p.RData")

patient_files_all <- list.files(path = "/home/jules/Documents/phd/Data/vcf_FER/RESULTS/", full.names = TRUE)
patient_files_all

patient_variants <- patient_files_all %>% mclapply(function(patient_file) {
    patient_data <- read.table(patient_file, header = TRUE, sep = "\t", fill = TRUE, quote = "")
    # separate variants falling in +- strand genes
    patient_data$seq_depth <- patient_data[, 17] %>% as.numeric()
    patient_data$pLI <- cleaning_pLI(patient_data$pLI)
    patient_data <- patient_data %>%
        mutate(Gene.refGeneWithVer = strsplit(as.character(Gene.refGeneWithVer), ";")) %>%
        unnest(Gene.refGeneWithVer) %>%
        as.data.frame()
    patient_data_filter <- filter(
        patient_data,
        DB_freq <= 0.001,
        seq_depth > 20,
        # Func.refGene %in% c(
        #     "exonic",
        #     "splicing",
        #     "exonic;splicing",
        #     "ncRNA_exonic;splicing",
        #     "ncRNA_splicing",
        #     "ncRNA_exonic"
        # ),
        (PRED_value >= 0.8 | CLNSIG %in% c("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic")),
        !PRED_effect %in% c("nonframeshift substitution", "nonframeshift insertion", "nonframeshift deletion"),
        INHOUSE_FREQ <= 0.0001,
        pLI > 0.1,
        nchar(Ref) < 30,
        nchar(Alt) < 30
    )

    return(dplyr::select(
        patient_data_filter,
        Chr,
        Start,
        End,
        Ref,
        Alt,
        Gene.refGeneWithVer,
        Func.refGene,
        seq_depth,
        DB_freq,
        CLNSIG,
        pLI,
        PRED_value,
        PRED_effect,
        INHOUSE_FREQ,
        INHOUSE_COUNT,
        INHOUSE_EXP,
    ))
    gc()
}, mc.cores = 10)

save(patient_variants, file = "article/transit_data/iop_patient_variants.RData")

# Compute number of unique genes per patient
unique_gene_counts <- patient_variants %>%
    sapply(function(patient) {
        patient$Gene.refGeneWithVer %>%
            unique() %>%
            length()
    })
unique_gene_counts
# here
net_genes <- gene_converter(rownames(data_3p$coex_case), from = "ENSEMBL", to = "SYMBOL")
net_genes
unique_gene_counts_net <- patient_variants %>%
    sapply(function(patient) {
        filter(
            patient,
            Gene.refGeneWithVer %in% net_genes
        )$Gene.refGeneWithVer %>%
            unique() %>%
            length()
    })


# Convert to data frame for ggplot
df_genes <- data.frame(
    patient = as.character(c(1:21)),
    n_genes = unique_gene_counts,
    n_genes_net = unique_gene_counts_net
)

p <- ggplot(df_genes, aes(x = n_genes)) +
    geom_histogram(binwidth = 2, fill = "steelblue", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(
        x = "Number of seed per patient",
        y = "Count of patients",
        title = "Distribution of RWR seed genes across patients"
    )

print(p)

seeds <- patient_variants %>%
    sapply(function(patient) {
        return(patient$Gene.refGeneWithVer)
    })
names(seeds) <- paste0("Patient_", seq_along(seeds))
seeds[[1]]
# Match seeds with clinic data and check gene presence
clinic$seed_match <- sapply(seq_len(nrow(clinic)), function(i) {
    patient_genes <- strsplit(clinic$genes[i], "/")[[1]] %>% trimws()
    patient_seeds <- seeds[[clinic$index[i]]]

    matches <- patient_genes %in% patient_seeds
    match_count <- sum(matches)
    total_count <- length(patient_genes)

    if (match_count == 0) {
        return(NA)
    } else if (match_count == total_count) {
        return("full")
    } else {
        return("partial")
    }
})



clinic %>%
    dplyr::select(N..ADN, NGS.HPE., genes, seed_match) %>%
    View()









unique_genes <- patient_files_all %>%
    mclapply(function(patient_file) {
        patient_data <- read.table(patient_file, header = TRUE, sep = "\t", fill = TRUE, quote = "")
        patient_data <- patient_data %>%
            mutate(Gene.refGeneWithVer = strsplit(as.character(Gene.refGeneWithVer), ";")) %>%
            unnest(Gene.refGeneWithVer) %>%
            as.data.frame()
        genes <- unique(patient_data$Gene.refGeneWithVer)
        return(genes)
    }, mc.cores = 10) %>%
    unlist() %>%
    unique()
unique_genes


n <- 1000 # desired length
set.seed(42) # for reproducibility
y <- sample(unique_gene_counts, size = n, replace = TRUE)
random_cohort <- y %>%
    lapply(function(lenght) {
        sample(unique_genes, size = lenght, replace = FALSE)
    })

random_cohort_net_counts <- sapply(random_cohort, function(genes) {
    intersect(genes, net_genes) %>%
        length()
})
df_random <- data.frame(
    patient = as.character(c(1:1000)),
    n_random = y,
    n_random_net = random_cohort_net_counts
)

# Build patient dataframe: sample, gene counts before (n_genes) and after network filter (n_genes_net)
patient_samples <- basename(patient_files_all)
patient_df <- data.frame(
    sample = patient_samples,
    n_before = unique_gene_counts,
    n_after = unique_gene_counts_net,
    stringsAsFactors = FALSE
)

# Build random dataframe: sample, gene counts before (n_random = y) and after network filter (n_random_net)
random_df <- data.frame(
    sample = paste0("random_", seq_along(y)),
    n_before = y,
    n_after = random_cohort_net_counts,
    stringsAsFactors = FALSE
)

# (optional) inspect
head(patient_df)
head(random_df)
combined_df <- bind_rows(
    patient_df %>% mutate(type = "patient"),
    random_df %>% mutate(type = "random")
) %>%
    mutate(n_before = as.numeric(n_before), n_after = as.numeric(n_after))

p_comparison <- ggplot(combined_df, aes(x = n_before, y = n_after, color = type)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_smooth(aes(fill = type), method = "lm", se = FALSE, linetype = "solid", alpha = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c(patient = "steelblue", random = "orange")) +
    scale_fill_manual(values = c(patient = "steelblue", random = "orange")) +
    theme_minimal() +
    labs(
        x = "Number of genes before filter",
        y = "Number of genes after network filter",
        color = "Cohort",
        fill = "Cohort",
        title = "Patient vs Random: genes before vs after network filter"
    ) +
    coord_equal()

print(p_comparison)

random_cohort
save(random_cohort, file = "article/transit_data/random_seed_sets.RData")
