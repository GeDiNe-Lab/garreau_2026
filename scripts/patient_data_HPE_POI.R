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

# Function to score variant effect based on PRED_value, PRED_effect, CLNSIG, and PRED_splice
score_variant_effect <- function(PRED_value, PRED_effect, CLNSIG, PRED_splice) {
    high_impact <- c("frameshift deletion", "frameshift insertion", "stopgain")
    moderate_high <- c("stoploss", "startloss")
    moderate <- c("nonframeshift deletion", "nonframeshift insertion", "nonframeshift substitution")

    effect_score <- dplyr::case_when(
        PRED_value == 1 & PRED_effect %in% high_impact ~ 0.9,
        PRED_value == 1 & PRED_effect %in% moderate_high ~ 0.8,
        PRED_value == 1 & PRED_effect %in% moderate ~ 0,
        PRED_value < 1 ~ pmin(pmax(PRED_value, 0), 1), # continuous missense meta-score
        TRUE ~ 0.3 # PRED_value == 1 with an unmapped/unknown categorical effect
    )

    splice_no <- c("0/1", "0/2", "0/3", "0/4")
    splice_low <- c("1/2", "1/3", "1/4", "2/4")
    splice_moderate <- c("2/2", "2/3", "3/4")
    splice_high <- c("3/3", "4/4")

    effect_score <- pmax(effect_score, dplyr::case_when(
        PRED_effect %in% splice_high ~ 0.9,
        PRED_effect %in% splice_moderate ~ 0.8,
        PRED_effect %in% splice_low ~ 0.5,
        PRED_effect %in% splice_no ~ 0,
        TRUE ~ 0.3 # PRED_value == 1 with an unmapped/unknown categorical effect
    ))

    is_clinvar_patho <- CLNSIG %in% c("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic")
    effect_score <- pmax(effect_score, ifelse(is_clinvar_patho, 0.95, 0))

    effect_score
}

# Processing HPE patient data
HPE_patient_files_all <- list.files(path = "/home/jules/Documents/phd/Data/vcf_HPE/RESULTS/", full.names = TRUE)
# extract substring between "HPE_WESA_" and ".annovar"
patient_ids <- str_extract(basename(HPE_patient_files_all), "(?<=HPE_WESA_)[^\\.]+")
patient_ids %>% write.csv("article/transit_data/patient_ids.csv", row.names = FALSE)

# Filter patient data based on quality criteria and score variant effect
HPE_patient_variants <- HPE_patient_files_all %>% mclapply(function(patient_file) {
    patient_data <- read.table(patient_file, header = TRUE, sep = "\t", fill = TRUE, quote = "")
    # separate variants falling in +- strand genes
    patient_data$seq_depth <- patient_data[, 17] %>% as.numeric()
    patient_data$pLI <- cleaning_pLI(patient_data$pLI)
    patient_data <- patient_data %>%
        mutate(Gene.refGeneWithVer = strsplit(as.character(Gene.refGeneWithVer), ";")) %>%
        unnest(Gene.refGeneWithVer) %>%
        as.data.frame()

    patient_data_qual <- filter(
        patient_data,
        DB_freq <= 0.001,
        seq_depth > 20,
        pLI >= 0.1,
        INHOUSE_FREQ <= 0.0001,
        nchar(Ref) < 30,
        nchar(Alt) < 30
    )

    patient_data_qual$PRED_value_reworked <- c(1:nrow(patient_data_qual)) %>% sapply(function(i) {
        score_variant_effect(
            patient_data_qual[i, "PRED_value"],
            patient_data_qual[i, "PRED_effect"],
            patient_data_qual[i, "CLNSIG"],
            patient_data_qual[i, "PRED_splice"]
        )
    })

    patient_data_final <- filter(
        patient_data_qual,
        PRED_value_reworked >= 0.7
    )

    return(dplyr::select(
        patient_data_final,
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
        PRED_splice,
        Splice_inter,
        PRED_value_reworked,
        INHOUSE_FREQ,
        INHOUSE_COUNT,
        INHOUSE_EXP,
    ))
    gc()
}, mc.cores = 15)

save(HPE_patient_variants, file = "article/transit_data/hpe_patient_variants.RData")

# POI patient data processing
POI_patient_files_all <- list.files(path = "/home/jules/Documents/phd/Data/vcf_IOP/RESULTS/", full.names = TRUE)
col <- colnames(read.table(POI_patient_files_all[1], header = TRUE, sep = "\t", fill = TRUE, quote = ""))

POI_patient_variants <- POI_patient_files_all %>% mclapply(function(patient_file) {
    patient_data <- read.table(patient_file, header = TRUE, sep = "\t", fill = TRUE, quote = "")
    # separate variants falling in +- strand genes
    patient_data$seq_depth <- patient_data[, 17] %>% as.numeric()
    patient_data$pLI <- cleaning_pLI(patient_data$pLI)
    patient_data <- patient_data %>%
        mutate(Gene.refGeneWithVer = strsplit(as.character(Gene.refGeneWithVer), ";")) %>%
        unnest(Gene.refGeneWithVer) %>%
        as.data.frame()

    patient_data_qual <- filter(
        patient_data,
        DB_freq <= 0.001,
        seq_depth > 20,
        pLI >= 0.1,
        INHOUSE_FREQ <= 0.0001,
        nchar(Ref) < 30,
        nchar(Alt) < 30
    )

    patient_data_qual$PRED_value_reworked <- c(1:nrow(patient_data_qual)) %>% sapply(function(i) {
        score_variant_effect(
            patient_data_qual[i, "PRED_value"],
            patient_data_qual[i, "PRED_effect"],
            patient_data_qual[i, "CLNSIG"],
            patient_data_qual[i, "PRED_splice"]
        )
    })

    patient_data_final <- filter(
        patient_data_qual,
        PRED_value_reworked >= 0.7
    )

    return(dplyr::select(
        patient_data_final,
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
        PRED_splice,
        Splice_inter,
        PRED_value_reworked,
        INHOUSE_FREQ,
        INHOUSE_COUNT,
        INHOUSE_EXP,
    ))
    gc()
}, mc.cores = 15)

save(POI_patient_variants, file = "article/transit_data/poi_patient_variants.RData")
