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
load("article/transit_data/patient_variants.RData")

patient_ids <- str_extract(basename(patient_files_all), "(?<=HPE_WESA_)[^\\.]+")
patient_ids %>% write.csv("article/transit_data/patient_ids.csv", row.names = FALSE)

patient_ids

clean_clinic <- read.csv("/home/jules/Documents/phd/Data/vcf_HPE/clinic_clean.csv", header = TRUE, stringsAsFactors = FALSE)

clean_clinic_f <- filter(clean_clinic, patient %in% patient_ids, genes != "")

clean_clinic_f %>% View()

patient_to_check <- patient_variants[which(patient_ids %in% clean_clinic_f$patient)]


diag_check <- seq_len(nrow(clean_clinic_f)) %>%
    lapply(function(i) {
        idx <- match(clean_clinic_f$patient[i], patient_ids)
        diag_genes <- strsplit(clean_clinic_f$genes[i], "/")[[1]] %>% trimws()

        present <- diag_genes %in% patient_variants[[idx]]$Gene.refGeneWithVer

        data.frame(
            gene = diag_genes,
            absent = !present,
            patient_id = clean_clinic_f$patient[i],
            patient_number = idx,
            stringsAsFactors = FALSE
        )
    }) %>%
    bind_rows()

# Genes diag absents de patient_variants, par patient
diag_check %>% filter(absent)
patient_variants[[13]] %>% View()
