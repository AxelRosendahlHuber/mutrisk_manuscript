# New script making mirror plots
library(data.table)
library(tidyverse)
library(cowplot)
library(MutationalPatterns)
source("code/functions/analysis_variables.R")
getwd()

# load data sources
GENIE_data = fread("processed_data/GENIE_17/GENIE_17_genie_tissue_type.txt.gz")
metadata_files = c("processed_data/blood/blood_metadata.tsv", "processed_data/colon/colon_metadata.tsv",
                   "processed_data/lung/lung_metadata.tsv")

names(metadata_files) = str_split_i(metadata_files, "\\/", 2)
metadata = lapply(metadata_files, \(x) fread(x)[,c("sampleID", "category", "age", "donor")]) |>
  rbindlist(idcol = "tissue")

# Load gene_of_interest boostdm
boostdm = fread("processed_data/boostdm/boostdm_genie_cosmic/pancancer_boostDM_intersect.txt.gz")
boostdm_TP53 = boostdm |> filter(gene_name == "TP53")


# load the mutation rates
expected_rate_list = list()
ratio_list = list()
for (tissue in c("colon", "blood", "lung")) {
  expected_rate_list[[tissue]] = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
  ratio_list[[tissue]] = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))
}
expected_rates = rbindlist(expected_rate_list, idcol = "tissue", use.names = TRUE)
ratios = rbindlist(ratio_list, idcol = "tissue", use.names = TRUE)

# filters
ratios = ratios |> filter(gene_name == 'TP53')
expected_rates = expected_rates |> filter(category %in% c("normal", "non-smoker"))
mutated_rates = left_join(expected_rates, ratios) |>
  left_join(metadata) |>
  group_by(donor, mut_type, tissue, ratio) |>
  summarize(across(c(mle, cilow, cihigh), mean)) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio))

boostdM_driver_counts = boostdm_TP53[driver == TRUE , .N, by = c("mut_type")]


# make mirror plots:

mutated_rates |> inner_join(boostdM_driver_counts)


mutated_rates$sampleID |> n_distinct()