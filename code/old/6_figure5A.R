library(data.table)
library(tidyverse)
library(cowplot)
source("code/0_functions/analysis_variables.R")

ncells = 1e5

######
## MutRisk estimates for neutral mutation accumulation in healthy tissues: Blood
######
# load metadata files for the different samples
metadata_files = c("processed_data/blood/blood_metadata.tsv", "processed_data/colon/colon_metadata.tsv",
                   "processed_data/lung/lung_metadata.tsv")

names(metadata_files) = str_split_i(metadata_files, "\\/", 2)
metadata = lapply(metadata_files, \(x) fread(x)[,c("sampleID", "category", "age", "donor")]) |>
  rbindlist(idcol = "tissue")

# load the mutation rates
expected_rate_list = list()
ratio_list = list()
for (tissue in c("colon", "blood", "lung")) {
  expected_rate_list[[tissue]] = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
  ratio_list[[tissue]] = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))
}
expected_rates = rbindlist(expected_rate_list, idcol = "tissue", use.names = TRUE)
ratios = rbindlist(ratio_list, idcol = "tissue", use.names = TRUE)

# mutation rates
mutation_rates = expected_rates |>
  mutate(tissue_category = paste0(tissue, "_", category)) |>
  left_join(metadata)

# Mean mutation rate for each trinucleotide
# Check if the blood rates actually make sense - this seems too low - this must be because of the cord blood donors
mean_rates = mutation_rates |>
  group_by(tissue_category, mut_type) |>
  summarize(mle = mean(mle),
            mean_age = mean(age)) |>
  mutate(tissue_category_age = paste0(tissue, "\n(", format(mean_age, digits = 3, nsmall = 1 ), ") years"))

# make specific plots for specific positions:
# Blood make specific mirror genes for specific sites:
boostdm_ch = fread("processed_data/boostdm/boostdm_genie_cosmic/CH_boostDM_cancer.txt.gz")

# UKBiobank DNMT3A mutations:
UKB_DNMT3A_muts = fread("raw_data/UKBiobank/UkBiobank_DNMT3A_mut_age.csv")
UKB_DNMT3A_counts = UKB_DNMT3A_muts[, .N, by = c("aa_change", "REF", "ALT")]  |>
  mutate(position = parse_number(aa_change),
         type = paste0(REF, ">", ALT),
         type = case_match(type, .default = type,
                           "A>T" ~ "T>A", "A>G" ~ "T>C", "A>C" ~ "T>G",
                           "G>T" ~ "C>A", "G>A" ~ "C>T", "G>C" ~ "C>G")) |>
  mutate(mrate = N, tissue_category = "UKBiobank CH") |>
  select(position, type, tissue_category, mrate)

# plot not used in the main script (can be removed)
# # first attempt ukbiobank plots:
# DNMT3A_plot = UKB_DNMT3A_counts |>
#   ggplot(aes(x = position, y = N, fill = type)) +
#   geom_col() +
#   scale_fill_manual(values = COLORS6) +
#   theme_cowplot() +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
#   labs(x  = NULL , y = "Number of CH mutations\nobserved in the UKBiobank cohort", title = "DNMT3A", fill = NULL)


genes_ch = "DNMT3A"
mean_rates_blood = mean_rates |> filter(tissue_category == "blood_normal")

# only plot the driver genes
DNMT3A = boostdm_ch |> filter(gene_name == "DNMT3A")

mutations_blood_DNMT3A = left_join(DNMT3A, mean_rates_blood, relationship = "many-to-many", by = "mut_type") |>
  left_join(triplet_match_substmodel, by = "mut_type") |>
  group_by(position, tissue_category, type) |>
  summarize(mrate = sum(mle) * ncells) |>
  select(position, type, tissue_category, mrate)

library(ggh4x)
df_mirror = bind_rows(mutations_blood_DNMT3A, UKB_DNMT3A_counts) |>
  mutate(
    tissue_category = ifelse(tissue_category == "blood_normal", "Expected mutrate\nblood", tissue_category),
    tissue_category = factor(tissue_category, levels = c("UKBiobank CH", "Expected mutrate\nblood")),
    mrate = ifelse(tissue_category == "Expected mutrate\nblood", 0-mrate, mrate)) |>
  ungroup()

# way to make the plot extend both upper and lower axes
df_point = df_mirror |>
  group_by(tissue_category, position) |>
  summarize(mrate = sum(mrate), .groups = "drop_last") |>
  summarize(mrate = max(abs(mrate)) * 1.1) |>
  mutate(position = 500,
         mrate = ifelse(tissue_category == "Expected mutrate\nblood", 0-mrate, mrate)) |>
  ungroup()

F5A = ggplot(df_point, aes(x = position, y = mrate)) +
  geom_point(color = "white") +
  geom_col(data = df_mirror, aes(fill = type)) +
  geom_text(data = data.frame(tissue_category = factor("UKBiobank CH"), position = 50, mrate = 1500, label = "DNMT3A"),
            aes(label = label)) +
  facet_grid2(tissue_category ~ . , scales = "free") +
  scale_fill_manual(values = COLORS6) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = "none", panel.spacing.y = unit(0, "mm")) +
  labs(y = "Number expected/\nobserved muts",  x = "AA position") +
  scale_y_continuous(expand=expansion(mult=c(0,0)), breaks = scales::breaks_extended(n = 3), labels = abs)
F5A

saveRDS(F5A, "manuscript/figure_panels/figure_5/figure_5A.rds")
