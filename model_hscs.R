library(viridisLite)
library(tidyverse)
library(data.table)
library(parallel)
source("code/functions/analysis_variables.R")
library(Rcpp)

# load Watson DNMT3A tables (requires interactively selecting the table
DNMT3A_watson_file = "processed_data/watson_DNMT3A_supp.tsv"
if (!file.exists(DNMT3A_watson_file)) {
  pdf = tabulapdf::extract_areas("~/Downloads/aay9333_watson_sm.pdf", pages = c(54, 55))
  nrow(pdf[[1]])

  # DNMT3A mutation rates
  page54 = pdf[[1]][-1:-3,]
  page55 = pdf[[2]][-1:-4,] |> `colnames<-`(pdf[[2]][1,])
  page55_select = page55[]
  dim(page54)
  dim(page55)
  page55[,5] = as.numeric(gsub(" .*", "", page55[[5]]))
  colnames(page55)[5] = colnames(page54)[5]
  pages_DNMT3A = rbind(page54[, 1:5], page55[,1:5])
  pages_DNMT3A = pages_DNMT3A |>
    mutate(s = as.numeric(s),
         `Site-specific` = as.numeric(`Site-specific`),
         aachange = gsub("p.", "", DNMT3A))

  fwrite(pages_DNMT3A, DNMT3A_watson_file)
} else {pages_DNMT3A = fread(DNMT3A_watson_file)}

# get the mutation rates for all respective DNMT3A mutations
# This function is quite often used: Check if this is possible to make into function:
# laos important to check if it is possible to check the mutation rate for a specific mutation
tissue = "blood"
ncells = 1e5
expected_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz")) # load the relative mutation ratio for each sample for the DNMT3A gene:
ratios = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))
sig_donor_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_sig_donor_rates.tsv.gz"))
metadata = fread("processed_data/blood/blood_metadata.tsv")

CH_bDM = fread("processed_data/boostdm/boostdm_genie_cosmic/CH_boostDM_cancer.txt.gz")
DNMT3A_sites = CH_bDM[aachange %in% pages_DNMT3A$aachange & gene_name == "DNMT3A", .N, c("gene_name", "mut_type", "aachange", "driver")]

expected_rate_aachange = expected_rates |>
  inner_join(DNMT3A_sites, by = "mut_type", relationship = "many-to-many") |>
  left_join(metadata, by = c("sampleID", "coverage", "category")) |>
  left_join(ratios, by = c("category", "gene_name")) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells * N)) |>
  group_by(donor, aachange, category, sampleID) |>
  summarise(across(c(mle, cilow, cihigh), sum),
            age = mean(age), .groups = "drop_last") |>
  summarise(across(c(mle, cilow, cihigh, age), mean), .groups = "drop")

output = data.frame(aachange = unique(expected_rate_aachange$aachange),
                    intercept = NA, slope = NA)
rownames(output) = output$aachange

for (select_aachange in output$aachange) {
  data_df = expected_rate_aachange |> filter(aachange %in% select_aachange)
  model = lm(mle ~ age, data = data_df) |> summary()

  intercept = model$coefficients[1,1]/ ncells
  output[select_aachange, 2] = intercept
  output[select_aachange, 3] = model$coefficients[2,1]/ ncells

  # # if negative intercept, set fixed intercept at 0
  if (intercept <= 0) {
  model_fixed_icpt = lm(mle ~ age + 0, data = data_df) |> summary()
  output[select_aachange, 2] = 0
  output[select_aachange, 3] = model_fixed_icpt$coefficients[1,1]/ ncells
  }

}

# do some plotting to get more information about the full differences:
ggplot(output, aes(x = intercept, y = slope )) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = aachange))

pages_DNMT3A |>
  ggplot(aes(x = `Site-specific`, y = s)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label = DNMT3A))

mrate_watson = pages_DNMT3A |>
  left_join(output, join_by(aachange))  |>
  as.data.frame() |> `rownames<-`(pages_DNMT3A$aachange)

ggplot(mrate_watson, aes(x = slope, y = `Site-specific`*1e-9)) +
  geom_point(aes(color = s)) +
  ggrepel::geom_text_repel(data = mrate_watson |> filter(Observed > 10),
                           aes(label = aachange), max.overlaps = Inf) +
  geom_abline(slope = 1) +
  scale_color_viridis_c()

output_list = list()
for (select_aachange in mrate_watson$aachange) {
  print(select_aachange)
  output_list[[select_aachange]] = full_modeling(nreps = 1e3, B = 1e5, s =  mrate_watson[select_aachange, "s"]*0.01,
                mrate = mrate_watson[select_aachange, "slope"],
                intercept = mrate_watson[select_aachange, "intercept"], max_time = 90)
}

clonal_growth = rbindlist(output_list, idcol = "aachange")
clonal_growth |>
  arrange(fraction)

# start modeling 1e6 individuals to get a smoother line
clonal_growth |>
  ggplot(aes(x = age, y = fraction, group = aachange)) +
  geom_line() +
  ggrepel::geom_text_repel(data= clonal_growth |> filter(age == max(age)), aes(label = aachange),hjust = -2, max.overlaps = 1) +
  cowplot::theme_cowplot() +
  labs(title = "CH across individuals - measured by growth rate")

# total DNMT3A mutations
clonal_growth_total = clonal_growth |>
  group_by(age) |>
  summarize(fraction = sum(fraction)) |>
  mutate(type = "CH_modeling")

ggplot(clonal_growth_total, aes(x = age, y = fraction)) +
  geom_line() +
  cowplot::theme_cowplot() +
  labs(title = "fraction of individuals with DNMT3A mutation")

# check the overlap of BoostDM sites with the site sfor which we have and 's'
ch_driver_aa = CH_bDM |> filter(gene_name == "DNMT3A" & driver == TRUE) |> pull(aachange) |> unique()
s_DNMT3A = clonal_growth |> pull(aachange) |> unique()

euler_data = eulerr::euler(list(`BoostDM drivers` = ch_driver_aa, `Watson et al., \n mutations\n N = 102` = s_DNMT3A))
plot(euler_data, quantities = TRUE)

# compare simulated data against the ukbiobank data
UKB_DNMT3A = fread("raw_data/UKBiobank/UKB_age_frequencies_DNMT3A.tsv")
UKB_DNMT3A = UKB_DNMT3A |>
  filter(Individuals > 5000) |> # filter for the years when at least 5000 patients are in the cohort
  mutate(R882H = `R/H` / Individuals,
         R882C = `R/C` / Individuals,
         R882S = `R/S` / Individuals,
         R882P = `R/P` / Individuals) |>
  dplyr::rename(age = Age)

# for now compare the clonal growth of the main four mutations
UKB_DNMT3A_data = UKB_DNMT3A |>
  select(-contains("/"), -Individuals, -DNMT3A) |>
  pivot_longer(contains("R882"), names_to = "aachange", values_to = "fraction") |>
  mutate(type = "UKB_DNMT3A")

growth_DNMT3A_data = rbind(clonal_growth, UKB_DNMT3A_data) |>
  mutate(linegroup = paste0(aachange, type)) |>
  filter(aachange %in% c("R882H", "R882C", "R882S", "R882P"))

ggplot(growth_DNMT3A_data, aes(x = age, y = fraction, group = linegroup, color = type)) +
          geom_line() +
  scale_color_manual(values = c(blood_colors[[1]], "black")) +
  facet_wrap(aachange ~ ., scales = "free_y") +
  theme_bw()


ggplot(growth_DNMT3A_data |> filter(aachange %in% c("R882C", "R882H")), aes(x = age, y = fraction, group = linegroup, color = type)) +
  geom_line() +
  scale_color_manual(values = c(blood_colors[[1]], "black")) +
  facet_wrap(aachange ~ ., scales = "free_y") +
  theme_bw()

# growth all
growth_DNMT3A_data |>
  group_by(type, age) |>
  summarize(fraction = sum(fraction)) |>
  ggplot(aes(x = age, y = fraction, color = type)) +
  geom_line() +
  scale_color_manual(values = c(blood_colors[[1]], "black")) +
  theme_bw() + labs(title = "R882H/C/P/S mutations compared to UKbiobank incidence")


# also compare all mutations against all sites in DNMT3A

# specifically select against sites in DNMT3A:
UKB = fread("~/Downloads/UkBiobank_DNMT3A_mut_age.csv")
UKB = UKB |> dplyr::count(aa_change, Age) |>
  dplyr::rename(age = Age) |>
  left_join(UKB_DNMT3A |> select(age, Individuals)) |>
  filter(Individuals > 5000)

UKB_fraction = UKB |>
  mutate(fraction = n / Individuals)

UKB_fraction_all = UKB_fraction |>
  mutate(type = "UKB all DNMT3A") |>
  group_by(age, type) |>
  summarize(fraction = sum(fraction)) |>
  select(age, fraction, type)


rbind(clonal_growth_total, UKB_fraction_all) |>
  ggplot(aes(x = age, y = fraction, color = type)) +
  geom_line() +
  scale_color_manual(values = c(blood_colors[[1]], "black")) +
  theme_bw() + labs(title = "CH modeling compared to UKbiobank incidence")

# ideas for the bladder project: plot the mutations across the samples -
# due to what reason are they induces?
