# UKB incidence of CH mutations:
# Script focusing on the UKbiobank analyses
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(UpSetR)
library(mutrisk)
source("code/functions/analysis_variables.R")

# CH frequency
UKB_DNMT3A = fread("raw_data/UKBiobank/UKB_age_frequencies_DNMT3A.tsv")
UKB_DNMT3A = UKB_DNMT3A |>
  filter(Individuals > 5000) |> # filter for the years when at least 5000 patients are in the cohort
  select(Age, Individuals, "R/H", "R/C", "R/S", "R/P") |>
  mutate(R882H = `R/H` / Individuals,
         R882C = `R/C` / Individuals,
         R882S = `R/S` / Individuals,
         R882P = `R/P` / Individuals) |>
  dplyr::rename(age = Age)

# tumor incidence
n_cohort = ggplot(UKB_DNMT3A, aes(x = age)) +
  geom_col(aes(y = Individuals)) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y = "number of individuals in cohort", title = "Number of individuals in the cohort") +
  theme(plot.title = element_text(hjust = 0.5))

# tumor incidence
CH_incidence = UKB_DNMT3A |>
  pivot_longer(cols = contains("/"),
               names_to = "mutation", values_to = "n_mutated") |>
                 ggplot(aes(x = age)) +
  geom_col(aes(y = n_mutated, fill = mutation))  +
  facet_grid(. ~ mutation) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y = "CH incidence", title = "CH Incidence") +
  theme(plot.title = element_text(hjust = 0.5))

# tumor incidence
UKB_CH_incidence = UKB_DNMT3A |>
  pivot_longer(starts_with("R882"),
               values_to = "incidence") |>
  select(age, name, incidence) |>
  mutate(type = "UKB_incidence")

UK_CH_incidence_plot = UKB_CH_incidence |>
  ggplot(aes(x = age, fill = name)) +
  geom_col(aes(y = incidence)) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = label_percent()) +
  labs(y = "Risk percentage", title = "Yearly incidence rates") +
  theme(plot.title = element_text(hjust = 0.5))


# get the mutation rate data
(n_cohort + UK_CH_incidence_plot) / CH_incidence

# BoostDM-CH info on DNMT3A
ch_bDM = fread("processed_data/boostdm/boostdm_genie_cosmic/CH_boostDM_cancer.txt.gz") |>
  filter(gene_name == "DNMT3A")

# load metadata
tissue = "blood"
metadata = fread(paste0("processed_data/", tissue, "/", tissue, "_metadata.tsv")) |>   distinct() |>
  mutate(category = factor(category, levels = c("normal", "chemotherapy"))) |>
  filter(category == "normal")

# load the mutation rates
ncells = tissue_ncells_ci$mid_estimate[1]
expected_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
ratios = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))


# Get the actual values for the graphic to be correct
# Expected number of cells with double mutations:
DNMT3A_counts_boostdm = ch_bDM[aachange %in% c("R882H", "R882C", "R882S", "R882P"),
                               .N, by = c("gene_name", "aachange",  "driver", "mut_type")]


expected_rate_DNMT3A = expected_rates |>
  inner_join(DNMT3A_counts_boostdm, by = "mut_type") |>
  inner_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells)) |>
  group_by(donor, category, aachange) |>
  summarise(across(c(mle, age), mean)) |>
  setDT()

vogelgram_mut_risks_long = expected_rate_DNMT3A |>
  dplyr::rename(name = aachange, incidence = mle) |>
  select(age, name, incidence) |>
  mutate(type = "estimate")

# get the 'vogelgram' indidence rates
df_vogelgram_incidence = rbind(UKB_CH_incidence, vogelgram_mut_risks_long) |>
  mutate(aachange = case_when(grepl("S$", name) ~ "R882S",
                              grepl("H$", name) ~ "R882H",
                              grepl("C$", name) ~ "R882C",
                              grepl("P$", name) ~ "R882P"))

base_plot = ggplot(df_vogelgram_incidence) +
  geom_point(aes(x = age, y = incidence,
                 color = aachange, alpha = type)) +
  scale_alpha_manual(values = c(1, 1)) +
  cowplot::theme_cowplot()

point_plot = base_plot +
  scale_alpha_manual(values = c(0.5, 1))
point_plot

base_plot +
  facet_wrap(. ~ type, scale = "free")

point_plot +
  facet_wrap(. ~ name, scale = "free")

plot_rates_log = point_plot +
  scale_y_log10()
plot_rates_log

plot_rates_log_facet = plot_rates_log +
  facet_grid(. ~ aachange)
plot_rates_log_facet

fold_increase = df_vogelgram_incidence |>
  filter(age > 40) |>
  group_by(aachange, type) |>
  summarize(mean = mean(incidence),
            mean_age = mean(age)) |>
  select(-mean_age) |>
  pivot_wider(values_from = mean, names_from = type) |>
  mutate(fold_change = paste0("est./obs.\n", round(estimate/UKB_incidence, 2)))

plot_rates_log_facet +
  geom_text(data = fold_increase, aes(x = 40, y = 0.025,
                                      label = fold_change))


# get the yearly mutation rates


rates = expected_rate_DNMT3A |>
  mutate(mrate = mle / ncells) |>
  select(-mle) |>
  pivot_wider(names_from = aachange, values_from = mrate)

# rates
ggplot(rates, aes(x = age, y= R882H)) +
  geom_point()

df_rate = data.frame(names = c("R882H", "R882C", "R882S", "R882P"),
  yearly_mutation_rate = c(lm(R882H  ~ age, data = rates)$coef[[2]],
  lm(R882C ~ age, data = rates)$coef[[2]],
  lm(R882S ~ age, data = rates)$coef[[2]],
  lm(R882P ~ age, data = rates)$coef[[2]]))

df_rate = df_rate |>
  mutate(label = paste0(format(yearly_mutation_rate/1e-9, digits = 3), "x10^-9"))


ggplot(df_rate, aes(x = names, y = yearly_mutation_rate, fill = names)) +
  geom_col() +
  cowplot::theme_cowplot() +
  geom_text(aes(label = format(yearly_mutation_rate, digits = 3), y = yearly_mutation_rate),
            vjust = -0.2, position = position_dodge(0.9))


