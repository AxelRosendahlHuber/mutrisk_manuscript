# Script focusing on the UKbiobank analyses
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(UpSetR)
library(mutrisk)
source("code/functions/analysis_variables.R")

## UKBiobank analyses
# load the bowel cancer data:
crc_freq = fread("raw_data/UKBiobank/colorectal_cancer_frequency_UKB.csv")

# calculate the incidence rates
ukbiobank_crc = data.frame(age = 0:max(crc_freq$current_age)) |>
  mutate(n_alive = sapply(age, \(x) sum(crc_freq$current_age >= x)),
         n_tumor = sapply(age, \(x) sum(crc_freq$var_Colorectal_age == x, na.rm = TRUE)),
         n_no_tumor = n_alive - n_tumor,
         risk = n_tumor / n_alive,
         CRC_cumulative_risk = cumsum(risk)) |>
  filter(n_alive > 5000) # filter for the years when at least 5000 patients are in the cohort

# tumor incidence
n_cohort = ggplot(ukbiobank_crc, aes(x = age)) +
  geom_col(aes(y = n_alive)) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y = "number of individuals in cohort", title = "Number of individuals in the cohort") +
  theme(plot.title = element_text(hjust = 0.5))

# tumor incidence
CRC_incidence = ggplot(ukbiobank_crc, aes(x = age)) +
  geom_col(aes(y = n_tumor))  +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y = "CRC incidence", title = "CRC Incidence") +
  theme(plot.title = element_text(hjust = 0.5))

# tumor incidence
yearly_incidence_rates = ggplot(ukbiobank_crc, aes(x = age)) +
  geom_col(aes(y = risk)) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = label_percent()) +
  labs(y = "Risk percentage", title = "Yearly incidence rates") +
  theme(plot.title = element_text(hjust = 0.5))
# Cumulative incidence
cumulative_incidence = ggplot(ukbiobank_crc, aes(x = age)) +
  geom_col(aes(y = CRC_cumulative_risk)) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = label_percent()) +
  labs(y = "Risk percentage", title = "Cumulative Incidence") +
  theme(plot.title = element_text(hjust = 0.5))

# combine plots in an overview plot:
n_cohort + CRC_incidence + yearly_incidence_rates + cumulative_incidence

# load metadata
tissue = "colon"
metadata = fread(paste0("processed_data/", tissue, "/", tissue, "_metadata.tsv")) |>   distinct() |>
  mutate(category = factor(category, levels = c("normal", "IBD", "POLD1", "POLE")))

# load the mutation rates
ncells = tissue_ncells_ci$mid_estimate[2]
expected_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
ratios = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))

# load the mutation data:
cancer_bDM = fread("processed_data/boostdm/boostdm_genie_cosmic/colon_boostDM_cancer.txt.gz")

# Get the actual values for the graphic to be correct
APC_1450_hotspot = cancer_bDM[gene_name == "APC" & aachange == "R1450*" , c("gene_name", "mut_type", "aachange", "position", "driver")]
expected_rate_APC_1450 = expected_rates |>
  inner_join(APC_1450_hotspot, by = "mut_type") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells)) |>
  group_by(donor, category) |>
  summarise(across(c(mle, cilow, cihigh, age), mean)) |>
  setDT()

KRAS_driver_muts = cancer_bDM[gene_name == "KRAS" & driver == TRUE & position == 12,
                              c("gene_name", "mut_type", "aachange", "position")]
KRAS_G12V_hotspot = KRAS_driver_muts

expected_rate_KRAS_driver_sc = expected_rates |>
  inner_join(KRAS_G12V_hotspot, by = "mut_type") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio)) |>
  group_by(donor, category) |>
  summarise(across(c(mle, cilow, cihigh, age), mean)) |>
  arrange(category, donor, age)

# Expected number of cells with double mutations:
apc_counts_boostdm = cancer_bDM[gene_name == "APC", .N, by = c("gene_name", "mut_type",  "driver")]

apc_double_driver = expected_rates |>
  left_join(ratios |> filter(gene_name == "APC")) |>
  left_join(metadata) |>
  inner_join(apc_counts_boostdm |> filter(driver), by = "mut_type", relationship = "many-to-many") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio)) |>
  group_by(category, donor, age,  mut_type) |>
  summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop_last") |>
  summarize(across(c(mle, cilow, cihigh), sum)) |>
  mutate(across(c(mle, cilow, cihigh), ~ ((.^2) / 2) * ncells ))


# APC rates for a single cell
apc_double_driver_sc = apc_double_driver |>
  mutate(across(c(mle, cilow, cihigh), ~ ./ ncells)) |>
  arrange(category, donor, age)

# make a dataframe with the relative mutation rates for all the different changes:
vogelgram_mut_risks_sc = tibble(
  category = apc_double_driver_sc$category,
  donor = apc_double_driver_sc$donor,
  age = apc_double_driver_sc$age,
  `APC 2x SNV` = apc_double_driver_sc$mle,
  KRAS_driver = expected_rate_KRAS_driver_sc$mle) |>
  mutate(`APC 2x +KRAS` = `APC 2x SNV` * KRAS_driver)

vogelgram_mut_risks = vogelgram_mut_risks_sc |>
  mutate(across(-c(category, donor, age), ~ . * ncells))

vogelgram_mut_risks_long = vogelgram_mut_risks |>
  filter(category == "normal") |>
  pivot_longer(-c(category, donor, age), names_to = "mutation", values_to = "expected_mutated_cells")

ukbiobank_crc_plot = ukbiobank_crc |>
  filter(n_tumor != 0) |>
  select(age, CRC_cumulative_risk) |>
  pivot_longer(CRC_cumulative_risk, names_to = "name", values_to = "incidence")

vogelgram_risks_long = vogelgram_mut_risks_long |>
  dplyr::rename(name = mutation, incidence = expected_mutated_cells) |>
  select(age, name, incidence)

# plot the mutation risk for individual cells:
vogelgram_risks_long |>
  mutate(incidence) |>
  ggplot(aes(x = incidence, fill = name)) + geom_histogram(binwidth = 0.2)

# get the 'vogelgram' indidence rates
df_vogelgram_incidence = rbind(ukbiobank_crc_plot, vogelgram_risks_long) |>
  mutate(incidence_100k = incidence * 1e5,
         name = factor(name, levels = c("KRAS_driver", "APC 2x SNV", "APC 2x +KRAS", "CRC_cumulative_risk")))

vogel_plot = ggplot(df_vogelgram_incidence, aes(x = age)) +
  geom_point(aes(x = age, y = incidence_100k, color = name)) +
  annotate(geom = "text", x = 10, y = 12,
           label = "Cumulative Risk\n100,000 individuals\nUK Biobank",vjust = 0, hjust = 0, size = 3.5) +
  scale_y_log10(guide = "axis_logticks",
                breaks = c(1e-4, 1e-2,  1, 100, 1e4, 1e6),
                labels = ~ ifelse(.x < 1, scales::comma(.x), scales::comma(.x, accuracy = 1))) +
  scale_color_manual(values = c("orange", "red", "darkred", "black")) +
  theme_cowplot()  +
  theme(legend.position = "inside", legend.position.inside = c(0.6, 0.1), legend.background = element_rect(fill = alpha("white", 0.7), color = "black"),
        legend.margin = margin(1,1,1,1,"mm")) +
  coord_cartesian(clip = 'off') +
  labs(y = "Incidence per 100,000 individuals/\nMutated cells per 100,000 individuals", x = "Age (years)", color = NULL)
vogel_plot

# changes to plot to make:

# 1. Every idividual only counted once
  # this makes it possible to go for the 'fraction' of mutated reads
# 2. Reduce the cumulative risk by the risk for having those specific mutations
# 3. No log scale


# to get the probability of mutation, divide the incidence
get_prob_mutated_N = function(risk, ncells, N = 1) {
  1 - pbinom(N - 1, size = ncells, prob = risk)
}



df_vogelgram_incidence = df_vogelgram_incidence |>
  mutate(fraction_1mut = case_when(name != "CRC_cumulative_risk" ~ get_prob_mutated(incidence/ncells, ncells), .default = incidence)) |>
  mutate(dif = incidence - fraction_1mut)
df_vogelgram_incidence$incidence - df_vogelgram_incidence$fraction_1mut


vogel_plot = ggplot(df_vogelgram_incidence, aes(x = age)) +
  geom_point(aes(x = age, y = fraction_1mut, color = name)) +
  # scale_y_log10(guide = "axis_logticks",
  #               breaks = c(1e-4, 1e-2,  1, 100, 1e4, 1e6),
  #               labels = ~ ifelse(.x < 1, scales::comma(.x), scales::comma(.x, accuracy = 1))) +
  scale_color_manual(values = c("orange", "red", "darkred", "black")) +
  theme_cowplot()  +
   coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0, 1.1)) +
  labs(y = "Fraction of individuals", x = "Age (years)", color = NULL)
vogel_plot

vogel_plot = ggplot(df_vogelgram_incidence, aes(x = age)) +
  geom_point(aes(x = age, y = incidence2, color = name)) +
  scale_y_log10(guide = "axis_logticks",
                breaks = c(1e-4, 1e-2,  1, 100, 1e4, 1e6),
                labels = ~ ifelse(.x < 1, scales::comma(.x), scales::comma(.x, accuracy = 1))) +
  scale_color_manual(values = c("orange", "red", "darkred", "black")) +
  theme_cowplot()  +
  coord_cartesian(clip = 'off') +
  labs(y = "Fraction of individuals", x = "Age (years)", color = NULL)
vogel_plot


# incidence rates of more than 1 mut
df_vogelgram_incidence |>
  ggplot(aes(x = age, color = name)) +
  geom_point(aes(y = fraction_1mut), alpha = 0.3) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c("orange", "red", "darkred", "black"))


df_vogelgram_multiple_muts = df_vogelgram_incidence |>
  mutate(fraction_2mut = get_prob_mutated_N(incidence/ncells, ncells, 2),
         fraction_3mut = get_prob_mutated_N(incidence/ncells, ncells, 3),
         fraction_5mut = get_prob_mutated_N(incidence/ncells, ncells, 5),
         fraction_10mut = get_prob_mutated_N(incidence/ncells, ncells, 10))

# incidence rates of more than 1 mut
vogelgram_1_mut = df_vogelgram_multiple_muts |>
  filter(name != "CRC_cumulative_risk")  |>
  ggplot(aes(x = age, color = name)) +
  geom_point(aes(y = fraction_1mut)) +
  cowplot::theme_cowplot() +
  labs(y = "fraction of mutated individuals", title = "fraction individuals with cell with specific abberation")
vogelgram_1_mut

vogelgram_2_mut = vogelgram_1_mut +
    geom_point(aes(y = fraction_2mut)) +
    geom_segment(aes(x = age, y = fraction_1mut, yend = fraction_2mut)) +
    labs(title = "Difference fraction mutated with 1, or with 2 mutations")
vogelgram_2_mut


vogelgram_3_mut = vogelgram_1_mut +
  geom_point(aes(y = fraction_3mut)) +
  geom_segment(aes(x = age, y = fraction_1mut, yend = fraction_3mut)) +
  labs(title = "Difference fraction mutated with 1, or with 3 mutations")
vogelgram_3_mut


vogelgram_5_mut = vogelgram_1_mut +
  geom_point(aes(y = fraction_5mut)) +
  geom_segment(aes(x = age, y = fraction_1mut, yend = fraction_5mut)) +
  labs(title = "Difference fraction mutated with 1, or with 5 mutations")
vogelgram_5_mut

vogelgram_10_mut = vogelgram_1_mut +
  geom_point(aes(y = fraction_10mut)) +
  geom_segment(aes(x = age, y = fraction_1mut, yend = fraction_10mut)) +
  labs(title = "Difference fraction mutated with 1, or with 10 mutations")
vogelgram_10_mut

plots = vogelgram_1_mut + vogelgram_2_mut + vogelgram_3_mut + vogelgram_5_mut
plots_log = vogelgram_1_mut + vogelgram_2_mut + vogelgram_3_mut + vogelgram_5_mut & scale_y_log10()
