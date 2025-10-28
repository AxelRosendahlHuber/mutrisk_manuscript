# Script focusing on the UKbiobank analyses
library(GenomicRanges)
library(rtracklayer)
library(cowplot)
library(UpSetR)
library(mutrisk)
source("code/functions/analysis_variables.R")

# internal function
## UKBiobank analyses
# load the bowel cancer data:
tissue = "colon"
crc_freq = fread("raw_data/UKBiobank/colorectal_cancer_frequency_UKB.csv")
colors = tissue_colors[[tissue]]

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

# groupp the UKBiobank cohort in groups of 10 years
ukbiobank_crc = ukbiobank_crc |>
  mutate(age_group = cut(age, seq(0, 80, 10), labels = seq(5, 75, 10))) |>
  group_by(age_group) |>
  summarize(risk = mean(risk),
            CRC_cumulative_risk = mean(CRC_cumulative_risk)) |>
  filter(!is.na(age_group)) |>
  dplyr::rename(age = age_group)
ukbiobank_crc$age = as.numeric(as.character(ukbiobank_crc$age))

# Compare the UK Biobank data to the globocan incidence
globocan_data = fread("raw_data/globocan_incidence/dataset-cumulative-risk-by-age-in-inc-both-sexes-age-0-84-in-2017-colorectum.csv")
globocan_data = globocan_data |> distinct() |>
  dplyr::rename(country = `Country label`, CRC_cumulative_risk = `Cumulative risk`) |>
  select(country, CRC_cumulative_risk) |>
  mutate(CRC_cumulative_risk = CRC_cumulative_risk / 100)

ukbiobank_crc_cumulative = ukbiobank_crc |>
  filter(age == 79) |>
  mutate(country = "UK - UKBiobank*\n UKB 79 years") |>
  select(country, CRC_cumulative_risk)

rbind(ukbiobank_crc_cumulative, globocan_data) |>
  ggplot(aes(x = country, y = CRC_cumulative_risk)) +
  geom_col() +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = NULL, y = "CRC cumulative risk", title = "Cumulative risk across countries",
       subtitle = "CRC cumulative risk at 84 years") +
  geom_text(aes(label = paste0(format(CRC_cumulative_risk*100, digits = 3), "%"), y = CRC_cumulative_risk),
            vjust = -0.2, position = position_dodge(0.9)) +
  theme_cowplot()

# load metadata
metadata = fread(paste0("processed_data/", tissue, "/", tissue, "_metadata.tsv")) |>   distinct() |>
  mutate(category = factor(category, levels = c("normal", "IBD", "POLD1", "POLE")))

# load the mutation rates
ncells = tissue_ncells_ci$mid_estimate[1]
expected_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
expected_rates_normal = expected_rates |>   filter(category == "normal")
ratios = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz")) |>
  filter(gene_name %in% c("APC", "KRAS", "TP53"))

# load the mutation data:
cancer_bDM = fread("processed_data/boostdm/boostdm_genie_cosmic/colon_boostDM_cancer.txt.gz")

# Calculate the driver mutation rates
KRAS_single_snv_muts = cancer_bDM[gene_name == "KRAS" & driver == TRUE, .N,
                                  c("gene_name", "mut_type", "aachange", "position")]

KRAS_single_snv = expected_rates_normal |>
  inner_join(KRAS_single_snv_muts, by = "mut_type", relationship = "many-to-many") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio)) |>
  group_by(donor, category) |>
  summarise(across(c(mle, cilow, cihigh, age), mean), groups = "drop")

# Expected number of cells with double mutations:
apc_counts_boostdm = cancer_bDM[gene_name == "APC" & driver == TRUE, .N, by = c("gene_name", "mut_type",  "driver")]

apc_single_snv = expected_rates_normal |>
  left_join(ratios |> filter(gene_name == "APC")) |>
  left_join(metadata) |>
  inner_join(apc_counts_boostdm |> filter(driver), by = "mut_type", relationship = "many-to-many") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio)) |>
  group_by(category, donor, age,  mut_type) |>
  summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop_last") |>
  summarize(across(c(mle, cilow, cihigh), sum), .groups = "drop")

double_apc = apc_single_snv |>
  mutate(across(c(mle, cilow, cihigh), ~ ((.^2) / 2)))

apc_kras = apc_single_snv |>
  mutate(mle = mle * KRAS_single_snv$mle,
         cilow = cilow * KRAS_single_snv$cilow,
         cihigh = cihigh * KRAS_single_snv$cihigh)

double_apc_kras = double_apc |>
  mutate(mle = mle * KRAS_single_snv$mle,
         cilow = cilow * KRAS_single_snv$cilow,
         cihigh = cihigh * KRAS_single_snv$cihigh)


# Plot the driver mutation rates
plot_driver_muts = function(driver_rates, y_axis = "INSERT TITLE") {
  driver_rates |>
    filter(category == "normal") |>
    mutate(across(c(mle, cilow, cihigh), ~ . * ncells)) |>
    ggplot(aes(x = age, y = mle, color = category)) +
    geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
    scale_y_continuous(labels = scales::label_comma()) +
    scale_color_manual(values = colors) +
    theme_cowplot() +
    labs(y = y_axis, x = "Age (years)") +
    theme(legend.position = "none")
}

F4B = plot_driver_muts(apc_single_snv, y_axis = "Number of cells with\nan APC driver mutation")
F4C = plot_driver_muts(KRAS_single_snv, y_axis = "Number of cells with\nKRAS driver mutations")
F4D = plot_driver_muts(double_apc, y_axis = "Number of cells with\ndouble APC driver mutations")
F4E = plot_driver_muts(apc_kras, y_axis = "Number of cells with\nAPC + KRAS driver mutations")
F4F = plot_driver_muts(double_apc_kras, y_axis = "Number of cells with\n doubleAPC + KRAS driver mutations") +
  scale_y_continuous(labels = scales::label_number_auto())

# make the figures for
figs = list(F4B = F4B, F4C = F4C, F4D = F4D, F4E = F4E, F4F = F4F)
annotated_figs = lapply(names(figs), \(x) prep_plot(figs[[x]], substr(x, 3,3)))

F3_middle = wrap_plots(annotated_figs, nrow = 1 )

# TP53 mutations - check if this needs to be here or in Figure 3 script (the TP53 script)
TP53_single_driver = expected_rates |>
  left_join(ratios |> filter(gene_name == "TP53")) |>
  left_join(metadata) |>
  inner_join(apc_counts_boostdm |> filter(driver), by = "mut_type", relationship = "many-to-many") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio)) |>
  group_by(category, donor, age,  mut_type) |>
  summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop_last") |>
  summarize(across(c(mle, cilow, cihigh), sum))


# make this plot for all tissues
TP53_single_driver_plot = TP53_single_driver |>
#  filter(category %in% c("normal", "POLD1")) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ncells),
         category = factor(category, levels = c("normal", "IBD", "POLD1", "POLE"))) |>
  ggplot(aes(x = age, y = mle, color = category)) +
  geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_color_manual(values = colors) +
  facet_grid(. ~ category ,axes = "all_y") +
  theme_cowplot() +
  labs(y = 'number of cells with TP53 mutations', x = "Age (years)")
TP53_single_driver_plot

# get the UKBiobank incidence
ukbiobank_crc_plot = ukbiobank_crc |>
  select(age, CRC_cumulative_risk) |>
  pivot_longer(CRC_cumulative_risk, names_to = "name", values_to = "incidence") |>
  mutate(category = "CRC Cumulative Incidence")

# Fraction of CRC mutated for APC or KRAS SNV combinations
mutated_fractions = fread("processed_data/GENIE_17/CRC_mutation_fractions.txt")

# make combinations of the different ages and CRC-mutation combinations
risk = expand.grid(percentages = mutated_fractions$percentages, CRC_cumulative_risk = ukbiobank_crc$CRC_cumulative_risk)
names = expand.grid(mutation_combination = mutated_fractions$mutation_combination, age = ukbiobank_crc$age)
mutated_fraction_CRC = cbind(names, risk) |> as_tibble() |>
 mutate(CRC_mut_fraction = percentages * CRC_cumulative_risk)

mutated_fraction_CRC |> filter(mutation_combination == "APC_double")
mutated_fraction_CRC |> filter(age == 75)

ggplot(mutated_fraction_CRC, aes(x = age, y = CRC_mut_fraction)) +
  geom_point() +
  geom_line() +
  ggh4x::facet_wrap2(mutation_combination ~ . , nrow = 1, axes = "y" ) +
  scale_y_continuous(breaks = scales::breaks_pretty(2), labels = label_percent()) +
  theme_cowplot() +
  labs(x = "Age (years)", y = "Lifetime risk for CRC with\nspecifc mutation")

#########################################################
# overlay the number of expected polyps in the epithelium
#########################################################
ad_incidence = fread("raw_data/polyp_incidence/Adenomas_by10_intestines.csv")
age_min_column = rep(c(40,50, 60, 70, 80), 2 )
age_max_column = rep(c(50, 60, 70, 80, 90), 2)

# first get the lesions by age, make horizontal bars to add in the text
# Mean Number of Adenomas per Ten Intestines Examined, by Age, Sex, and Adenoma Size
# Multiply by the fraction of adenoma's mutated for APC

# 135 CNADs sequenced, of which
73 / 135  # 54% of adenomas - check where this number is coming from - must be from the paper but
ad_incidence = ad_incidence |>
  filter(!Sex %in% c("Males Total", "Females Total")) |>
  mutate(min_age = age_min_column,
         max_age = age_max_column,
         age = (min_age + max_age) / 2,
         fraction_adenoma_apc = (`All sizes` * 0.54) / 10)

# take the mean across samples
ad_incidence_mean = ad_incidence |>
  group_by(age) |>
  summarize(fraction_adenoma_apc = mean(fraction_adenoma_apc))

ggplot(ad_incidence_mean, aes(x = age, y = fraction_adenoma_apc)) +
  geom_point() +
  geom_line() +
  theme_cowplot() +
  scale_y_continuous(breaks = scales::breaks_extended(5), labels = label_percent(), limits = c(0, .75)) +
  scale_x_continuous(limits = c(20, 85)) +
  labs(y = "Lifetime risk for Adenoma\nwith APC mutation",
       x = "Age (years)")