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
ratios = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz")) |>
  filter(gene_name %in% c("APC", "KRAS"))

# load the mutation data:
cancer_bDM = fread("processed_data/boostdm/boostdm_genie_cosmic/colon_boostDM_cancer.txt.gz")

# get KRAS driver mutation rate
KRAS_single_snv_muts = cancer_bDM[gene_name == "KRAS" & driver == TRUE, .N,
                                  c("gene_name", "mut_type", "aachange", "position")]
KRAS_G12V_hotspot = KRAS_single_snv_muts

expected_rate_KRAS_single_snv_sc = expected_rates |>
  inner_join(KRAS_G12V_hotspot, by = "mut_type", relationship = "many-to-many") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio)) |>
  group_by(donor, category) |>
  summarise(across(c(mle, cilow, cihigh, age), mean)) |>
  arrange(category, donor, age)

# Expected number of cells with double mutations:
apc_counts_boostdm = cancer_bDM[gene_name == "APC" & driver == TRUE, .N, by = c("gene_name", "mut_type",  "driver")]

apc_single_driver = expected_rates |>
  left_join(ratios |> filter(gene_name == "APC")) |>
  left_join(metadata) |>
  inner_join(apc_counts_boostdm |> filter(driver), by = "mut_type", relationship = "many-to-many") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio)) |>
  group_by(category, donor, age,  mut_type) |>
  summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop_last") |>
  summarize(across(c(mle, cilow, cihigh), sum))

apc_double_driver = apc_single_driver |>
  mutate(across(c(mle, cilow, cihigh), ~ ((.^2) / 2)))

APC_single_mut_plot = apc_single_driver |>
  filter(category == "normal") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ncells)) |>
  ggplot(aes(x = age, y = mle, color = category)) +
  geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
  scale_y_continuous(labels = scales::label_comma()) +
  scale_color_manual(values = colors) +
  theme_cowplot() +
  labs(y = 'number of cells with APC mutations', x = "Age (years)") +
  theme(legend.position = "none")

APC_double_mut_plot = apc_double_driver |>
  filter(category == "normal") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ncells)) |>
  ggplot(aes(x = age, y = mle, color = category)) +
  geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
  theme_cowplot() +
  scale_color_manual(values = colors) +
  labs(y = 'number of cells w/\nAPC mutations on both alleles', x = "Age (years)") +
  theme(legend.position = "none")

figure2BC = APC_single_mut_plot + APC_double_mut_plot
figure2BC
ggsave("plots/colon/Figure2BC.png", figure2BC, width = 8.5, height = 4)

# make a dataframe with the relative mutation rates for all the different changes:
vogelgram_mut_risks_sc = tibble(
  category = apc_double_driver$category,
  donor = apc_double_driver$donor,
  age = apc_double_driver$age,
  `APC_single_snv` = apc_single_driver$mle,
  `APC_double` = apc_double_driver$mle,
  KRAS_single_snv = expected_rate_KRAS_single_snv_sc$mle) |>
  mutate(`KRAS_APC_double` = `APC_double` * KRAS_single_snv,
         `KRAS_APC_single` = `APC_single_snv` * KRAS_single_snv)

vogelgram_mut_risks = vogelgram_mut_risks_sc |>
  mutate(across(-c(category, donor, age), ~ . * ncells))

vogelgram_mut_risks_long = vogelgram_mut_risks |>
  pivot_longer(-c(category, donor, age), names_to = "mutation", values_to = "expected_mutated_cells")

ukbiobank_crc_plot = ukbiobank_crc |>
  select(age, CRC_cumulative_risk) |>
  pivot_longer(CRC_cumulative_risk, names_to = "name", values_to = "incidence") |>
  mutate(category = "CRC Cumulative Incidence")

vogelgram_risks_long = vogelgram_mut_risks_long |>
  dplyr::rename(name = mutation, incidence = expected_mutated_cells) |>
  select(category, age, name, incidence)

# plot the mutation risk for individual cells:
df_vogelgram_incidence = rbind(vogelgram_risks_long, ukbiobank_crc_plot) |>
  mutate(name = factor(name, levels = c("APC_single_snv", "KRAS_single_snv", "KRAS_APC_single",
                                        "APC_double", "KRAS_APC_double", "CRC_cumulative_risk")))

df_vogelgram_incidence = df_vogelgram_incidence |>
  filter(name != "APC_single_snv")

# remake the vogelgram plot so that it is faceted for the different tissues, and has a different multiplier
# for each of the occurrences

# new script to compare the incidence of the vogelgram plot
mutated_fractions = fread("processed_data/GENIE_17/CRC_mutation_fractions.txt") |>
  as.data.frame() |>
  column_to_rownames("mutation_combination")

# make plots for normal individuals:
# make in a for loop multiple plots - with for each of them the %age of CRC indicated
i = 1 # for testing purposes only
plot_multiple_fraction_list = plot_fraction_list = plot_list = list(spacer = plot_spacer())
for (i in 1:nrow(mutated_fractions)) {
  print(i)

  i_name = rownames(mutated_fractions)[i]
  ukbiobank_crc_plot_i = ukbiobank_crc_plot |>
    mutate(incidence = incidence * mutated_fractions[i, 1])

  vogelgram_risks_i = vogelgram_risks_long |>
    filter(name == i_name) |>
    filter(category == "normal")

  # plot the mutation risk for individual cells:
  df_vogelgram_incidence = rbind(vogelgram_risks_i, ukbiobank_crc_plot_i) |>
    mutate(name = factor(name, levels = c("APC_single_snv", "KRAS_single_snv", "KRAS_APC_single",
                                          "APC_double", "KRAS_APC_double", "CRC_cumulative_risk")))

  # add fold change to the plots:
  vogelgram_risks_i_old = df_vogelgram_incidence |> filter(age > 40 & name != "CRC_cumulative_risk")
  mean_incidence_ukbiobank = ukbiobank_crc_plot_i |>
    filter(age %in% vogelgram_risks_i_old$age) |>
    pull(incidence) |> mean()

  mean_incidence_clones = vogelgram_risks_i_old |>
    pull(incidence) |> mean()

  ratio_mut_cells_crc = mean_incidence_clones / mean_incidence_ukbiobank

  plot_list[[i_name]] = df_vogelgram_incidence |>
    ggplot(aes(x = age)) +
    geom_point(aes(x = age, y = incidence, color = name)) +
    scale_y_log10(guide = "axis_logticks",
                  breaks = c(1e-4, 1e-2,  1, 100, 1e4, 1e6),
                  labels = ~ ifelse(.x < 1, scales::comma(.x), scales::comma(.x, accuracy = 1))) +
    scale_color_manual(values = c(c("darkgoldenrod1",  "darkorange3", "red", "darkred")[i], "black")) +
    theme_cowplot()  +
    theme(legend.position = "inside", legend.position.inside = c(0.5, 0.1), legend.background = element_rect(fill = alpha("white", 0.7)),
          legend.margin = margin(1,1,1,1,"mm")) +
    coord_cartesian(clip = 'off') +
    labs(y = "Incidence", x = "Age (years)", color = NULL, title = i_name,
         subtitle = paste0("Percentage of CRC = ",format(mutated_fractions[i, 1]*100, digits = 3),
                           "%\nRatio mutated cells:CRC = ", format(ratio_mut_cells_crc, digits = 3, big.mark = ",")))

  # plot the individual fractions summing up to 1
  df_vogelgram_fraction = df_vogelgram_incidence |>
    mutate(incidence = ifelse(name != "CRC_cumulative_risk",
                              get_prob_mutated_N(incidence/ncells, ncells = ncells),
                              incidence))

  plot_fraction_list[[i_name]] =  df_vogelgram_fraction |>
    ggplot(aes(x = age)) +
    geom_point(aes(x = age, y = incidence, color = name)) +
    scale_y_log10(guide = "axis_logticks",
                  breaks = c(1e-4, 1e-2,  1, 100, 1e4, 1e6),
                  labels = ~ ifelse(.x < 1, scales::comma(.x), scales::comma(.x, accuracy = 1))) +
    scale_color_manual(values = c(c("darkgoldenrod1",  "darkorange3", "red", "darkred")[i], "black")) +
    theme_cowplot()  +
    theme(legend.position = "inside", legend.position.inside = c(0.5, 0.2), legend.background = element_rect(fill = alpha("white", 0.7)),
          legend.margin = margin(1,1,1,1,"mm")) +
    coord_cartesian(clip = 'off') +
    labs(y = "Incidence mutation/CRCR cancer", x = "Age (years)", color = NULL, title = i_name,
         subtitle = paste0("Percentage of CRC = ",format(mutated_fractions[i, 1]*100, digits = 3),
                           "%\nRatio mutated cells:CRC = ", format(ratio_mut_cells_crc, digits = 3, big.mark = ",")))

  # plot which fraction of cells has multiple mutations:
  df_vogelgram_multiple_fraction = df_vogelgram_incidence |>
    filter(name != "CRC_cumulative_risk")

  fraction_list = list()
  for (n in c(1, 2,5, 10, 100)) {
    fraction_list[[as.character(n)]] = df_vogelgram_multiple_fraction |>
      mutate(fraction = get_prob_mutated_N(incidence/ncells, ncells, n))
  }
  df_vogelgram_multiple_fraction = rbindlist(fraction_list, idcol = "N")

  df_all_multiple_fraction = df_vogelgram_incidence |>
    filter(name == "CRC_cumulative_risk") |>
    mutate(N = "CRC cumulative risk",
           fraction = incidence) |>
    bind_rows(df_vogelgram_multiple_fraction) |>
    mutate(N = as.character(N),
      N = factor(N, levels = c("1", "2", "5", "10", "100", "CRC cumulative risk")))

  plot_multiple_fraction_list[[i_name]] =  df_all_multiple_fraction |>
    ggplot(aes(x = age)) +
    geom_point(aes(x = age, y = fraction, color = N)) +
    scale_color_manual(values = c(c("darkgoldenrod1",  "darkorange3", "firebrick1", "red", "darkred"), "black")) +
    theme_cowplot()  +
    theme(legend.position = "inside", legend.position.inside = c(0.5, 0.2), legend.background = element_rect(fill = alpha("white", 0.7)),
          legend.margin = margin(1,1,1,1,"mm")) +
    coord_cartesian(clip = 'off') +
    labs(y = "Incidence mutation/CRCR cancer", x = "Age (years)", color = NULL, title = i_name,
         subtitle = paste0("Percentage of CRC = ",format(mutated_fractions[i, 1]*100, digits = 3),
                           "%\nRatio mutated cells:CRC = ", format(ratio_mut_cells_crc, digits = 3, big.mark = ",")))



}


#########################################################
# overlay the number of expected polyps in the epithelium
#########################################################
hp_incidence = fread("raw_data/polyp_incidence/HP_by10_intestines.csv")
ad_incidence = fread("raw_data/polyp_incidence/Adenomas_by10_intestines.csv")
polypoid_lesions = fread("raw_data/polyp_incidence/Polyploid_lesions.csv")
age_min_column = rep(c(40,50, 60, 70, 80), 2 )
age_max_column = rep(c(50, 60, 70, 80, 90), 2)

# first get the lesions by age, make horizontal bars to add in the text
# Mean Number of Adenomas per Ten Intestines Examined, by Age, Sex, and Adenoma Size
# Multiply by the fraction of adenoma's mutated for APC

# 135 CNADs sequenced, of which
73 / 135  # 54% of adenomas
ad_incidence = ad_incidence |>
  filter(!Sex %in% c("Males Total", "Females Total")) |>
  mutate(min_age = age_min_column,
         max_age = age_max_column,
         fraction_adenoma_apc = (`All sizes` * 0.54) / 10)
ad_incidence_mean = ad_incidence |>
  group_by(min_age, max_age) |>
  summarize(fraction_adenoma_apc = mean(fraction_adenoma_apc))

incidence_list = list()
for (i in rownames(mutated_fractions)) {
  ad = ad_incidence_mean
  ad$fraction_adenoma_apc = ad_incidence_mean$fraction_adenoma_apc * mutated_fractions[i, 1]
  incidence_list[[i]] = ad
}
ad_incidence_mean = rbindlist(incidence_list, idcol = "mutation")

ad_incidence_mean = ad_incidence_mean |>
  mutate(fraction_adenoma_apc = get_prob_mutated_N(fraction_adenoma_apc, 1))

# Prevalence of Hyperplastic Polyps (HP) and Mean Number of HP per Ten   Intestines Examined, by Age and Sex
# As these values are estimates of the actual number of individuals with these HP numbers,
# it is better to continue for now with the adenoma numbers
polypoid_lesions = polypoid_lesions |>
  filter(!Sex %in% c("Males Total", "Females Total")) |>
  mutate(min_age = age_min_column) |>
  mutate(max_age = age_max_column) |>
  mutate(sum_polyps = `1` +  (`2–4` * 3) + (`5–9` * 7) + (`10+` * 10))

###########################################
# Study the incidence rates of more than 1 mut
# using the probabilities - for 1, 10 and 100 muts
##########################################
# less mutagenic sites:
vogelgram_mut_prob_list = list()
for (i in c(1, 4, 10, 30, 100, 300, 1000, 3000, 25000)) {
  name = as.character(i)
  vogelgram_mut_prob_list[[i]] = vogelgram_mut_risks_sc |>
    mutate(across(-c(category, donor, age), ~ get_prob_mutated_N(., ncells, i)))
}

vogelgram_mut_probs = rbindlist(vogelgram_mut_prob_list, idcol = "Nmut")
vogelgram_mut_probs_long = vogelgram_mut_probs |>
  pivot_longer(contains("_"), names_to = "mutation", values_to = "mut_probability") |>
  mutate(nmut = factor(Nmut, sort(unique(Nmut))),
         mutation = factor(mutation, levels = c("APC_single_snv", "KRAS_single_snv", "APC_double", "KRAS_APC_single",
                                                "KRAS_APC_double"))) |>
  filter(category == "normal")



vogelgram_mut_plot = vogelgram_mut_probs_long |>
  filter(mutation != "APC_single_snv" & Nmut < 100 |
           mutation == "APC_single_snv" & Nmut > 100) |>  # filter for the less mutagenic sites
  mutate(mutname = gsub("_", " ", mutation),
         mutname = gsub("single", "1", mutname),
         mutname = gsub("snv", "SNV", mutname),
         mutname = case_match(mutname, "APC double" ~ "APC 2 SNV",
                              "KRAS APC double" ~ "KRAS 1 SNV APC 2 SNV",
                              "KRAS APC 1" ~ "KRAS 1 SNV APC 1 SNV",
                              .default = mutname))

# filter plots:
ad_incidence_mean  = ad_incidence_mean |>
  mutate(mutname = gsub("_", " ", mutation),
         mutname = gsub("single", "1", mutname),
         mutname = gsub("snv", "SNV", mutname),
         mutname = case_match(mutname, "APC double" ~ "APC 2 SNV",
                              "KRAS APC double" ~ "KRAS 1 SNV APC 2 SNV",
                              "KRAS APC 1" ~ "KRAS 1 SNV APC 1 SNV",
                              .default = mutname))

ggplot(vogelgram_mut_plot, aes(x = age, y = mut_probability)) +
  geom_smooth(aes( color = nmut), formula = y ~ x,se = FALSE,
              method = "glm", fullrange = TRUE,
              method.args = list(family = quasibinomial(link = "probit"))) +
#  geom_point(aes( color = nmut)) +
  geom_smooth(data = ad_incidence_mean,
              aes(x = (min_age + max_age)/2,
                  y = fraction_adenoma_apc),
              formula = y ~ x,se = FALSE,
              method = "glm", fullrange = TRUE, color = "black",
              method.args = list(family = quasibinomial(link = "probit"))) +
    geom_point(data = ad_incidence_mean,
               aes(x = (min_age + max_age)/2,
                   y = fraction_adenoma_apc), color = "black") +
  scale_color_manual(values = c('#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850')) +
  facet_wrap(mutname ~ . , scale = "free_y", axes = "all") +
  theme_cowplot() +
  labs(color = "Number mutated cells", y = "Probability of event occurring\nin individual", x = "Age (Years)") +
  theme(strip.background = element_blank())

# adding both to the plot:
ukbiobank_list = list()
for (i in rownames(mutated_fractions)) {
  ukbiobank_list[[i]] = ukbiobank_crc_plot |>
    mutate(incidence = incidence * mutated_fractions[i, 1])
}

ukbiobank = rbindlist(ukbiobank_list, idcol = "mutation") |>
  mutate(mutname = gsub("_", " ", mutation),
         mutname = gsub("single", "1", mutname),
         mutname = gsub("snv", "SNV", mutname),
         mutname = case_match(mutname, "APC double" ~ "APC 2 SNV",
                              "KRAS APC double" ~ "KRAS 1 SNV APC 2 SNV",
                              "KRAS APC 1" ~ "KRAS 1 SNV APC 1 SNV",
                              .default = mutname))

vogelgram_mut_plot |>
  mutate(Nmut = factor(Nmut, levels = (unique(Nmut)))) |>
  filter(Nmut %in% c(1)) |>
  ggplot(aes(x = age, y = mut_probability)) +
  geom_smooth(aes( color = nmut), formula = y ~ x,se = FALSE,
              method = "glm", fullrange = TRUE,
              method.args = list(family = quasibinomial(link = "probit"))) +

  geom_point(data = ukbiobank, aes(x = age, y = incidence)) +
  geom_smooth(data = ukbiobank,
              aes(x = age, y = incidence),
              formula = y ~ x,se = FALSE,
              method = "glm", fullrange = TRUE, color = "black",
              method.args = list(family = quasibinomial(link = "probit"))) +
  scale_color_manual(values = c('#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850')) +
  facet_wrap(mutname ~ . , scale = "free_y", axes = "all") +
  theme_cowplot() +
  labs(color = "Number mutated cells", y = "Probability of event occurring\nin individual", x = "Age (Years)") +
  theme(strip.background = element_blank())
ggsave("plots/")

