# UK biobank analyses
# what happens when we add extra colibactin to the samples?

# Script focusing on the UKbiobank analyses
library(cowplot)
library(UpSetR)
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

# Compare the UK Biobank data to the globocan incidence
globocan_data = fread("raw_data/globocan_incidence/dataset-cumulative-risk-by-age-in-inc-both-sexes-age-0-84-in-2017-colorectum.csv")
globocan_data = globocan_data |> distinct() |>
  dplyr::rename(country = `Country label`, CRC_cumulative_risk = `Cumulative risk`) |>
  select(country, CRC_cumulative_risk) |>
  mutate(CRC_cumulative_risk = CRC_cumulative_risk / 100)
ukbiobank_crc_cumulative = ukbiobank_crc |>
  filter(age == 79) |>
  mutate(country = "UK - UKBiobank*\n UKB 82 years\ninstead of 82") |>
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
tissue = "colon"
metadata = fread(paste0("processed_data/", tissue, "/", tissue, "_metadata.tsv")) |>   distinct() |>
  mutate(category = factor(category, levels = c("normal", "IBD", "POLD1", "POLE")))
# select only the information for specific donors
metadata_donors = metadata |>
  select(category, donor, age) |>
  distinct()

# load the mutation rates
ncells = tissue_ncells_ci$mid_estimate[2]
expected_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
ratios = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))
ratios = ratios |>
  filter(gene_name %in% c("TP53", "KRAS", "APC"))

# load signature rates
expected_rates_sig = fread("processed_data/colon/colon_sig_donor_rates.tsv.gz")
normal_donors = metadata |> filter(category == "normal") |> pull(donor) |> unique()
expected_rates_sig = expected_rates_sig |> filter(donor %in% normal_donors )

# add the max colibactin exposure to all other cells
max_sbs89_sbs88_exposures = expected_rates_sig |>
  filter(signature %in% c("SBS88", "SBS89")) |>
  group_by(signature, donor) |>
  summarise(mle = sum(mle), .groups = "drop_last") |>
  filter(mle == max(mle))

max_sbs89_sbs88_exposures = expected_rates_sig |>
  filter(signature %in% c("SBS88", "SBS89")) |>
  group_by(signature, donor) |>
  summarise(mle = sum(mle), .groups = "drop_last") |>
  arrange(desc(mle))

# donor O174 and donor PD37449 have the highest exposure of SBS89 and SBS88 respectively
SBS88_exposure = expected_rates_sig |>
  filter(signature == "SBS88" & donor == "PD37449")
SBS89_exposure = expected_rates_sig |>
  filter(signature == "SBS89" & donor == "O174")

# replace all values from SBS88 and SBS89 by the highest values:
# new "bm" - bacterial mutagenesis dataframe
expected_rates_sig_SBS88 = expected_rates_sig |>
  pivot_wider(names_from = signature, values_from = mle) |>
  mutate(SBS88 = rep(SBS88_exposure$mle, n_distinct(donor))) |>
  pivot_longer(cols = starts_with("SBS"), names_to = "signature", values_to = "mle")

expected_rates_sig_SBS89 = expected_rates_sig |>
  pivot_wider(names_from = signature, values_from = mle) |>
  mutate(SBS89 = rep(SBS89_exposure$mle, n_distinct(donor))) |>
  pivot_longer(cols = starts_with("SBS"), names_to = "signature", values_to = "mle")

# load the mutation data:
cancer_bDM = fread("processed_data/boostdm/boostdm_genie_cosmic/colon_boostDM_cancer.txt.gz")

# merge the tree sets of mutation data:
list = list(normal = expected_rates_sig,
            SBS88_max = expected_rates_sig_SBS88,
            SBS89_max = expected_rates_sig_SBS89)

total_rates_sig = rbindlist(list, idcol = "set", use.names = TRUE) |>
  mutate(mle = ifelse(is.na(mle), 0, mle))

total_rates = total_rates_sig |>
  group_by(set, donor, mut_type) |>
  summarize(mle = sum(mle))

# linear mutation accumulation in normal tissue
total_rates |> left_join(metadata |> select(category, donor, age) |> distinct()) |>
  filter(set == "normal") |>
  group_by(donor, category, age) |>
    summarize(mle = mean(mle)) |>
  ggplot(aes(x = age, y = mle)) +
  geom_point() +
  theme_cowplot() +
  labs(y = "mean mutation density/site")

# get KRAS driver mutation rate
KRAS_drivers = cancer_bDM[gene_name == "KRAS" & driver == TRUE , .N,
                                  c("gene_name", "mut_type", "aachange", "position")]

expected_rate_KRAS_single_snv_sc = total_rates |>
  inner_join(KRAS_drivers, relationship = "many-to-many", by = "mut_type") |>
  left_join(metadata_donors, relationship = "many-to-many", by = "donor") |>
  left_join(ratios, by = c("category", "gene_name")) |>
  mutate(mle = mle * ratio * N) |>
  group_by(set, donor, category) |>
  summarise(age = mean(age),
            mle = sum(mle), .groups = "drop") |>
  arrange(set, category, donor, age)

# Expected number of cells with double mutations:
apc_counts_boostdm = cancer_bDM[gene_name == "APC" & driver == TRUE, .N, by = c("gene_name", "mut_type")]

apc_single_driver = total_rates |>
  inner_join(apc_counts_boostdm, by = "mut_type") |>
  left_join(metadata_donors, relationship ="many-to-many", by = "donor") |>
  left_join(ratios, by = c("category", "gene_name")) |>
  mutate(mle = mle * ratio * N) |>
  group_by(set, donor, category) |>
  summarise(age = mean(age),
            mle = sum(mle)) |>
  arrange(set, category, donor, age)

# extrapolate the risk for double driver
apc_double_driver = apc_single_driver |>
  mutate(mle =   ((mle^2) / 2))

# make a dataframe with the relative mutation rates for all the different changes:
vogelgram_mut_risks_sc = tibble(
  set = apc_double_driver$set,
  category = apc_double_driver$category,
  donor = apc_double_driver$donor,
  age = apc_double_driver$age,
  `APC_single_snv` = apc_single_driver$mle,
  `APC_double` = apc_double_driver$mle,
  KRAS_single_snv = expected_rate_KRAS_single_snv_sc$mle) |>
  mutate(`KRAS_APC_double` = `APC_double` * KRAS_single_snv,
         `KRAS_APC_single` = `APC_single_snv` * KRAS_single_snv)

vogelgram_mut_risks = vogelgram_mut_risks_sc |>
  mutate(across(-c(set, category, donor, age), ~ . * ncells))

vogelgram_mut_risks_long = vogelgram_mut_risks |>
  filter(category == "normal") |>
  pivot_longer(-c(set, category, donor, age), names_to = "mutation", values_to = "expected_mutated_cells")

ukbiobank_crc_plot = ukbiobank_crc |>
  filter(n_tumor != 0) |>
  select(age, CRC_cumulative_risk) |>
  pivot_longer(CRC_cumulative_risk, names_to = "mutation", values_to = "incidence")

vogelgram_risks_long = vogelgram_mut_risks_long |>
  filter(category == "normal") |>
  dplyr::rename(incidence = expected_mutated_cells) |>
  select(set, category, age, donor, mutation, incidence)

vogelgram_risks_long |>
  ggplot(aes(x = age, y = incidence, color = mutation, shape = set)) +
  geom_point() +
  theme_cowplot()

# make plot with all changes
plot_colibactin = vogelgram_risks_long |>
  ggplot(aes(x = age, y = incidence, color = mutation, shape = set)) +
  geom_point() + scale_y_continuous(labels = scales::label_comma())
plot_colibactin

plot_colibactin + scale_y_log10()
plot_colibactin +   facet_wrap(mutation ~ ., scales = "free_y")

plot_colibactin = vogelgram_risks_long |>
  ggplot(aes(x = age, y = incidence, color = mutation, shape = set)) +
  geom_point()
plot_colibactin

# make plot which shows the net-increase, and fold-increase
risks = vogelgram_risks_long |>
  select(-category) |>
  pivot_wider(names_from = set, values_from = incidence)
relative_risks_mutrate = risks |>
  mutate(rel_SBS88 = SBS88_max / normal,
         rel_SBS89 = SBS89_max / normal)

relative_risks_mutrate |>
  ggplot(aes(x = age, y = rel_SBS88)) +
  geom_point() +
  facet_wrap(mutation ~ .) +
  theme_bw() +
  geom_hline(yintercept = 1)

relative_risks_mutrate |>
  ggplot(aes(x = age, y = rel_SBS89)) +
  geom_point() +
  facet_wrap(mutation ~ .) +
  theme_bw() +
  geom_hline(yintercept = 1)

table(metadata |> filter(category == "normal") |> select(donor, age) |> distinct() |> pull(age))

## get-signature specific rates:
total_rates_sig

er_KRAS_single_snv_sig = total_rates_sig |>
  inner_join(KRAS_drivers, relationship = "many-to-many", by = "mut_type") |>
  left_join(metadata_donors, relationship = "many-to-many", by = "donor") |>
  left_join(ratios, by = c("category", "gene_name")) |>
  mutate(mle = mle * ratio * N) |>
  group_by(set, donor, category, signature) |>
  summarise(age = mean(age),
            mle = sum(mle), .groups = "drop") |>
  arrange(set, category, donor, age) |>
  filter(mle != 0)

# Expected number of cells with double mutations:
apc_single_driver = total_rates_sig |>
  inner_join(apc_counts_boostdm, by = "mut_type") |>
  left_join(metadata_donors, relationship = "many-to-many", by = "donor") |>
  left_join(ratios, by = c("category", "gene_name")) |>
  mutate(mle = mle * ratio * N) |>
  group_by(set, donor, category, signature) |>
  summarise(age = mean(age),
            mle = sum(mle)) |>
  arrange(set, category, donor, age) |>
  filter(mle != 0)

# extrapolate the risk for double driver
apc_double_driver = apc_single_driver |>
  mutate(mle = ((mle^2) / 2))

vogelgram_mut_risks_sc = tibble(
  set = apc_double_driver$set,
  category = apc_double_driver$category,
  donor = apc_double_driver$donor,
  age = apc_double_driver$age,
  `APC_single_snv` = apc_single_driver$mle,
  `APC_double` = apc_double_driver$mle,
  signature = apc_double_driver$signature,
  KRAS_single_snv = er_KRAS_single_snv_sig$mle) |>
  mutate(`KRAS_APC_double` = `APC_double` * KRAS_single_snv,
         `KRAS_APC_single` = `APC_single_snv` * KRAS_single_snv)

vogelgram_mut_risks_long = vogelgram_mut_risks_sc |>
  pivot_longer(cols = contains("_"), names_to = "mutation", values_to = "incidence")

vogelgram_mut_risks_long |>
  mutate(signature_col = ifelse(signature %in% c("SBS89", "SBS88"), ""))

# okay this does not work for signatures: by multiplying the risks it is not possible to achieve this

# barplots which show the increase during age

# plot the mutation risk for individual cells:
df_vogelgram_incidence = rbind(vogelgram_risks_long, ukbiobank_crc_plot) |>
  mutate(name = factor(name, levels = c("APC_single_snv", "KRAS_single_snv", "KRAS_APC_single",
                                        "APC_double", "KRAS_APC_double", "CRC_cumulative_risk")))


df_vogelgram_incidence = df_vogelgram_incidence |>
  filter(name != "APC_single_snv")

vogel_plot = df_vogelgram_incidence |>
  filter(name != "APC_single_snv") |>
  ggplot(aes(x = age)) +
  geom_point(aes(x = age, y = incidence, color = name)) +
  scale_y_log10(guide = "axis_logticks",
                breaks = c(1e-4, 1e-2,  1, 100, 1e4, 1e6),
                labels = ~ ifelse(.x < 1, scales::comma(.x), scales::comma(.x, accuracy = 1))) +
  scale_color_manual(values = c("darkgoldenrod1",  "darkorange3", "red", "darkred", "black")) +
  theme_cowplot()  +
  theme(legend.position = "inside", legend.position.inside = c(0.6, 0.1), legend.background = element_rect(fill = alpha("white", 0.7), color = "black"),
        legend.margin = margin(1,1,1,1,"mm")) +
  coord_cartesian(clip = 'off') +
  labs(y = "Estimated risk", x = "Age (years)", color = NULL)
vogel_plot


df_vogelgram5.5_incidence = df_vogelgram_incidence |>
  mutate(incidence = case_when(
    name == "CRC_cumulative_risk" ~ incidence * 0.0552, .default = incidence))
vogel_plot_5.5 = df_vogelgram5.5_incidence |>
  filter(name != "APC_single_snv") |>
  ggplot(aes(x = age)) +
  geom_point(aes(x = age, y = incidence, color = name)) +
  scale_y_log10(guide = "axis_logticks",
                breaks = c(1e-4, 1e-2,  1, 100, 1e4, 1e6),
                labels = ~ ifelse(.x < 1, scales::comma(.x), scales::comma(.x, accuracy = 1))) +
  scale_color_manual(values = c("darkgoldenrod1",  "darkorange3", "red", "darkred", "black")) +
  theme_cowplot()  +
  theme(legend.position = "inside", legend.position.inside = c(0.6, 0.1), legend.background = element_rect(fill = alpha("white", 0.7), color = "black"),
        legend.margin = margin(1,1,1,1,"mm")) +
  coord_cartesian(clip = 'off') +
  labs(y = "Incidence", x = "Age (years)", color = NULL)
vogel_plot_5.5