# Figure_4_blood new
source("code/functions/analysis_variables.R")
library(patchwork)
library(cowplot)
library(geomtextpath)
tissue = "blood"
blood_colors = "#ff725c"

# load blood metadata
ncells = 1e5
tissue = "blood"
metadata = fread(paste0("processed_data/", tissue, "/", tissue, "_metadata.tsv")) |>
  distinct()

# load the mutation rates
expected_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))

# load the relative mutation ratio for each sample for the DNMT3A gene:
ratios = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))
sig_donor_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_sig_donor_rates.tsv.gz"))

# load the GENIE data
genie_blood = fread("processed_data/GENIE_17/GENIE_17_processed.txt.gz") |>
  filter(CANCER_TYPE %in% c("Leukemia","B-Lymphoblastic Leukemia/Lymphoma","Myeloproliferative Neoplasms",
                            "Myelodysplastic Syndromes","Mature T and NK Neoplasms",
                            "Myelodysplastic/Myeloproliferative Neoplasms"))

# load boostdm_ch-genie-cosmic intersections
CH_bDM = fread("processed_data/boostdm/boostdm_genie_cosmic/CH_boostDM_cancer.txt.gz")

##### PLOTTING #####
# check which mutation is the R882H frequently mutated hotspot
####
# changes to make the script more accommodating for multiple different tissues:
select_gene = "KRAS"

# calculate expected rates
calc_exp_muts = function(expected_rates, mut_positions, metadata, ratios, ncells) {
  expected_rates |>
    inner_join(mut_positions, by = "mut_type", relationship = "many-to-many") |>
    left_join(metadata, by = c("sampleID", "coverage", "category")) |>
    left_join(ratios, by = c("category", "gene_name")) |>
    mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells * N)) |>
    group_by(donor, category, sampleID) |>
    summarise(across(c(mle, cilow, cihigh), sum),
              age = mean(age), .groups = "drop_last") |>
    summarise(across(c(mle, cilow, cihigh, age), mean), .groups = "drop")
}

DNMT3A_R882H_hotspot = CH_bDM[aachange == "R882H" & gene_name == "DNMT3A", .N, c("gene_name", "mut_type", "aachange", "position", "driver")]
DNMT3A_drivers = CH_bDM[gene_name == "DNMT3A" & driver == TRUE , .N, c("gene_name", "mut_type", "aachange", "position", "driver")]

DNMT3A_R882H_rates = calc_exp_muts(expected_rates, DNMT3A_R882H_hotspot, metadata= metadata, ratios = ratios, ncells = ncells)
DNMT3A_driver_rates = calc_exp_muts(expected_rates, DNMT3A_drivers, metadata= metadata, ratios = ratios, ncells = ncells)

F5C = ggplot(DNMT3A_driver_rates,
             aes(x = age, y = mle)) +
  geom_pointrange(aes(ymin = mle/4, ymax = mle*13), color = blood_colors)  +
  labs(y = "DNMT3A driver mutations", x = "Age (years)") +
  theme_cowplot()

F5D = ggplot(DNMT3A_R882H_rates,
       aes(x = age, y = mle)) +
  geom_pointrange(aes(ymin = mle/4, ymax = mle*13), color = blood_colors)  +
  labs(y = "DNMT3A R882H mutations", x = "Age (years)") +
  theme_cowplot()

# get all the driver mutations for the watson figure:
watson_variants = c("R882C", "R729W", "R326C", "R320*", "R882H", "R736H",
     "Y735C", "R736C", "W860R", "R771*", "R598*", "P904L")

DNMT3A_watson_drivers = CH_bDM[gene_name == "DNMT3A" & driver == TRUE  &
                                 aachange %in% watson_variants, .N, c("gene_name", "mut_type", "aachange", "position", "driver")]

mutation_list = list(
  DNMT3A_R882H = calc_exp_muts(expected_rates, DNMT3A_R882H_hotspot, metadata, ratios, ncells),
  DNMT3A_drivers = calc_exp_muts(expected_rates, DNMT3A_drivers, metadata, ratios, ncells),
  DNMT3A_watson_drivers = calc_exp_muts(expected_rates, DNMT3A_watson_drivers, metadata, ratios, ncells)) |>
  rbindlist(idcol = "name") |>
  mutate(name = factor(name, levels = c("DNMT3A_drivers", "DNMT3A_R882H", "DNMT3A_watson_drivers")),
        prob_mle = get_prob_mutated_N(risk = mle/ncells ,ncells = ncells, N =  1),
        prob_cilow = get_prob_mutated_N(risk = (mle)/ncells ,ncells = ncells/4, N =  1),
        prob_cihigh = get_prob_mutated_N(risk = (mle)/ncells ,ncells = ncells*4, N =  1),
        prob_mut_5 = get_prob_mutated_N(risk = mle/ncells ,ncells = ncells, N =  5),
        prob_mut_10 = get_prob_mutated_N(risk = mle/ncells ,ncells = ncells, N = 13))

F5E = ggplot(mutation_list |> filter(name %in% c("DNMT3A_R882H", "DNMT3A_drivers")),
       aes(x = age, y = prob_mle)) +
  geom_ribbon(aes(ymin = prob_cilow, ymax = prob_cihigh), alpha = 0.5, fill = blood_colors) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_smooth(formula = y ~ x,se = FALSE,
              method = "glm", fullrange = TRUE,
              method.args = list(family = quasibinomial(link = "probit")),
              color = blood_colors) +
  geom_point() +
  scale_fill_manual(values = blood_colors) +
  scale_color_manual(values = blood_colors) +
  theme_cowplot() +
  labs(y = "number of cells/individual",
       x = "Age (years)",fill = NULL) +
  theme(legend.position = "none")
F5E

# specifically select against sites in DNMT3A:
UKB_DNMT3A = fread("raw_data/UKBiobank/UKB_age_frequencies_DNMT3A.tsv")
UKB_DNMT3A = UKB_DNMT3A |>
  filter(Individuals > 5000) |> # filter for the years when at least 5000 patients are in the cohort
  mutate(R882H = `R/H` / Individuals,
         R882C = `R/C` / Individuals,
         R882S = `R/S` / Individuals,
         R882P = `R/P` / Individuals)


UKB = fread("raw_data/UKBiobank/UkBiobank_DNMT3A_mut_age.csv")
UKB = UKB |> dplyr::count(aa_change, Age) |>
  left_join(UKB_DNMT3A |> select(Age, Individuals)) |>
  filter(Individuals > 5000)

UKB = UKB |> dplyr::count(aa_change, Age) |>
  left_join(UKB_DNMT3A |> select(Age, Individuals)) |>
  filter(Individuals > 5000)

UKB_fraction = UKB |>
  mutate(fraction = n / Individuals)

UKB_fraction_all = UKB_fraction |>
  mutate(type = "UKB all DNMT3A") |>
  group_by(Age, type) |>
  summarize(fraction = sum(fraction)) |>
  select(Age, fraction, type)

curve_data = mutation_list |> filter(name == "DNMT3A_drivers") |>
  dplyr::rename(mutation = name, Age = age) |>
  pivot_longer(cols = c("prob_mle", "prob_mut_5", "prob_mut_10")) |>
  group_by(name) |>
  mutate(label = if_else(value == max(value), as.character(name), NA_character_),
         label = case_when(label == "prob_mle" ~ "> 1 DNMT3A driver",
                           label == "prob_mut_5" ~ "> 5 DNMT3A drivers",
                           label == "prob_mut_10" ~ "> 10 DNMT3A drivers",
                            .default = label))

F5F = ggplot(curve_data, aes(x = Age, y = value)) +
  geom_smooth(aes(group = name), formula = y ~ x,se = FALSE,
              method = "glm",               method.args = list(family = quasibinomial(link = "probit")),
              color = blood_colors) +
  geom_point(data = UKB_fraction_all, mapping = aes(y = fraction))  +
  ggrepel::geom_text_repel(aes(label = label),
                   nudge_x = 3, min.segment.length = 100,
                   na.rm = TRUE) +
  scale_fill_manual(values = ) +
  scale_color_manual(values = blood_colors) +
  scale_x_continuous(limits = c(0, 90)) +
  theme_cowplot() +
  labs(y = "number of cells/individual",
       x = "Age (years)",fill = NULL) +
  theme(legend.position = "none")
F5F


# make figure
F5C = prep_plot(F5C,label = "C", t = 8, r = 8, l = 8, b = 8)
F5D = prep_plot(F5D,label = "D", t = 8, r = 8, l = 8, b = 8)
F5E = prep_plot(F5E,label = "E", t = 8, r = 8, l = 8, b = 8)
F5F = prep_plot(F5F,label = "F", t = 8, r = 8, l = 8, b = 8)

F5C + F5D + F5E + F5F

# for the fitness effect: we need to consider that we can only model mutations with a specific fitness effect
# perform two types of analyses:s
# 1. with the DNMT3A mutation variant
# 2. with all DNMT3A variants. Check how we need to calculate the fitness. Probably taking the sum of all of the curves would be the best solution

# also needed would be the list of all the variants in the Watson analysis reported to be mutated

# 1: R882H Plot:
mutation_list_ci_prob = mutation_list_ci |>
  mutate(cilow =get_prob_mutated_N( cilow / 2.5e4, 2.5e4),
         mle =get_prob_mutated_N(mle / 1e5,  1e5),
         cihigh =get_prob_mutated_N( cihigh / 1.3e6, 1.3e6))

figure_4b = ggplot(mutation_list_ci_prob,
                   aes(x = age, y = mle, fill = category)) +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh, color = category),
                width = 0,
                show.legend = FALSE) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_point(shape = 21, size = 2.4, color = "white", stroke = 0.3) +
  scale_fill_manual(values = blood_colors) +
  scale_color_manual(values = blood_colors) +
  theme_cowplot() +
  labs(y = "number of cells/individual",
       x = "Age (years)",fill = NULL) +
  theme(legend.position = "none")


# get the 'weigthed' exposure terms
# individual rates for driver mutations:
DNMT3A_driver_list = split(DNMT3A_watson_drivers, DNMT3A_watson_drivers$aachange)
DNMT3A_driver_list = DNMT3A_driver_list[c("R882H", "R882C")]
name = names(DNMT3A_driver_list)[1]
list_drivers = list()

for (name in names(DNMT3A_driver_list)) {
  muts = DNMT3A_driver_list[[name]]
  list_drivers[[name]] = calc_exp_muts(expected_rates, muts, metadata, ratios, ncells)
}

expected_drivers = rbindlist(list_drivers, idcol = 'aachange')

expansion_df = data.frame(percent = c(18.7, 16.0, 15.8,15, 14.8, 14.1, 12.9, 12.3, 12.1, 12, 11.9, 11.2),
                            aachange = watson_variants) |>
  mutate(s = 1 + percent / 100,
         age_expansion = log(2000) / log(s))


figure_5 =

