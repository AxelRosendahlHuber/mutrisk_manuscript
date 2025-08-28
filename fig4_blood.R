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
  mutate(name = factor(name, levels = c("DNMT3A_drivers", "DNMT3A_R882H", "DNMT3A_watson_drivers")))

# correct the mutation loads for cihigh and cilow by the number of cells
mutation_list_ci = mutation_list |>
  mutate(cihigh = mle * 13,
         cilow = mle / 4)
figure_4a = ggplot(mutation_list_ci,
                   aes(x = age, y = mle, fill = category)) +
  geom_smooth(method = "lm", color = "black") +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh, color = category),
                width = 0,
                show.legend = FALSE) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_point(shape = 21, size = 2.4, color = "white", stroke = 0.3) +
  ggpubr::stat_cor() +
  scale_fill_manual(values = blood_colors) +
  scale_color_manual(values = blood_colors) +
  theme_cowplot() +
  labs(y = "number of cells/individual",
       x = "Age (years)",fill = NULL) +
  theme(legend.position = "none")
figure_4a

# for the fitness effect: we need to consider that we can only model mutations with a specific fitness effect

# perform two types of analyses:
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
figure_4b



ggplot(mutation_list_ci_prob |> filter(name == "DNMT3A_watson_drivers"),
                   aes(x = age, y = mle, fill = category)) +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh, color = category),
                width = 0,
                show.legend = FALSE) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_point(shape = 21, size = 2.4, color = "white", stroke = 0.3) +
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
         age_expansion = log(2000) / log(s) )

log(2000) / log(1.187)

# calculate the matching
expansion_group_DNMT3A = expected_drivers |>
  arrange(donor) |>
  left_join(expansion_df) |>
  group_by(donor) |>
  mutate(expansion_weight = mle / mean(mle),
         expansion = age_expansion * expansion_weight) |>
  summarize(across(c(mle, cilow, cihigh), sum),
            expansion = mean(expansion)) |>
  mutate(mle = get_prob_mutated_N(mle / 1e5, 1e5),
         cilow = get_prob_mutated_N(mle / 1e5, 2.5e4),
         cihigh = get_prob_mutated_N(mle / 1e5, 1.3e6))

# check if there is no cross of the lines (indicating an impossible situtation)
expansion_group_DNMT3A |>
  left_join(metadata) |>
  mutate(detection_age = age + expansion) |>
  ggplot(aes(x = detection_age, y = mle, ymax = cihigh, ymin = cilow)) +
  geom_pointrange() +
  geom_smooth(method = 'lm', fullrange = TRUE) +
  geom_pointrange(data = mutation_list_ci_prob |> filter(name == "DNMT3A_R882H"),
                  aes(x = age + 55.07050), color = "red")

expansion_group_DNMT3A |>
  left_join(metadata) |>
  mutate(detection_age = age + expansion) |>
  ggplot(aes(x = detection_age, y = mle, ymax = cihigh, ymin = cilow)) +
  geom_pointrange(color = "grey80") +
  geom_point(aes(x = age), ) +
  theme_cowplot() +
  geom_segment(aes(x = age ,xend = detection_age, y = mle), color = "grey80", linetype = "dashed")  +
  labs( y = "probability of 20 DNMT3A drivers", x = "Age in years", title = "DNMT3A R882H & R882C probability")