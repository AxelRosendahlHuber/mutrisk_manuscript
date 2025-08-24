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


# get the 'weigthed' exposure terms
# individual rates for driver mutations:
DNMT3A_driver_list = split(DNMT3A_watson_drivers, DNMT3A_watson_drivers$aachange)
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
  geom_pointrange() +
  geom_point(aes(x = age), color = "grey80") +
  theme_cowplot() +
  geom_segment(aes(x = age ,xend = detection_age, y = mle), color = "grey80", linetype = "dashed")  +
  labs( y = "probability of one of 20 DNMT3A driver mutations", x = "Age in years")

# # Read the UKbiobank data:
# # make figures for CH, DNMT3A and DNMT3A R882H mutation
# DNMT3A_age_frequencies = fread("raw_data/UKBiobank/UKB_age_frequencies_DNMT3A.tsv")
# CH_age_frequencies = fread("raw_data/UKBiobank/UKB_age_frequencies.tsv")
#
# CH_frequencies = CH_age_frequencies |>
#   mutate(CH = rowSums(select(CH_age_frequencies, c(-Age, -Individuals)))) |>
#   mutate(fraction = CH / Individuals) |>
#   filter(Individuals > 100)
#
# # compare the CH frequencies over time with the predictions:
# ch_genes = intersect(CH_bDM$gene_name, colnames(CH_age_frequencies))
#
# CH_driver_counts = CH_bDM[driver == TRUE , .N, by = c("gene_name", "mut_type", "driver")]
# CH_driver_count_list = split(CH_driver_counts, CH_driver_counts$gene_name)
# CH_driver_muts = lapply(CH_driver_count_list, \(x) calc_exp_muts(expected_rates, x, metadata, ratios, ncells)) |>
#   rbindlist(idcol = "gene_name")
#
# CH_rates = CH_frequencies |>
#   select(-fraction) |>
#   pivot_longer(-c(Age, Individuals), names_to = "gene_name", values_to = "n_CH") |>
#   mutate(mle = n_CH / Individuals, age = Age) |>
#   filter(gene_name %in% ch_genes)
#
# # plot for all genes the induction:
# plot_cr_driver_muts = CH_driver_muts |>
#   filter(gene_name %in% ch_genes) |>
#   ggplot(aes(x = age, y = mle)) +
#   geom_point() +
#   geom_point(data = CH_rates) +
#   facet_wrap(gene_name ~ . , scales = "free_y") +
#   theme_cowplot() + panel_border() +
#   labs(y = "incidence driver mutations/CH", y = "Age (years)", title = "linear y-axis")
# plot_cr_driver_muts
#
# plot_cr_driver_muts +
#   scale_y_log10() + labs(title = "log10 y-axis")
#
# plot_cr_driver_muts = CH_driver_muts |>
#   filter(gene_name == "DNMT3A") |>
#   mutate(cilow = mle / 4, cihigh = mle * 13) |>  # simulate having 25.00 or 1.3M cells
#   ggplot(aes(x = age, y = mle)) +
#   geom_smooth(method = "lm", se = FALSE) +
#     geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
#   geom_point(data = CH_rates |> filter(gene_name == "DNMT3A")) +
#   facet_wrap(gene_name ~ . , scales = "free_y") +
#   theme_cowplot() + panel_border() +
#   labs(y = "incidence driver mutations/CH", y = "Age (years)",
#        subtitle = "Confidence interval: Low = 25,000 HSCs, high 1.3 million HSCs")
# plot_cr_driver_muts
#
#
# # extend the expansion timing according to the estimates proposed by Watson et al.,
# exp_R882C = log(2000) / log(1.187)
# exp_R882C
# exp_R882H = log(2000) / log(1.148)
# exp_R882H
#
# plot_cr_driver_muts = CH_driver_muts |>
#   filter(gene_name == "DNMT3A") |>
#   mutate(cilow = mle / 4, cihigh = mle * 13,
#          age = age + 46) |>  # simulate having 25.00 or 1.3M cells
#   ggplot(aes(x = age, y = mle)) +
#   geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
#   geom_smooth(method = "lm", se = FALSE) +
#   geom_point(data = CH_rates |> filter(gene_name == "DNMT3A")) +
#   facet_wrap(gene_name ~ . , scales = "free_y") +
#   theme_cowplot() + panel_border() +
#   labs(y = "incidence driver mutations/CH", y = "Age (years)",
#        title = "Allowing for 46 years of expansion time",
#        subtitle = "Confidence interval: Low = 25,000 HSCs, high 1.3 million HSCs")
# plot_cr_driver_muts
#
#
# plot_cr_driver_muts = CH_driver_muts |>
#   filter(gene_name == "DNMT3A") |>
#   mutate(cilow = get_prob_mutated(mle / 4 / ncells, ncells),
#          mle = get_prob_mutated(mle / ncells, ncells),
#          cihigh = get_prob_mutated(mle * 13 / ncells, ncells),
#                  age = age ) |>  # simulate having 25.00 or 1.3M cells
#   ggplot(aes(x = age, y = mle)) +
#   geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
#   geom_point(data = CH_rates |> filter(gene_name == "DNMT3A")) +
#   facet_wrap(gene_name ~ . , scales = "free_y") +
#   theme_cowplot() + panel_border() +
#   labs(y = "incidence driver mutations/CH", y = "Age (years)",
#        title = "Probability of mut/CH incidence",
#        subtitle = "Confidence interval: Low = 25,000 HSCs, high 1.3 million HSCs")
# plot_cr_driver_muts
#
#
# plot_cr_driver_muts = CH_driver_muts |>
#   filter(gene_name == "DNMT3A") |>
#   mutate(cilow = mle / 4, cihigh = mle * 1,
#          age = age + 46) |>  # simulate having 25.00 or 1.3M cells
#   ggplot(aes(x = age, y = mle)) +
#   geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
#   geom_point(data = CH_rates |> filter(gene_name == "DNMT3A")) +
#   facet_wrap(gene_name ~ . , scales = "free_y") +
#   theme_cowplot() + panel_border() +
#   labs(y = "incidence driver mutations/CH", y = "Age (years)",
#        title = "Probability of mut: Allowing for 46 years of expansion time",
#        subtitle = "Confidence interval: Low = 25,000 HSCs, high 1.3 million HSCs")
# plot_cr_driver_muts
#
# # plot mutation rates in all genes
# site_freqs_ch_drivers = CH_bDM[driver == TRUE, .N, by = c("mut_type", "gene_name")]
# mrate_CH_drivers = get_gene_rate(exp_rates = expected_rates, metadata = metadata,
#                                  site_freqs = site_freqs_ch_drivers,  ratios = ratios,  ncells = ncells)
#
#
#
# # Get the last point of each group
# df_last <- mrate_CH_drivers |>
#   group_by(gene_name) |>
#   filter(mle == max(mle)) |>
#   arrange(mle) |>
#   tail(15)
#
# CH_driver_rates = ggplot(mrate_CH_drivers |> filter(gene_name %in% df_last$gene_name),
#                          aes(x = age, y = mle, color = gene_name)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE) +
#   ggrepel::geom_text_repel(data = df_last, aes(label = gene_name), nudge_x = 7,
#                            size = 4, na.rm = TRUE, max.overlaps = Inf) +
#   expand_limits(x = max(df_last$age) + 10) + # Add space for labels
#   ggsci::scale_color_igv() +
#   theme_cowplot() +
#   theme(legend.position = "none") +
#   labs(y = "number of CH driver\nmutations", x = "Age (years)")
# ggsave("plots/blood/comparison_muts_CH_rates.png", CH_driver_rates,  width = 6, height = 5, bg = "white")
#
# #  TODO Update the blood script with the new estimates for high and low numbers of mutation rate
# # min number of stem cells in the blood
# min_ncells = 2.5e4
# max_ncells = 1.3e6
#
# # DNMT3A mutations all
# site_freqs_ch_drivers = CH_bDM[driver == TRUE & gene_name == "DNMT3A", .N, by = c("mut_type", "gene_name")]
# estimates_area = get_mut_est_conf(site_freqs = site_freqs_ch_drivers, exp_rates = expected_rates,
#                                   metadata = metadata, ratios = ratios, ncells = ncells,
#                                   min_ncells = min_ncells, max_ncells = max_ncells)
#
# plot_muts_area = function(estimates) {
#   estimates |>
#     filter(gene_name == "DNMT3A") |>
#     ggplot(aes(x = age, y = ncells)) +
#     geom_ribbon(aes(ymin = model_min, ymax = model_max ), alpha = 0.1) +
#     geom_line(aes(y = model_ncells), linewidth = 1) +
#     geom_line(aes(y = model_min), color = "grey30", linewidth = 1, linetype = "dashed") +
#     geom_line(aes(y = model_max), color = "grey30", linewidth = 1, linetype = "dashed") +
#     theme_cowplot()  +
#     labs(y = "mutated HSCs with DNMT3A driver", x = "Age (years)")
# }
#
#
# label_min = paste0("25,000 HSCs\nEst. muts 80y: ", round(estimates_area[["model_min"]][8], 2))
# label_mid = paste0("100,000 HSCs\nEst. muts 80y: ", round(estimates_area[["model_ncells"]][8], 2))
# label_max = paste0("1.3M HSCs\nEst. muts: 80y: ", round(estimates_area[["model_max"]][8], 2))
#
# # plot mutation rate estimates for different initial estimates - plot in the confidence interval
# plot_confi = plot_muts_area(estimates_area)  +
#   labs(y = "Number of cells with DNMT3A\ driver mutation") +
#   annotate("text", x = 80, y = estimates_area[["model_ncells"]][8], label = label_mid, hjust = 0, vjust = 0, size = 3.5) +
#   annotate("text", x = 80, y = estimates_area[["model_min"]][8], label = label_min, hjust = 0, vjust = 0.6, color = "grey30", size = 3.5) +
#   annotate("text", x = 80, y = estimates_area[["model_max"]][8], label = label_max, hjust = 0, color = "grey30", size = 3.5) +
#   theme(plot.margin = margin(5,32,5,5, unit = "mm")) +
#   coord_cartesian(clip = "off")
# plot_confi
#
# site_freqs_ch_drivers = CH_bDM[driver == TRUE & gene_name == "DNMT3A" & aachange == "R882H",.N,by = c("mut_type", "gene_name")]
# estimates_area = get_mut_est_conf(site_freqs = site_freqs_ch_drivers, exp_rates = expected_rates,
#                                   metadata = metadata, ratios = ratios, ncells = ncells,
#                                   min_ncells = min_ncells,max_ncells = max_ncells)
#
# # check the ages needed for
# delay100k <- log(2000) / log(1.18)
# delay25k <- log(500) / log(1.18) # delay is much closer

