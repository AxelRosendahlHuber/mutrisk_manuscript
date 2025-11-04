# Figure_4_blood new
source("code/functions/analysis_variables.R")
library(patchwork)
library(cowplot)
library(geomtextpath)
library(ggpubr)
library(gridExtra)
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



plot_figures = function(driver_sites, y_label) {
  driver_rates = calc_exp_muts(expected_rates, driver_sites, metadata= metadata, ratios = ratios, ncells = ncells)

  high_estimate = ggplot(driver_rates,
               aes(x = age, y = mle*13)) +
    geom_point(color = blood_colors)  +
    labs(y = y_label, x = "Age (years)", title = "1.3 million HSCs") +
    theme_cowplot()


  mid_estimate = ggplot(driver_rates,
                aes(x = age, y = mle)) +
    geom_point(aes(y = mle), color = blood_colors)  +
    labs(y = y_label, x = "Age (years)", title = "100,000 HSCs") +
    theme_cowplot()

  return(list(high_estimate = high_estimate, mid_estimate = mid_estimate))
}


# DNMT3A driver hotspot
DNMT3A_R882H_hotspot = CH_bDM[aachange == "R882H" & gene_name == "DNMT3A", .N, c("gene_name", "mut_type", "aachange", "position", "driver")]
DNMT3A_R882H_hotspot_plots = plot_figures(DNMT3A_R882H_hotspot, "Cells with DNMT3A R882H mutation")

DNMT3A_drivers = CH_bDM[gene_name == "DNMT3A" & driver == TRUE , .N, c("gene_name", "mut_type", "aachange", "position", "driver")]
DNMT3A = plot_figures(DNMT3A_drivers, "Cells with DNMT3A driver mutation")

TET2_drivers = CH_bDM[gene_name == "TET2" & driver == TRUE , .N, c("gene_name", "mut_type", "aachange", "position", "driver")]
TET2 = plot_figures(TET2_drivers, "Cells with TET2 driver mutation")

TP53_drivers = CH_bDM[gene_name == "TP53" & driver == TRUE , .N, c("gene_name", "mut_type", "aachange", "position", "driver")]
TP53 = plot_figures(TP53_drivers, "Cells with TP53 driver mutation")

# make the general figure:
plots = c(DNMT3A, TET2, TP53)
wrap_plots(plots, byrow = FALSE)


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

# also needed would be the list of all the variants in the Watson analysis reported to be mutated
F5A = readRDS("processed_data/plots/F5A.rds")
F5A = prep_plot(F5A, label = "A", t = 8, r = 8, l = 8, b = 8)


CH_genes = c("DNMT3A", "TET2", "TP53")


# UKbiobank individual counts
UKB_age_frequencies = fread("raw_data/UKBiobank/UKB_age_frequencies_DNMT3A.tsv") |>
  select(Age, Individuals)
gene = "DNMT3A"

# start for loop:
UKB_plot_list = list()
for (i in 1:3) {

  gene = CH_genes[i]
  UKB_gene_muts = fread(paste0("raw_data/UKBiobank/UkBiobank_", gene, "_mut_age.csv"))
  UKB_gene_drivers = CH_bDM[gene_name %in% gene & driver == TRUE]

  age_samples = UKB_gene_muts |>
    filter(aa_change %in% UKB_gene_drivers$aachange) |>
    count(Age)

  UKB_age_incidence = left_join(UKB_age_frequencies, age_samples) |>
    filter(Individuals > 2000) |>
    mutate(relative_incidence = n / Individuals)

  # consider grouping the data as in colon in bins of 10
  plt = ggplot(UKB_age_incidence, aes(x = Age, y = relative_incidence)) +
    geom_point() +
    scale_y_continuous(limits = c(0, NA), labels = label_percent()) +
    theme_cowplot() +
    labs(x = "Age (years)", y = paste0("Incidence of CH in UKB with\n", gene, " driver mutation"), subtitle = gene)

  UKB_plot_list[[gene]] = plt
}
wrap_plots(UKB_plot_list)


# also consider making mirror plots for the individiual samples

# TODO: Save the figure still.
ggsave("manuscript/Figure_5/Figure_5.png", width = 15, height = 9)



