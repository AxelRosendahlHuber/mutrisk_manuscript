# Make mirror plots/ 'histograms' with all the individual mutations
# Figure 2: Add in mirror plots with TP53 number of mutations for all tissues
# Figure 1: Add in mirror plots with TP53 mutation probability for all tissues (check which one to use)
# Re-make mirror plots for figure 3
# Mirror plots also required for figure 4 and figure 5. Should be relatively easy to make all of them

# New script making mirror plots
library(data.table)
library(tidyverse)
library(cowplot)
library(MutationalPatterns)
library(patchwork)
source("code/functions/analysis_variables.R")
getwd()

# load data sources
GENIE_data = fread("processed_data/GENIE_17/GENIE_17_genie_tissue_type.txt.gz")
metadata_files = c("processed_data/blood/blood_metadata.tsv", "processed_data/colon/colon_metadata.tsv",
                   "processed_data/lung/lung_metadata.tsv")

names(metadata_files) = str_split_i(metadata_files, "\\/", 2)
metadata = lapply(metadata_files, \(x) fread(x)[,c("sampleID", "category", "age", "donor")]) |>
  rbindlist(idcol = "tissue")

# Load gene_of_interest boostdm
boostdm = fread("processed_data/boostdm/boostdm_genie_cosmic/pancancer_boostDM_intersect.txt.gz") |>
  mutate(driver = ifelse(driver == TRUE, "driver", "non-driver"))# change names for overview


# load the mutation rates
expected_rate_list = list()
ratio_list = list()
for (tissue in c("colon", "blood", "lung")) {
  expected_rate_list[[tissue]] = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
  ratio_list[[tissue]] = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))
}
expected_rates = rbindlist(expected_rate_list, idcol = "tissue", use.names = TRUE)
ratios = rbindlist(ratio_list, idcol = "tissue", use.names = TRUE)

# filters
gene_of_interest = "TP53"


merge_mutrisk_drivers = function(boostdm, ratios, gene_of_interest, tissue_select = "colon",
                                 category_select = "normal", cell_probabilities = FALSE,
                                 individual = FALSE, older_individuals = TRUE) {

  older_individuals = metadata |> filter(tissue == tissue_select,
                                         category %in% category_select,
                                         age > 30)
  ratio_gene_tissue = ratios |> filter(gene_name == gene_of_interest &
                                         category %in% category_select,
                                       tissue == tissue_select) |> pull(ratio)
  expected_rates_select = expected_rates[category %in% category_select &
                                           tissue == tissue_select, ] |>
    left_join(older_individuals, by = c("tissue", "sampleID", "category")) |>
    filter(donor %in% older_individuals$donor)

  ncells_select = tissue_ncells_ci[tissue_ncells_ci$tissue == tissue_select, "mid_estimate"]

  if (cell_probabilities == TRUE) {
    ncells_select = 1
  }

  # group by donor individually
  mutated_rates = expected_rates_select |>
    group_by(donor, mut_type, tissue) |>
    summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop") |>
    mutate(across(c(mle, cilow, cihigh), ~ . * ratio_gene_tissue * ncells_select))

  # modify for specific individual, or for the entiriety
  if (individual == FALSE) {
    print("taking mean of the mutation rates")

    mutated_rates_select = mutated_rates |>
      group_by(mut_type, tissue) |>
      summarize(mle = mean(mle))

    individuals = older_individuals |>
      select(donor, age) |> distinct()

    label = paste(tissue_select,  category_select, "all donors, average age:", format(mean(individuals$age), digits = 3))


  } else if (individual %in% unique(mutated_rates$donor)) {
    mutated_rates_select = mutated_rates |>
      filter(donor == individual) |>
      select(mut_type, tissue, mle)

    label = paste(category_select, "donor", individual, " age:", older_individuals[donor == individual] |> pull(age))

  }  else if (individual == "all") {
    mutated_rates_select = mutated_rates
    label = "no_label"
  }   else {print("parameter 'individual' must either be FALSE, or donor-id")}

  boostdM_goi = boostdm[gene_name == gene_of_interest, c("mut_type", "position", "driver")]
  expected_gene_muts = boostdM_goi |>
    full_join(mutated_rates_select, relationship = "many-to-many", by = "mut_type")

  return(list(expected_gene_muts = expected_gene_muts, label = label))
}


make_gene_barplot = function(boostdm, ratios, gene_of_interest, tissue_select = "colon",
                             category_select = "normal",
                             cell_probabilities = FALSE, individual = FALSE, older_individuals = TRUE) {

  mr_drivers = merge_mutrisk_drivers(boostdm, ratios, gene_of_interest, tissue_select, category_select, cell_probabilities,
                        individual, older_individuals)

  expected_gene_muts = mr_drivers$expected_gene_muts
  label = mr_drivers$label

  y_label = "number of cells with mutation"
  if (cell_probabilities == TRUE) {
    ncells_select = 1
    y_label = "probability of mutation"
  }

  if (max(expected_gene_muts$position, na.rm = TRUE) > 2000) {
    expected_gene_muts = expected_gene_muts |>
      mutate(position = (position - 1) %/% 5 + 1,
             position = position * 5) |>
      group_by(position, tissue, mut_type, driver)  |>
      summarise(mle = sum(mle, na.rm = TRUE), .groups = "drop")
    x_label = "AA position (5AA bins)"
  } else { x_label = "AA position"}


  expected_gene_muts_label = left_join(expected_gene_muts, mutrisk:::triplet_match_substmodel)

  # way to make the plot extend both upper and lower axes
  ggplot(expected_gene_muts_label,
         aes(x = position, y = mle, fill = type)) +
    geom_col() +
    scale_fill_manual(values = mutrisk::COLORS6) +
    theme_cowplot() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = x_label, y = y_label, title = gene_of_interest, subtitle = label, fill = NULL)
}


make_gene_barplot(boostdm, ratios, gene_of_interest = "KRAS", tissue_select = "colon")
make_gene_barplot(boostdm, ratios, gene_of_interest = "KMT2C")

barplot_lung = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53", tissue_select = "lung", category_select = "non-smoker") + ggtitle(NULL)
barplot_blood = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53", tissue_select = "blood") + ggtitle(NULL) + ylab(NULL)
barplot_colon = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53", tissue_select = "colon") +   ylab(NULL)
wrap_plots(barplot_colon, barplot_lung, barplot_blood, ncol = 1, guides = "collect")

prob_barplot_lung = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53", tissue_select = "lung", cell_probabilities = TRUE)
prob_barplot_blood = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53", tissue_select = "blood", cell_probabilities = TRUE)
prob_barplot_colon = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53", tissue_select = "colon", cell_probabilities = TRUE)
wrap_plots(prob_barplot_colon, prob_barplot_lung, prob_barplot_blood, ncol = 1, guides = "collect")


# APC
make_gene_barplot(boostdm, ratios, gene_of_interest = "APC", tissue_select = "colon", cell_probabilities = FALSE) +
  facet_grid(driver ~ .) +
  scale_y_continuous(breaks = extended_breaks(4)) +
  theme(legend.position = "inside", legend.position.inside = c(0.05, 0.8))


# make function, inputting the boostdm driver mutations, and the expected rates.
# this will become the plot
barplot_colon +
  facet_grid(driver ~ .) +
  scale_y_continuous(breaks = extended_breaks(4)) +
  theme(legend.position = "inside", legend.position.inside = c(0.05, 0.8))

# plotting function making barplots
make_summary_barplots = function(boostdm, ratios, gene_of_interest,
                                 cell_probabilities = FALSE, individual = FALSE, older_individuals = TRUE) {

  tissue_samples = unlist(tissue_colors) |>
    names()

  # first thing to
  output_list = list()
  for (i in tissue_samples){
    tissue_select = str_split_1(i, "\\.")[1]
    category_select = str_split_1(i, "\\.")[2]
    # merge the samples in a for loop
    mr_rates = merge_mutrisk_drivers(boostdm, ratios, gene_of_interest, tissue_select = tissue_select,
                                            category_select = category_select, cell_probabilities = FALSE, individual = FALSE,
                                            older_individuals = TRUE)[[1]]
    output_list[[i]] = mr_rates |>
      group_by(driver, tissue) |>
      summarize(mle = sum(mle)) |>
      mutate(category = category_select)
  }

  output_rates = rbindlist(output_list, idcol = "tissue.color") |>
    mutate(driver = factor(driver, levels = c("non-driver", "driver")),
           tissue_category = gsub("\\.","_", tissue.color),
           category = factor(category, levels = rev(unique(color_df$category)))) |>
    mutate(label = paste0(driver, ": ", prettyNum(mle, big.mark = ",", digits = 2, scientific = FALSE)),
           label_x = ifelse(driver == "driver", 0, lag(mle)))

  plot_list = list()

  for (tissue_select in c("colon", "lung", "blood")) {
    tissue_rates = output_rates |> filter(tissue %in% tissue_select)
    plot_list[[tissue_select]] = ggplot(tissue_rates |> filter(tissue.color != "blood.chemotherapy"),
         aes(y = category, x = mle, fill = tissue.color)) +
      geom_col(aes(alpha = driver), color = "black") +
      geom_text(aes(label = label, x = label_x, color = driver), hjust = -0.2) +
      scale_fill_manual(values = unlist(tissue_colors)) +
      facet_grid(tissue ~., scales = "free", space = "free_y") +
      theme_cowplot() +
      scale_alpha_manual(values = c(0.5, 1)) +
      scale_color_manual(values = c("black", "white")) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.1)), labels = label_comma(), breaks = extended_breaks(n = 4)) +
      ggtitle(paste(gene_of_interest, tissue_select)) +
      labs(y = NULL, x = "Number of mutated cells") +
      theme(legend.position = "none", strip.text = element_blank())
  }

  wrap_plots(plot_list, ncol = 1, heights = c(4, 3, 1.5))
}

make_summary_barplots(boostdm, ratios, "TP53",
                      cell_probabilities = FALSE, individual = FALSE,
                      older_individuals = TRUE)


make_summary_barplots(boostdm, ratios, "APC",
                      cell_probabilities = FALSE, individual = FALSE,
                      older_individuals = TRUE)

make_summary_barplots(boostdm, ratios, "KRAS",
                      cell_probabilities = FALSE, individual = FALSE,
                      older_individuals = TRUE)

# use boostdm ch for DNMT3A
boostDM_ch = fread("processed_data/boostdm/boostdm_genie_cosmic/CH_boostDM_cancer.txt.gz") |>
  mutate(driver = ifelse(driver == TRUE, "driver", "non-driver"))# change names for overview
make_summary_barplots(boostdm = boostDM_ch, ratios = ratios, gene_of_interest = "DNMT3A",
                      cell_probabilities = FALSE, individual = FALSE,
                      older_individuals = TRUE)


# make dotplots with the individual variation:
all_rates = merge_mutrisk_drivers(boostdm, ratios, gene_of_interest = "TP53", tissue_select = "colon", individual = "all")[[1]] |>
  group_by(donor, driver) |>
  summarize(mle = sum(mle))

ggplot(all_rates, aes(x = driver, y = mle, alpha = driver)) +
  geom_boxplot(fill = colon_colors[1]) +
  ggbeeswarm::geom_beeswarm(cex = 1.7) +
  scale_alpha_manual(values = c(1, 0.5)) +
  scale_y_continuous(labels = label_comma(), limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  theme_cowplot() +
  labs(y = "number of cells with TP53 mutation", x = NULL) +
  theme(legend.position = "none")

# boxplots for all?
boxplot_list = list()
for (i in 1:nrow(color_df)) {

  category_select = color_df$category[i]
  tissue_select = color_df$tissue[i]

  # make dotplots with the individual variation:
  boxplot_list[[i]] = merge_mutrisk_drivers(boostdm, ratios, gene_of_interest = "TP53", tissue_select = tissue_select, category_select = category_select,
                                    individual = "all")[[1]] |>
    group_by(donor, driver) |>
    summarize(mle = sum(mle)) |>
    mutate(tissue = tissue_select,
           category = category_select)
}

boxplot_df = rbindlist(boxplot_list) |>
  mutate(tissue_category = paste0(tissue, "_", category))

ggplot(boxplot_df |> filter(category != "chemotherapy"), aes(x = category, y = mle, alpha = driver, fill = tissue_category)) +
  geom_boxplot(position = "dodge") +
  scale_fill_manual(values = tissue_category_colors) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_alpha_manual(values = c(1, 0.5)) +
  facet_wrap(. ~ tissue, scale = "free", space = "free_x") +
  theme_cowplot()


# UKBiobank DNMT3A mutations:
UKB_DNMT3A_muts = fread("raw_data/UKBiobank/UkBiobank_DNMT3A_mut_age.csv")
UKB_DNMT3A_counts = UKB_DNMT3A_muts[, .N, by = c("aa_change", "REF", "ALT")]  |>
  mutate(position = parse_number(aa_change),
         type = paste0(REF, ">", ALT),
         type = case_match(type, .default = type,
                           "A>T" ~ "T>A", "A>G" ~ "T>C", "A>C" ~ "T>G",
                           "G>T" ~ "C>A", "G>A" ~ "C>T", "G>C" ~ "C>G"))

DNMT3A_plot = UKB_DNMT3A_counts |>
  ggplot(aes(x = position, y = N, fill = type)) +
  geom_col() +
  scale_fill_manual(values = COLORS6) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x  = NULL , y = "Number of CH mutations\nobserved in the UKBiobank cohort", title = "DNMT3A", fill = NULL)

barplot_blood = make_gene_barplot(boostDM_ch, ratios, gene_of_interest = "DNMT3A", tissue_select = "blood") +
  scale_y_continuous(trans = "reverse") +
  labs(title = NULL, subtitle = NULL)

DNMT3A_plot / barplot_blood + plot_layout(guides = "collect")

