# Script to produce figure 3, and the barplots indicating mutation accumulatoin in figure 4
library(data.table)
library(tidyverse)
library(cowplot)
library(MutationalPatterns)
library(patchwork)
library(ggh4x)
source("code/0_functions/analysis_variables.R")
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
  merge_mutrisk_drivers = function(boostdm, ratios, gene_of_interest, tissue_select = "colon", tissue_name,
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

      label = paste(tissue_name, "- average age:", format(mean(individuals$age), digits = 3))


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


make_gene_barplot = function(boostdm, ratios, gene_of_interest,
                             tissue_select = "colon", tissue_name = NULL,
                             category_select = "normal",
                             cell_probabilities = FALSE, individual = FALSE, older_individuals = TRUE,
                             lollipop_dots = FALSE) {

  if (is.null(tissue_name)) {tissue_name = tissue_select}

  mr_drivers = merge_mutrisk_drivers(boostdm, ratios, gene_of_interest, tissue_select, tissue_name, category_select, cell_probabilities,
                        individual, older_individuals)

  expected_gene_muts = mr_drivers$expected_gene_muts
  label = mr_drivers$label

  y_label = "Number of cells with mutation"
  if (cell_probabilities == TRUE) {
    ncells_select = 1
    y_label = "Probability of mutation\n per cell(x10⁻⁶)"
  }

  if (max(expected_gene_muts$position, na.rm = TRUE) > Inf) { # for now set the level to Inf to allow for large genes
    expected_gene_muts = expected_gene_muts |>
      mutate(position = (position - 1) %/% 5 + 1,
             position = position * 5) |>
      group_by(position, tissue, mut_type, driver)  |>
      summarise(mle = sum(mle, na.rm = TRUE), .groups = "drop")
    x_label = "AA position (5AA bins)"
  } else { x_label = "AA position"}

  expected_gene_muts_label = left_join(expected_gene_muts, mutrisk:::triplet_match_substmodel)

  # way to make the plot extend both upper and lower axes
  pl = ggplot(expected_gene_muts_label,
         aes(x = position, y = mle)) +
    geom_col(aes(fill = type)) +
    scale_fill_manual(values = mutrisk::COLORS6) +
    theme_cowplot() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = label_comma()) +
    scale_x_continuous(expand = expansion(mult = c(0.01,0.01))) +
    labs(x = x_label, y = y_label, title = gene_of_interest, subtitle = label, fill = NULL)

  if (cell_probabilities == TRUE) {
    pl = pl + scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = function(x) x * 1e6)
  }

  pl
}

# Make a barplot showing the probabilitites for TP53 (poster usage)
prob_barplot_lung = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53", tissue_select = "lung", category_select = "non-smoker",
                                      individual = "PD34215",cell_probabilities = TRUE) + labs(title = "TP53", subtitle = NULL, y = NULL)
prob_barplot_blood = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53", tissue_select = "blood",
                                       individual = "KX008", cell_probabilities = TRUE) + labs(title = "TP53", subtitle = NULL, y = NULL)
prob_barplot_colon = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53", tissue_select = "colon",
                                       individual = "O340", cell_probabilities = TRUE) + labs(title = "TP53", subtitle = NULL)

# Make a barplot indicating the number of mutations across TP53 across the three tissues (colon, lung, blood)
barplot_colon = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53", tissue_select = "colon",
                                  tissue_name = "Colon", cell_probabilities = FALSE) + labs(y = NULL)
barplot_lung = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53", tissue_select = "lung",
                                 tissue_name = "Lung", category_select = "non-smoker",
                                 cell_probabilities = FALSE)
barplot_blood = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53", tissue_select = "blood",
                                  tissue_name = "Blood", cell_probabilities = FALSE) + labs(y = NULL)
F3A = wrap_plots(barplot_colon, barplot_lung, barplot_blood, ncol = 3, guides = "collect") &
  labs(title = "TP53 mutations across normal tissues")
saveRDS(F3A, "manuscript/figure_panels/figure_3/figure_3A.rds")

# Supplementary Figure 2 - TP53 for all tissues
tissue_categories = ratios |> select(tissue, category) |> distinct()
tissue_categories = tissue_categories[-5]
plot_list = list()
for (i in 1:nrow(tissue_categories)) {

  tissue_name = paste(as.character(tissue_categories[i]), collapse = "\n")
  plot_list[[i]] = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53",
                                     tissue_select = tissue_categories$tissue[i],
                                     category_select = tissue_categories$category[i],
                                    tissue_name = tissue_name, cell_probabilities = FALSE) +
    labs(title = NULL, y = "Number of cells \nwith mutation") + theme_classic()
}

plot_list = c(plot_list)
figure_S2 = wrap_plots(plot_list, nrow = 4) + plot_layout(guides = "collect") +
  plot_annotation(title = 'TP53: Expected number of mutated cells')
ggsave("manuscript/Supplementary_Figures/Figure_S2/Figure_S2.png", figure_S2, width = 14, height = 12)
ggsave("manuscript/Supplementary_Figures/Figure_S2/Figure_S2.pdf", figure_S2, width = 14, height = 12)

# APC colon barplot
APC_colon_normal = make_gene_barplot(boostdm, ratios, gene_of_interest = "APC",
                                     tissue_select = "colon", category_select = "normal", cell_probabilities = FALSE) +
  ggh4x::facet_grid2(driver ~ ., strip = strip_themed(background_y = elem_list_rect(fill = c("#C03830", "#707071")),
                                                      text_y = elem_list_text(colour = c("white"), face = "bold")), axes = "all",
                     remove_labels = "x")



colon_normal = make_gene_barplot(boostdm, ratios, gene_of_interest = "TP53",
                                 tissue_select = "colon", category_select = "normal", cell_probabilities = FALSE) +
  ggh4x::facet_grid2(driver ~ ., strip = strip_themed(background_y = elem_list_rect(fill = c("#C03830", "#707071")),
                                                      text_y = elem_list_text(colour = c("white"), face = "bold")), axes = "all",
                     remove_labels = "x")
ggsave("plots/colon/TP53_driver_non-driver.png", colon_normal, width = 10, height = 4.5, bg = "white")

# APC colon barplot
APC_colon_normal = make_gene_barplot(boostdm, ratios, gene_of_interest = "APC",
                                 tissue_select = "colon", category_select = "normal", cell_probabilities = FALSE) +
  ggh4x::facet_grid2(driver ~ ., strip = strip_themed(background_y = elem_list_rect(fill = c("#C03830", "#707071")),
                                                      text_y = elem_list_text(colour = c("white"), face = "bold")), axes = "all",
                     remove_labels = "x")

# Add dots to the colon barplot to make the individual bars more visible
df_dots = APC_colon_normal@data |>
  group_by(position, driver) |>
  summarize(mle = sum(mle))

APC_colon_normal = APC_colon_normal + geom_point(data = df_dots, aes(x = position, y = mle), size = 1.5) +
  scale_y_continuous(limits = c(NA, 2250), expand = expansion(mult = c(0, 0.1)))
ggsave("plots/colon/APC_driver_non-driver.png", APC_colon_normal, width = 12, height = 5, bg = "white")


# Figure 4A
F4A1 = make_gene_barplot(boostdm, ratios, gene_of_interest = "APC", tissue_select = "colon", cell_probabilities = FALSE) +
  scale_y_continuous(breaks = extended_breaks(4), expand = expansion(mult = c(0, 0.1))) +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = "none",
        legend.text = element_text(size = rel(0.8)),
        legend.title = element_text(size = rel(0.8)),
        legend.key.size = unit(0.8, "lines"), legend.background = element_blank())

F4A2 = make_gene_barplot(boostdm, ratios, gene_of_interest = "KRAS", tissue_select = "colon", cell_probabilities = FALSE) +
  scale_y_continuous(breaks = extended_breaks(4), expand = expansion(mult = c(0, 0.1))) +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 1),
        legend.text = element_text(size = rel(0.8)),
        legend.title = element_text(size = rel(0.8)),
        legend.key.size = unit(0.8, "lines"), legend.background = element_blank()) +
  scale_x_continuous(expand = c(0,0))

saveRDS(list(F4A1, F4A2), "manuscript/figure_panels/figure_4/figures_AB.rds")


# plot the number of mutations for TP53 as individual points
dotplot_list = list()
for (i in 1:nrow(color_df)) {

  category_select = color_df$category[i]
  tissue_select = color_df$tissue[i]

  # make dotplots with the individual variation:
  dotplot_list[[i]] = merge_mutrisk_drivers(boostdm, ratios, gene_of_interest = "TP53",
                                            tissue_select = tissue_select, category_select = category_select,
                                    individual = "all")[[1]] |>
    group_by(donor, driver) |>
    summarize(across(c(mle, cilow, cihigh), sum)) |>
    mutate(tissue = tissue_select, category = category_select)
}

dotplot_df = rbindlist(dotplot_list) |>
  mutate(tissue = factor(tissue, levels = c("colon", "lung", "blood")),
                         tissue_category = paste0(tissue, "_", category))

df_total_muts = dotplot_df |>
  filter(category != "chemotherapy") |>
  left_join(metadata |> select(-sampleID) |> distinct()) |>
  mutate(category = factor(category, levels = unique(color_df$category)))

tissue_plots = tissue_plots_raw =  list()
tissue_order = c("colon", "lung", "blood")
for (i in 1:3) {

  select_tissue = tissue_order[[i]]
  df_tissue = df_total_muts |>
    filter(tissue %in% select_tissue) |>
    filter(driver == "driver")

  plt = ggplot(df_tissue, aes(x = age, y = mle, color = tissue_category)) +
    geom_pointrange(aes(ymin = cilow, ymax = cihigh)) +
    facet_nested(. ~ tissue + category, axes = "y", remove_labels = "y") +
    scale_color_manual(values = tissue_category_colors) +
    scale_y_continuous(labels = scales::label_comma(), limits = c(0,NA)) +
    theme_cowplot() +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          legend.position = "none", ggh4x.facet.nestline = element_line()) +
    labs(x = "Age (years)", y = "Number of cells with\nTP53 driver mutation")
  tissue_plots_raw[[select_tissue]] = plt
  tissue_plots[[select_tissue]] = prep_plot(plt, LETTERS[i+2])
}

tissue_plots_raw[[2]] = tissue_plots_raw[[2]] + labs(y = NULL)
tissue_plots_raw[[3]] = tissue_plots_raw[[3]] + labs(y = NULL)

saveRDS(tissue_plots_raw, "manuscript/figure_panels/figure_3/figure_3CDE.rds")

# Figures for poster (can be removed if needed)
# figure exploration for poster:
poster = wrap_plots(tissue_plots_raw, nrow = 1, widths = c(4, 3, 1.2))
F3C = wrap_plots(tissue_plots, nrow = 1, widths = c(4, 3, 1.2))
F3C

# alternative poster figure [two rows of figures ]
F3C_bottom = wrap_plots(tissue_plots[-1], widths = c(2.8, 1))
F3C = tissue_plots[[1]] / F3C_bottom
F3C


