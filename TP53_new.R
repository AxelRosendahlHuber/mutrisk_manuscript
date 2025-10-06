# Get the mutation rates for driver genes. This for Figure 3, 4 and possibly 5.
library(data.table)
library(tidyverse)
library(cowplot)
library(MutationalPatterns)
source("code/functions/analysis_variables.R")
getwd()

# Select all gene_of_interest mutations in gene_of_interest
gene_oi = "TP53" # gene of interest TP53 (can be changed into different genes, for instance = "SF3B1", "SRSF2"

# Load GENIE mutation data
genie_data = fread("raw_data/GENIE_17/data_mutations_extended.txt")
nucs = c("A", "C", "G", "T")
genie_data = genie_data |>
  filter(Reference_Allele %in% nucs & Tumor_Seq_Allele2 %in% nucs) |>
  mutate(type = paste0(Reference_Allele, ">", Tumor_Seq_Allele2)) |>
  mutate(type = case_when(type == "A>T" ~ "T>A",
                          type == "A>C" ~ "T>G",
                          type == "A>G" ~ "T>C",
                          type == "G>C" ~ "C>G",
                          type == "G>A" ~ "C>T",
                          type == "G>T" ~ "C>A", .default = type))

genie_metadata = fread("raw_data/GENIE_17/GENIE_release_17.0/data_clinical_sample.txt")
genie_samples = genie_metadata[`Cancer Type` %in% c("Colorectal Cancer", "Non-Small Cell Lung Cancer", "Leukemia"),] |>
  select(`Sample Identifier`, `Age at Which Sequencing was Reported (Years)`, `Cancer Type`) |>
  `colnames<-`(c("Tumor_Sample_Barcode", "age", "tissue_type"))
genie_tissue_type_all = inner_join(genie_samples, genie_data)
fwrite(genie_tissue_type_all, "processed_data/GENIE_17/GENIE_17_genie_tissue_type.txt.gz")

# Count number of tumors in the GENIE data
sample_counts = genie_samples |>
  group_by(tissue_type) |>
  summarize(n_tumors = dplyr::n()) |>
  mutate(tissue_type = factor(tissue_type, levels = c("Leukemia", "Colorectal Cancer", "Non-Small Cell Lung Cancer", "Melanoma")))

plot_gene_oi = function(genie_data, gene_oi) {

  print(gene_oi)
  genie_gene_oi = genie_tissue_type_all[Hugo_Symbol == gene_oi, ]

  gene_of_interest_counts = genie_gene_oi |>
    group_by(tissue_type) |>
    summarize(n_gene_of_interest_mut = dplyr::n())

  labels = merge(sample_counts, gene_of_interest_counts) |>
    mutate(n_tumors = format(n_tumors, big.mark = ","),
           n_gene_of_interest_mut  = format(n_gene_of_interest_mut, big.mark = ","),
           label = paste0(tissue_type, "\nnumber tumors: ", n_tumors, "\nnumber ",  gene_oi, " muts: ", n_gene_of_interest_mut))

  genie_position_counts = genie_gene_oi |>
    group_by(Protein_position, type, tissue_type) |>
    mutate(tissue_type = factor(tissue_type, levels = levels(sample_counts$tissue_type))) |>
    count()

  genie_total_plot = ggplot(genie_position_counts, aes(x = Protein_position, y = n, fill = type)) +
    geom_col() +
    ggpp::geom_text_npc(data = labels, aes(label = label, npcx = 0.01, npcy = 0.9), size = 3) +
    facet_grid(tissue_type ~ . , scales = "free_y") +
    scale_fill_manual(values = COLORS6) +
    theme_cowplot() +
    panel_border() +
    theme(legend.position = "none", strip.text = element_blank(), axis.text.y = element_text(size = 9)) +
    labs(y = "number of observed mutations", subtitle = paste(gene_oi, "mutations across tissue_types - GENIE data"), x = "AA position") +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)), labels = label_comma())
  ggsave(paste0("plots/genie_genes/GENIE_", gene_oi, "_tissue_types.png"), genie_total_plot, width = 10, height = 6, bg = "white")
  genie_total_plot
  # get the consequence
  genie_gene_oi |> count(tissue_type, Consequence)
}

most_mutated = genie_tissue_type_all |> count(Hugo_Symbol) |> arrange(desc(n)) |>
  head(5) |> pull(Hugo_Symbol)

gene_oi = "TP53"
sapply(most_mutated, plot_gene_oi, genie_data = genie_data)
plot_gene_oi(genie_data, "RB1")

######
## MutRisk estimates for neutral mutation accumulation in healthy tissues
######
# load metadata files for the different samples
metadata_files = c("processed_data/blood/blood_metadata.tsv", "processed_data/colon/colon_metadata.tsv",
                   "processed_data/lung/lung_metadata.tsv")

names(metadata_files) = str_split_i(metadata_files, "\\/", 2)
metadata = lapply(metadata_files, \(x) fread(x)[,c("sampleID", "category", "age", "donor")]) |>
  rbindlist(idcol = "tissue")

# load the mutation rates
expected_rate_list = list()
ratio_list = list()
for (tissue in c("colon", "blood", "lung")) {
  expected_rate_list[[tissue]] = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
  ratio_list[[tissue]] = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))
}
expected_rates = rbindlist(expected_rate_list, idcol = "tissue", use.names = TRUE)
ratios = rbindlist(ratio_list, idcol = "tissue", use.names = TRUE)

# Mean mutation rate for each trinucleotide
# Check if the blood rates actually make sense - this seems too low - this must be because of the cord blood donors
mean_rates = mutation_rates |>
  group_by(tissue_type, mut_type) |>
  summarize(mle = mean(mle),
            mean_age = mean(age)) |>
  mutate(tissue_type_age = paste0(tissue_type, "\n(", format(mean_age, digits = 3, nsmall = 1 ), ") years"))

# get rates for a plot
# add the mean rates for the cohort -> make sure that we weigh each of the donors equally..
mean_rates_table = mutation_rates |>
  group_by(tissue_type, mut_type, donor, age) |>
  summarize(mle = mean(mle)) |>
  group_by(tissue_type) |>
  summarize(`Mean Age` = round(mean(age), 1)) |>
  mutate(Tissue = gsub("\\\n.*", "",  tissue_type)) |>
  select(Tissue, `Mean Age`) |>
  distinct()
mean_rates_table$`Cell Type` = c("Hematopoietic stem cells", "Crypt stem cells", "Lung basal cells", "Melanocytes")
table_1 = wrap_table(mean_rates_table, space = unit(0.1, "npc"))

# determine the ratio of exonic to genomic mutations
# load the number of WGS genomes
select_categories = c("smoker", "normal")
wgs_muts = list.files("processed_data/", pattern = "cell_muts.tsv", full.names = TRUE, recursive = TRUE)
wgs_muts = wgs_muts[c(2,3,5,6)]
names(wgs_muts) = names(files)
muts = lapply(wgs_muts, \(x) fread(x)[,c("sampleID", "chr", "pos", "ref", "alt")]) |>
  rbindlist(idcol = "tissue_type")

wgs_muts = muts |>
  group_by(tissue_type) |>
  summarize(n_wgs = dplyr::n())

# load the number of exonic mutations in the cohort
annot_files = c("processed_data/blood/blood_normal_dnds.rds", "processed_data/colon/colon_normal_dnds.rds",
                "processed_data/lung/lung_smoker_dnds.rds", "processed_data/skin/skin_melanocyte_dnds.rds")
names(annot_files) = names(files)[1:3] # exclude skin - as this is whole-exome sequencing we cannot calculate a WGS/WES ratio
readRDS(annot_files[1])$annotmuts |>
  filter(gene != "intronic")
exonicmuts = lapply(annot_files, \(x) readRDS(x)[["annotmuts"]]) |>
  rbindlist(idcol = "tissue_type", fill = TRUE) |>
  filter(gene != "intronic") |>
  group_by(tissue_type) |> summarize(n_exonic = dplyr::n())

merged_muts = merge(wgs_muts, exonicmuts) |>
  mutate(fraction_exome_genome = n_exonic / n_wgs,
         label = paste0(format(fraction_exome_genome * 100, big.mark = ",", digits = 2), "%"))

# plot which indicates the fraction of mutations in the genome and exome (potentially supplementary figure)
ratio_exome_genome = ggplot(merged_muts, aes(x = tissue_type, y = fraction_exome_genome)) +
  geom_col() +
  geom_text(aes(label = label, y = fraction_exome_genome), vjust = -0.2, position = position_dodge(0.9)) +
  theme_cowplot() +
  scale_y_continuous(expand=expansion(mult=c(0,0.1)), labels = scales::label_percent()) +
  labs(x = NULL, y = "percentage of exonic mutations out of whole-genome mutations")
ratio_exome_genome

# get the number of cells
tissue_ncells_ci

# Load gene_of_interest boostdm
boostdm = fread("processed_data/boostdm/boostdm_genie_cosmic/pancancer_boostDM_intersect.txt.gz")

# filter for tissue
plot_mirror = function(genie_exp, df_point, tissue_select, labels = labels, labels_exp = labels_exp) {

  if (max(genie_exp$Protein_position, na.rm = TRUE) > 2000) {
    genie_exp = genie_exp |>
      mutate(Protein_position = (Protein_position - 1) %/% 5 + 1,
             Protein_position = Protein_position *5) |>
      group_by(Protein_position, tissue_type, tissue, model_type, type)  |>
      summarise(n = sum(n, na.rm = TRUE), .groups = "drop")
    x_label = "AA position (5AA bins)"
    } else { x_label = "AA position"}


  # way to make the plot extend both upper and lower axes
  df_point = genie_exp |>
     group_by(tissue, tissue_type, model_type, Protein_position) |>
    summarize(n = sum(n)) |>
    summarize(n = max(abs(n)) * 1.1) |>
    mutate(Protein_position = mean(genie_exp$Protein_position),
           n = ifelse(model_type == "expected", 0-n, n))

  labels_plot = left_join(labels, genie_exp, by = "tissue_type") |>
    select(tissue_type, tissue, label, model_type) |> distinct()

  labels_exp = left_join(labels_exp, genie_exp, by = "tissue_type") |>
    select(tissue_type, tissue, label, model_type) |> distinct()

  plot = ggplot(data = df_point |> filter(tissue == tissue_select), aes(x = Protein_position, y = n)) +
    geom_point(alpha = 0) +
    geom_col(data = genie_exp |>   filter(tissue == tissue_select), aes(fill = type)) +
    ggpp::geom_text_npc(data = labels_plot |> filter(tissue == tissue_select), aes(label = label, npcx = 0.05, npcy = 0.9), size = 3) +
    ggpp::geom_text_npc(data = labels_exp |>  filter(tissue == tissue_select), aes(label = label, npcx = 0.05, npcy = 0.1), size = 3) +
    facet_grid(model_type ~ tissue  , scales = "free_y")  +
    scale_fill_manual(values = COLORS6) +
    theme_cowplot() +
    panel_border() +
    theme(legend.position = "none", axis.text.y = element_text(size = 8), panel.spacing.y = unit(0, "mm"), axis.title = element_text(size = 11)) +
    labs(y = "n. expected/observed muts",  x = x_label) +
    scale_y_continuous(expand=expansion(mult=c(0,0)), labels =  function(x) label_comma()(abs(x)))
  plot
}

mirror_plot_gene_oi = function(boostdm, mean_rates, genie_data, gene_oi, drivers_only = FALSE) {

  ## genie plots
  genie_gene_oi = genie_data[Hugo_Symbol == gene_oi, ]

  gene_of_interest_counts = genie_gene_oi |>
    group_by(tissue_type) |>
    summarize(n_gene_of_interest_mut = dplyr::n())

  labels = merge(sample_counts, gene_of_interest_counts) |>
    mutate(n_tumors = format(n_tumors, big.mark = ","),
           n_gene_of_interest_mut  = format(n_gene_of_interest_mut, big.mark = ","),
           label = paste0(tissue_type, " ", n_tumors, " cases\nnumber ",  gene_oi, " muts: ", n_gene_of_interest_mut))

    # optional - filter drivers
  if (drivers_only == TRUE) {
    gene_of_interest_boostdm = boostdm[gene_name == gene_oi,]
    genie_gene_oi = genie_gene_oi |>
      dplyr::rename(pos = Start_Position, alt = Tumor_Seq_Allele2, ref = Reference_Allele) |>
      inner_join(gene_of_interest_boostdm) |>
      filter(driver == TRUE)
  }

  genie_position_counts = genie_gene_oi |>
    group_by(Protein_position, type, tissue_type) |>
    mutate(tissue_type = factor(tissue_type, levels = levels(sample_counts$tissue_type))) |>
    count()

  # correct the trinucleotide rate for a specific gene for the relative mutation load of that gene
  gene_of_interest_boostdm = boostdm[gene_name == gene_oi,]
  ratio_mutagenesis = ratios[gene_name == gene_oi, ]
  mean_rate_gene = inner_join(mean_rates, ratio_mutagenesis, by = "tissue_type")

  mutations = left_join(gene_of_interest_boostdm, mean_rate_gene, relationship = "many-to-many", by = "mut_type") |>
    left_join(triplet_match_substmodel, by = "mut_type") |>
    group_by(position, tissue_type, type) |>
    summarize(mle = sum(mle))

  mutated_cells = mutations |>
    mutate(tissue = gsub("\\n.*", "", tissue_type)) |>
    left_join(tissue_ncells_ci, by = "tissue") |>
    mutate(mle = mle * mid_estimate)  |>
    filter(!is.na(position))

  driver_mutations = left_join(gene_of_interest_boostdm[driver == TRUE, ], mean_rate_gene, relationship = "many-to-many", by = "mut_type") |>
    left_join(triplet_match_substmodel, by = "mut_type") |>
    group_by(position, tissue_type, type) |>
    summarize(mle = sum(mle))

  driver_mutated_cells = driver_mutations |>
    mutate(tissue = gsub("\\n.*", "", tissue_type)) |>
    left_join(tissue_ncells_ci, by = "tissue") |>
    mutate(mle = mle * mid_estimate)  |>
    filter(!is.na(position))

  driver_labels_exp = driver_mutated_cells |>
    group_by(tissue_type) |>
    summarize(n_mutated_cells = sum(mle)) |>
    mutate(label = paste("expected cells with driver:", format(round(n_mutated_cells), big.mark = ",")))

  labels_exp = mutated_cells |>
    group_by(tissue_type) |>
    summarize(n_mutated_cells = sum(mle)) |>
    mutate(label = paste("expected cells with any mutation:", format(round(n_mutated_cells), big.mark = ",")))

  labels_exp$label = paste0(driver_labels_exp$label, "\n", labels_exp$label)

  ### MIRROR PLOTS #####
  # combine data frames
  genie_exp = mutated_cells |>
    dplyr::rename(Protein_position = position, n = mle) |>  # rename cols so that values match with GENIE data
    bind_rows(genie_position_counts) |>
    filter(!is.na(Protein_position)) |>
    mutate(model_type = ifelse(grepl("\n", tissue_type, fixed = TRUE), "expected", "GENIE"),
           model_type = factor(model_type, levels = c("GENIE", "expected")),
           n = ifelse(model_type == "expected", 0-n, n)) |> # bind with the genie rows
    mutate(tissue = case_when(grepl("Colo|colo", tissue_type) ~ "colon",
                              grepl("lung|Lung", tissue_type) ~ "lung",
                              grepl("blood|Leukemia", tissue_type) ~ "blood",
                              grepl("Mel|mel", tissue_type) ~ "skin"))

  # plot the different plots as mirror plots
  mirror_plots = lapply(unique(genie_exp$tissue), plot_mirror, genie_exp = genie_exp, df_point = df_point, labels = labels, labels_exp = labels_exp)
  names(mirror_plots) = unique(genie_exp$tissue)
  names(mirror_plots) = names(mirror_plots)

  plots = wrap_plots(mirror_plots) + plot_annotation(title = gene_oi)
  print(gene_oi)
  return(plots)
}

genes_boostdm = most_mutated[most_mutated %in% boostdm$gene_name]

lapply(genes_boostdm, \(x) {
  plots = mirror_plot_gene_oi(boostdm = boostdm, genie_data = genie_tissue_type_all, mean_rates = mean_rates, gene_oi = x)
  ggsave(paste0("plots/genie_genes/", x, "_expected.png"), plots, width = 10, height = 10)})

pdf("plots/genie_genes/allgenes.pdf")
for (gene_oi in genes_boostdm) {
  plot = mirror_plot_gene_oi(boostdm = boostdm, genie_data = genie_tissue_type_all, mean_rates = mean_rates, gene_oi = gene_oi)
  print(plot)
}
dev.off()

# only plot the driver genes
lapply(genes_boostdm, \(x) {
  plots = mirror_plot_gene_oi(boostdm = boostdm, genie_data = genie_tissue_type_all, mean_rates = mean_rates, gene_oi = x, drivers_only = TRUE)
  ggsave(paste0("plots/genie_genes/", x, "_expected_drivers.png"), plots, width = 10, height = 10)})

pdf("plots/genie_genes/allgenes_drivers.pdf", width = 10, height = 8)
for (gene_oi in genes_boostdm) {
  plot = mirror_plot_gene_oi(boostdm = boostdm, genie_data = genie_tissue_type_all, mean_rates = mean_rates, gene_oi = gene_oi, drivers_only = TRUE)
  print(plot)
}
dev.off()

# targeted plots:
mirror_plot_gene_oi(gene_oi = "BRAF",
                    boostdm = boostdm,
                    genie_data = genie_tissue_type_all,
                    mean_rates = mean_rates, drivers_only = TRUE)

# make it possible to select hotspots:
mirror_plot_gene_oi(gene_oi = "BRAF",
                    boostdm = boostdm |> filter(gene_name == "BRAF" & position == 600),
                    genie_data = genie_tissue_type_all |> filter(Protein_position == 600),
                    mean_rates = mean_rates, drivers_only = TRUE)

mirror_plot_gene_oi(gene_oi = "BRAF",
                    boostdm = boostdm |> filter(gene_name == "BRAF" & aachange == "V600E"),
                    genie_data = genie_tissue_type_all |> filter(Protein_position == 600),
                    mean_rates = mean_rates, drivers_only = TRUE)


mirror_plot_gene_oi(gene_oi = "BRAF",
                    boostdm = boostdm |> filter(position %in% 580:620) |>
                      mutate(driver = ifelse(aachange == "V600E", TRUE, FALSE)),
                    genie_data = genie_tissue_type_all |> filter(tissue_type == "Melanoma"),
                    mean_rates = mean_rates |> filter(tissue_type == "skin\nmelanocyte"), drivers_only = TRUE)



mirror_plot_gene_oi(gene_oi = "BRAF",
                    boostdm = boostdm |>
                      mutate(driver = ifelse(aachange == "V600E", TRUE, FALSE)),
                    genie_data = genie_tissue_type_all |> filter(tissue_type == "Melanoma"),
                    mean_rates = mean_rates |> filter(tissue_type == "skin\nmelanocyte"), drivers_only = TRUE)



mirror_plot_gene_oi(gene_oi = "KRAS",
                    boostdm = boostdm |> filter(gene_name == "KRAS" & position %in% 11:13),
                    genie_data = genie_tissue_type_all |> filter(Protein_position %in% 11:13),
                    mean_rates = mean_rates, drivers_only = TRUE)

mirror_plot_gene_oi(gene_oi = "KRAS",
                    boostdm = boostdm,
                    genie_data = genie_tissue_type_all |> filter(tissue_type %in% c("Colorectal Cancer", "Non-Small Cell Lung Cancer")),
                    mean_rates = mean_rates |> filter(tissue_type %in% c("colon\nnormal", "lung\nsmoker")),
                    drivers_only = TRUE)


# make specific plots for specific positions:



# Blood make specific mirror genes for specific sites:
boostdm_ch = fread("processed_data/boostdm/boostdm_genie_cosmic/CH_boostDM_cancer.txt.gz")

genes_ch = "DNMT3A"
mean_rates_blood = mean_rates |> filter(tissue_type == "blood\nnormal")
# only plot the driver genes
mirror_plot_gene_oi(boostdm = boostdm_ch, genie_data = genie_tissue_type_all, mean_rates = mean_rates, gene_oi = "DNMT3A")

mirror_plot_gene_oi(boostdm = boostdm_ch, genie_data = genie_tissue_type_all, mean_rates = mean_rates_blood, gene_oi = "DNMT3A")

plot_DNMT3A = mirror_plot_gene_oi(boostdm = boostdm_ch, genie_data = genie_tissue_type_all |> filter(tissue_type == "Leukemia"),
                    mean_rates = mean_rates_blood, gene_oi = "DNMT3A")



# Todo in the future: Make a final 'figure' plot including all other changes
# if not, make at least a plot for the top part of the plot.

# load the bowel cancer data:
crc_freq = fread("raw_data/UKBiobank/colorectal_cancer_frequency_UKB.csv")

# calculate the incidence rates
ukbiobank_crc = data.frame(age = 0:max(crc_freq$current_age)) |>
  mutate(n_alive = sapply(age, \(x) sum(crc_freq$current_age >= x)),
         n_tumor = sapply(age, \(x) sum(crc_freq$var_Colorectal_age == x, na.rm = TRUE)),
         n_no_tumor = n_alive - n_tumor,
         risk = n_tumor / n_alive,
         cumulative_risk = cumsum(risk)) |>
  filter(n_alive > 5000)

# tumor incidence
n_cohort = ggplot(ukbiobank_crc, aes(x = age)) +
  geom_col(aes(y = n_alive)) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y = "number of individuals in cohort", title = "Number of individuals in the cohort") +
  theme(plot.title = element_text(hjust = 0.5))

# tumor incidence
tumor_incidence = ggplot(ukbiobank_crc, aes(x = age)) +
  geom_col(aes(y = n_tumor))  +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y = "Tumor incidence", title = "Tumor Incidence") +
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
  geom_col(aes(y = cumulative_risk)) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = label_percent()) +
  labs(y = "Risk percentage", title = "Cumulative Incidence") +
  theme(plot.title = element_text(hjust = 0.5))

n_cohort + tumor_incidence + yearly_incidence_rates + cumulative_incidence

# bowel_cancer = readRDS("processed_data/colon/bowel_cancer.rds")
ukbiobank_crc_plot = ukbiobank_crc |>
  filter(n_tumor != 0) |>
  select(age, cumulative_risk) |>
  pivot_longer(cumulative_risk, names_to = "name", values_to = "incidence")


vogelgram_risks_long = readRDS("processed_data/colon/vogelgram_plot.rds") |>
  dplyr::rename(name = mutation, incidence = expected_mutated_cells) |>
  select(age, name, incidence)

df_vogelgram_incidence = rbind(ukbiobank_crc_plot, vogelgram_risks_long) |>
  mutate(incidence_100k = incidence * 1e5,
         name = factor(name, levels = c("KRAS_driver", "APC 2x SNV", "APC 2x +KRAS", "cumulative_risk")))

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

mirror_plots_all = mirror_plots[[1]] + mirror_plots[[2]] + mirror_plots[[3]] + mirror_plots[[4]]

left = table_1 / exp_types_percent / exp_types  &
  theme(plot.margin = margin(8,8,8,8,"mm"))
figure_5 = (left | mirror_plots_all | (percent_mutated / vogel_plot)) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(1, 3, 1.3)) &
  theme(plot.margin = margin(8,7,8,7,"mm"), plot.tag = element_text(face = 'bold'))
ggsave("plots/manuscript/main_figures/figure_5.svg", figure_5, width = 23, height = 13, bg = "white", dpi = 600)
ggsave("plots/manuscript/main_figures/figure_5.png", figure_5, width = 23, height = 11, bg = "white", dpi = 600)

