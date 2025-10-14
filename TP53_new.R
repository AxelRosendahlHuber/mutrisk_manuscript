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
  `colnames<-`(c("Tumor_Sample_Barcode", "age", "tissue_category"))
genie_tissue_category_all = inner_join(genie_samples, genie_data)
fwrite(genie_tissue_category_all, "processed_data/GENIE_17/GENIE_17_genie_tissue_category.txt.gz")

# Count number of tumors in the GENIE data
sample_counts = genie_samples |>
  group_by(tissue_category) |>
  summarize(n_tumors = dplyr::n()) |>
  mutate(tissue_category = factor(tissue_category, levels = c("Leukemia", "Colorectal Cancer", "Non-Small Cell Lung Cancer", "Melanoma")))

plot_gene_oi = function(genie_data, gene_oi) {

  print(gene_oi)
  genie_gene_oi = genie_tissue_category_all[Hugo_Symbol == gene_oi, ]

  gene_of_interest_counts = genie_gene_oi |>
    group_by(tissue_category) |>
    summarize(n_gene_of_interest_mut = dplyr::n())

  labels = merge(sample_counts, gene_of_interest_counts) |>
    mutate(n_tumors = format(n_tumors, big.mark = ","),
           n_gene_of_interest_mut  = format(n_gene_of_interest_mut, big.mark = ","),
           label = paste0(tissue_category, "\nnumber tumors: ", n_tumors, "\nnumber ",  gene_oi, " muts: ", n_gene_of_interest_mut))

  genie_position_counts = genie_gene_oi |>
    group_by(Protein_position, type, tissue_category) |>
    mutate(tissue_category = factor(tissue_category, levels = levels(sample_counts$tissue_category))) |>
    count()

  genie_total_plot = ggplot(genie_position_counts, aes(x = Protein_position, y = n, fill = type)) +
    geom_col() +
    ggpp::geom_text_npc(data = labels, aes(label = label, npcx = 0.01, npcy = 0.9), size = 3) +
    facet_grid(tissue_category ~ . , scales = "free_y") +
    scale_fill_manual(values = COLORS6) +
    theme_cowplot() +
    panel_border() +
    theme(legend.position = "none", strip.text = element_blank(), axis.text.y = element_text(size = 9)) +
    labs(y = "number of observed mutations", subtitle = paste(gene_oi, "mutations across tissue_categorys - GENIE data"), x = "AA position") +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)), labels = label_comma())
  ggsave(paste0("plots/genie_genes/GENIE_", gene_oi, "_tissue_categorys.png"), genie_total_plot, width = 10, height = 6, bg = "white")
  genie_total_plot
  # get the consequence
  genie_gene_oi |> count(tissue_category, Consequence)
}

most_mutated = genie_tissue_category_all |> count(Hugo_Symbol) |> arrange(desc(n)) |>
  head(5) |> pull(Hugo_Symbol)

gene_oi = "TP53"

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

# mutation rates
mutation_rates = expected_rates |>
  mutate(tissue_category = paste0(tissue, "_", category)) |>
  left_join(metadata)

# Mean mutation rate for each trinucleotide
# Check if the blood rates actually make sense - this seems too low - this must be because of the cord blood donors
mean_rates = mutation_rates |>
  group_by(tissue_category, mut_type) |>
  summarize(mle = mean(mle),
            mean_age = mean(age)) |>
  mutate(tissue_category_age = paste0(tissue, "\n(", format(mean_age, digits = 3, nsmall = 1 ), ") years"))

# get rates for a plot
# add the mean rates for the cohort -> make sure that we weigh each of the donors equally..
mean_rates_table = mutation_rates |>
  group_by(tissue_category, mut_type, donor, age) |>
  summarize(mle = mean(mle)) |>
  group_by(tissue_category) |>
  summarize(`Mean Age` = round(mean(age), 1)) |>
  mutate(Tissue = gsub("\\\n.*", "",  tissue_category)) |>
  select(Tissue, `Mean Age`) |>
  distinct()
#mean_rates_table$`Cell Type` = c("Hematopoietic stem cells", "Crypt stem cells", "Lung basal cells", "Melanocytes")
table_1 = wrap_table(mean_rates_table, space = unit(0.1, "npc"))

# determine the ratio of exonic to genomic mutations
# load the number of WGS genomes
select_categories = c("smoker", "normal")
wgs_muts = list.files("processed_data/", pattern = "cell_muts.tsv", full.names = TRUE, recursive = TRUE)
names(wgs_muts) = c("blood", "colon", "lung")
muts = lapply(wgs_muts, \(x) fread(x)[,c("sampleID", "chr", "pos", "ref", "alt")]) |>
  rbindlist(idcol = "tissue_category")

wgs_muts = muts |>
  group_by(tissue_category) |>
  summarize(n_wgs = dplyr::n())

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
      group_by(Protein_position, tissue_category, tissue, model_type, type)  |>
      summarise(n = sum(n, na.rm = TRUE), .groups = "drop")
    x_label = "AA position (5AA bins)"
    } else { x_label = "AA position"}


  # way to make the plot extend both upper and lower axes
  df_point = genie_exp |>
     group_by(tissue, tissue_category, model_type, Protein_position) |>
    summarize(n = sum(n)) |>
    summarize(n = max(abs(n)) * 1.1) |>
    mutate(Protein_position = mean(genie_exp$Protein_position),
           n = ifelse(model_type == "expected", 0-n, n))

  labels_plot = left_join(labels, genie_exp, by = "tissue_category") |>
    select(tissue_category, tissue, label, model_type) |> distinct()

  labels_exp = left_join(labels_exp, genie_exp, by = "tissue_category") |>
    select(tissue_category, tissue, label, model_type) |> distinct()

  plot = ggplot(data = df_point |> filter(tissue == tissue_select), aes(x = Protein_position, y = n)) +
    geom_point(alpha = 0) +
    geom_col(data = genie_exp |>   filter(tissue == tissue_select), aes(fill = type)) +
    ggpp::geom_text_npc(data = labels_plot |> filter(tissue == tissue_select), aes(label = label, npcx = 0.05, npcy = 0.9), size = 3) +
    ggpp::geom_text_npc(data = labels_exp |>  filter(tissue == tissue_select), aes(label = label, npcx = 0.05, npcy = 0.1), size = 3) +
    facet_grid(factor(model_type, levels = c("GENIE", "expected")) ~ tissue  , scales = "free_y")  +
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
    group_by(tissue_category) |>
    summarize(n_gene_of_interest_mut = dplyr::n())

  labels = merge(sample_counts, gene_of_interest_counts) |>
    mutate(n_tumors = format(n_tumors, big.mark = ","),
           n_gene_of_interest_mut  = format(n_gene_of_interest_mut, big.mark = ","),
           label = paste0(tissue_category, " ", n_tumors, " cases\nnumber ",  gene_oi, " muts: ", n_gene_of_interest_mut))

    # optional - filter drivers
  if (drivers_only == TRUE) {
    gene_of_interest_boostdm = boostdm[gene_name == gene_oi,]
    genie_gene_oi = genie_gene_oi |>
      dplyr::rename(pos = Start_Position, alt = Tumor_Seq_Allele2, ref = Reference_Allele) |>
      inner_join(gene_of_interest_boostdm) |>
      filter(driver == TRUE)
  }

  genie_position_counts = genie_gene_oi |>
    group_by(Protein_position, type, tissue_category) |>
    mutate(tissue_category = factor(tissue_category, levels = levels(sample_counts$tissue_category))) |>
    count() |>
    mutate(model_type = "GENIE")

  # correct the trinucleotide rate for a specific gene for the relative mutation load of that gene
  gene_of_interest_boostdm = boostdm[gene_name == gene_oi,]
  ratio_mutagenesis = ratios[gene_name == gene_oi, ] |>
    mutate(tissue_category = paste0(tissue, "_", category))
  mean_rate_gene = inner_join(mean_rates, ratio_mutagenesis)

  mutations = left_join(gene_of_interest_boostdm, mean_rate_gene, relationship = "many-to-many", by = "mut_type") |>
    left_join(triplet_match_substmodel, by = "mut_type") |>
    group_by(position, tissue_category, type) |>
    summarize(mle = sum(mle))

  mutated_cells = mutations |>
    mutate(tissue = gsub("_.*", "", tissue_category)) |>
    left_join(tissue_ncells_ci) |>
    mutate(mle = mle * mid_estimate)  |>
    filter(!is.na(position))

  driver_mutations = left_join(gene_of_interest_boostdm[driver == TRUE, ], mean_rate_gene, relationship = "many-to-many", by = "mut_type") |>
    left_join(triplet_match_substmodel, by = "mut_type") |>
    group_by(position, tissue_category, type) |>
    summarize(mle = sum(mle))

  driver_mutated_cells = driver_mutations |>
    mutate(tissue = gsub("_.*", "", tissue_category)) |>
    left_join(tissue_ncells_ci, by = "tissue") |>
    mutate(mle = mle * mid_estimate)  |>
    filter(!is.na(position))

  driver_labels_exp = driver_mutated_cells |>
    group_by(tissue_category) |>
    summarize(n_mutated_cells = sum(mle)) |>
    mutate(label = paste("expected cells with driver:", format(round(n_mutated_cells), big.mark = ",")))

  labels_exp = mutated_cells |>
    group_by(tissue_category) |>
    summarize(n_mutated_cells = sum(mle)) |>
    mutate(label = paste("expected cells with any mutation:", format(round(n_mutated_cells), big.mark = ",")))

  labels_exp$label = paste0(driver_labels_exp$label, "\n", labels_exp$label)

  ### MIRROR PLOTS #####
  # combine data frames
  genie_exp = mutated_cells |>
    mutate(model_type = "expected") |>
    ungroup() |>
    dplyr::rename(Protein_position = position, n = mle) |>  # rename cols so that values match with GENIE data
    bind_rows(genie_position_counts) |>
    filter(!is.na(Protein_position)) |>
    mutate(n = ifelse(model_type == "expected", 0-n, n)) |> # bind with the genie rows
    mutate(tissue = case_when(grepl("Colo|colo", tissue_category) ~ "colon",
                              grepl("lung|Lung", tissue_category) ~ "lung",
                              grepl("blood|Leukemia", tissue_category) ~ "blood",
                              grepl("Mel|mel", tissue_category) ~ "skin"))

  # plot the different plots as mirror plots
  mirror_plots = lapply(unique(genie_exp$tissue), plot_mirror, genie_exp = genie_exp, df_point = df_point, labels = labels, labels_exp = labels_exp)
  names(mirror_plots) = unique(genie_exp$tissue)
  names(mirror_plots) = names(mirror_plots)

  plots = wrap_plots(mirror_plots) + plot_annotation(title = gene_oi)
  print(gene_oi)
  return(plots)
}

genes_boostdm = most_mutated[most_mutated %in% boostdm$gene_name]
genes_boostdm = "TP53"

lapply(genes_boostdm, \(x) {
  plots = mirror_plot_gene_oi(boostdm = boostdm, genie_data = genie_tissue_category_all, mean_rates = mean_rates, gene_oi = x)
  ggsave(paste0("plots/genie_genes/", x, "_expected.png"), plots, width = 10, height = 10)})

pdf("plots/genie_genes/allgenes.pdf")
for (gene_oi in genes_boostdm) {
  plot = mirror_plot_gene_oi(boostdm = boostdm, genie_data = genie_tissue_category_all, mean_rates = mean_rates, gene_oi = gene_oi)
  print(plot)
}
dev.off()

# only plot the driver genes
lapply(genes_boostdm, \(x) {
  plots = mirror_plot_gene_oi(boostdm = boostdm, genie_data = genie_tissue_category_all, mean_rates = mean_rates, gene_oi = x, drivers_only = TRUE)
  ggsave(paste0("plots/genie_genes/", x, "_expected_drivers.png"), plots, width = 10, height = 10)})

pdf("plots/genie_genes/allgenes_drivers.pdf", width = 10, height = 8)
for (gene_oi in genes_boostdm) {
  plot = mirror_plot_gene_oi(boostdm = boostdm, genie_data = genie_tissue_category_all, mean_rates = mean_rates, gene_oi = gene_oi, drivers_only = TRUE)
  print(plot)
}
dev.off()

# make specific plots for specific positions:
# Blood make specific mirror genes for specific sites:
boostdm_ch = fread("processed_data/boostdm/boostdm_genie_cosmic/CH_boostDM_cancer.txt.gz")

# UKBiobank DNMT3A mutations:
UKB_DNMT3A_muts = fread("raw_data/UKBiobank/UkBiobank_DNMT3A_mut_age.csv")
UKB_DNMT3A_counts = UKB_DNMT3A_muts[, .N, by = c("aa_change", "REF", "ALT")]  |>
  mutate(position = parse_number(aa_change),
         type = paste0(REF, ">", ALT),
         type = case_match(type, .default = type,
                           "A>T" ~ "T>A", "A>G" ~ "T>C", "A>C" ~ "T>G",
                           "G>T" ~ "C>A", "G>A" ~ "C>T", "G>C" ~ "C>G")) |>
  mutate(mrate = N, tissue_category = "UKBiobank CH") |>
  select(position, type, tissue_category, mrate)

# first attempt ukbiobank plots:
DNMT3A_plot = UKB_DNMT3A_counts |>
  ggplot(aes(x = position, y = N, fill = type)) +
  geom_col() +
  scale_fill_manual(values = COLORS6) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x  = NULL , y = "Number of CH mutations\nobserved in the UKBiobank cohort", title = "DNMT3A", fill = NULL)


genes_ch = "DNMT3A"
mean_rates_blood = mean_rates |> filter(tissue_category == "blood_normal")

# only plot the driver genes
DNMT3A = boostdm_ch |> filter(gene_name == "DNMT3A")

mutations_blood_DNMT3A = left_join(DNMT3A, mean_rates_blood, relationship = "many-to-many", by = "mut_type") |>
  left_join(triplet_match_substmodel, by = "mut_type") |>
  group_by(position, tissue_category, type) |>
  summarize(mrate = sum(mle)) |>
  select(position, type, tissue_category, mrate)

library(ggh4x)
df_mirror = bind_rows(mutations_blood_DNMT3A, UKB_DNMT3A_counts) |>
  mutate(
    tissue_category = ifelse(tissue_category == "blood_normal", "Expected mutrate\nblood", tissue_category),
    tissue_category = factor(tissue_category, levels = c("UKBiobank CH", "Expected mutrate\nblood")),
    mrate = ifelse(tissue_category == "Expected mutrate\nblood", 0-mrate, mrate)) |>
  ungroup()


# way to make the plot extend both upper and lower axes
df_point = df_mirror |>
  group_by(tissue_category, position) |>
  summarize(mrate = sum(mrate), .groups = "drop_last") |>
  summarize(mrate = max(abs(mrate)) * 1.1) |>
  mutate(position = 500,
         mrate = ifelse(tissue_category == "Expected mutrate\nblood", 0-mrate, mrate)) |>
  ungroup()

F5B = ggplot(df_point, aes(x = position, y = mrate)) +
  geom_point(color = "white") +
  geom_col(data = df_mirror, aes(fill = type)) +
  geom_text(data = data.frame(tissue_category = factor("UKBiobank CH"), position = 50, mrate = 1500, label = "TP53"),
            aes(label = label)) +
  facet_grid2(tissue_category ~ . , scales = "free") +
  scale_fill_manual(values = COLORS6) +
  theme_cowplot() +
  panel_border() +
  theme(legend.position = "none", panel.spacing.y = unit(0, "mm")) +
  labs(y = "Number expected/\nobserved muts",  x = "AA position") +
  scale_y_continuous(expand=expansion(mult=c(0,0)), breaks = scales::breaks_extended(n = 3))