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
select_aachange = "G12V"

KRAS_G12V_hotspot = CH_bDM[ gene_name == select_gene & aachange == select_aachange,
                                c("gene_name", "mut_type", "aachange", "position", "driver")]
expected_rate_KRAS_G12V = expected_rates |>
  inner_join(KRAS_G12V_hotspot, by = "mut_type") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells)) |>
  group_by(donor, category) |>
  summarise(across(c(mle, cilow, cihigh, age), mean))

# calculate the mutational frequencies for the mean KRAS rates:
KRAS_G12_hotspot = CH_bDM[gene_name == select_gene & position == 12,c("gene_name", "mut_type", "aachange", "position", "driver")]
mean_KRAS_rates = expected_rates |>
  inner_join(KRAS_G12_hotspot, by = "mut_type") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells)) |>
  group_by(category, aachange) |>
  summarise(across(c(mle, cilow, cihigh, age), mean))
fwrite(mean_KRAS_rates,paste0("processed_data/",tissue,"/",tissue,"_",select_gene,"_",select_aachange,"tsv"))

figure_2a_KRAS = ggplot(expected_rate_KRAS_G12V,  aes(x = age, y = mle, fill = category)) +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh, color = category), width = 0, show.legend = FALSE) +
  geom_point(shape = 21, size = 2.4, color = "white", stroke = 0.3) +
  scale_fill_manual(values = blood_colors) +
  scale_color_manual(values = blood_colors) +
  theme_cowplot() +
  labs(subtitle = "N cells with KRAS G12V", y = "number of cells", x = "Age (years)", fill = NULL) +
  theme(legend.position = "inside", legend.position.inside = c(0.05, 0.8))
figure_2a_KRAS

# Get the actual values for the graphic to be correct
DNMT3A_R882H_hotspot = CH_bDM[aachange == "R882C" & gene_name == "DNMT3A", c("gene_name", "mut_type", "aachange", "position", "driver")]
DNMT3A_R882H_hotspot |> left_join(triplet_match_substmodel)

DNMT3A_R882H_hotspot = CH_bDM[aachange == "R882H" & gene_name == "DNMT3A",c("gene_name", "mut_type", "aachange", "position", "driver")]
expected_rate_DNMT3A_R882H = expected_rates |>
  inner_join(DNMT3A_R882H_hotspot, by = "mut_type") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells)) |>
  group_by(donor, category) |>
  summarise(across(c(mle, cilow, cihigh, age), mean))

figure_2a = ggplot(expected_rate_DNMT3A_R882H,
                   aes(x = age, y = mle, fill = category)) +
  geom_smooth(method = "lm", color = "black") +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh, color = category),
                width = 0,
                show.legend = FALSE) +
  geom_point(shape = 21, size = 2.4, color = "white", stroke = 0.3) +
  ggpubr::stat_cor() +
  scale_fill_manual(values = blood_colors) +
  scale_color_manual(values = blood_colors) +
  theme_cowplot() +
  labs(subtitle = "N cells with DNMT3A R882H*",y = "number of cells",
       x = "Age (years)",fill = NULL) +
  theme(legend.position = "inside", legend.position.inside = c(0.05, 0.8))
figure_2a
ggsave("plots/blood/Figure_2A.png", width = 5, height = 4.5, bg = "white")

# Get the actual values for the graphic to be correct
DNMT3A_R882H_hotspot = CH_bDM[aachange == "R882C" & gene_name == "DNMT3A", c("gene_name", "mut_type", "aachange", "position", "driver")]
DNMT3A_R882H_hotspot |> left_join(triplet_match_substmodel)

DNMT3A_R882H_hotspot = CH_bDM[gene_name == "DNMT3A" & driver == TRUE  & aachange == "R882H" ,c("gene_name", "mut_type", "aachange", "position", "driver")]
expected_rate_DNMT3A_R882H = expected_rates |>
  inner_join(DNMT3A_R882H_hotspot, by = "mut_type") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells)) |>
  group_by(donor, category) |>
  summarise(across(c(mle, cilow, cihigh), sum),
            age = mean(age))

figure_2a = ggplot(expected_rate_DNMT3A_R882H,
                   aes(x = age, y = mle, fill = category)) +
  geom_smooth(method = "lm", color = "black") +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh, color = category),
                width = 0,
                show.legend = FALSE) +
  geom_point(shape = 21, size = 2.4, color = "white", stroke = 0.3) +
  ggpubr::stat_cor() +
  scale_fill_manual(values = blood_colors) +
  scale_color_manual(values = blood_colors) +
  theme_cowplot() +
  labs(subtitle = "N cells with DNMT3A driver mutation",y = "number of cells",
       x = "Age (years)",fill = NULL) +
  theme(legend.position = "inside", legend.position.inside = c(0.05, 0.8))
figure_2a
ggsave("plots/blood/Figure_2A.png", width = 5, height = 4.5, bg = "white")

###
gene_oi = "PPM1D"
TP53_drivers_hotspot = CH_bDM[gene_name == gene_oi,c("gene_name", "mut_type", "aachange", "position", "driver")]
expected_rate_TP53_drivers = expected_rates |>
  inner_join(TP53_drivers_hotspot, by = "mut_type") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells)) |>
  group_by(donor, category, age, sampleID) |>
  summarise(across(c(mle, cilow, cihigh), sum), .groups = "drop_last") |>
  summarise(across(c(mle, cilow, cihigh), mean))

TP53_drivers_plot = ggplot(expected_rate_TP53_drivers,
                           aes(x = age, y = mle, fill = category)) +
  geom_smooth(method = "lm", color = "black") +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh, color = category),
                width = 0,
                show.legend = FALSE) +
  geom_point(shape = 21, size = 2.4, color = "white", stroke = 0.3) +
  ggpubr::stat_cor() +
  scale_fill_manual(values = blood_colors) +
  scale_color_manual(values = blood_colors) +
  theme_cowplot() +
  labs(subtitle = paste("N cells with", gene_oi,  "driver mutations*") ,y = "number of cells",
       x = "Age (years)",fill = NULL) +
  theme(legend.position = "inside", legend.position.inside = c(0.05, 0.8))
TP53_drivers_plot

# print the DNMT3A hotspot rate for
single_cell_rate = expected_rate_DNMT3A_R882H |>
  mutate(across(c(mle, cilow, cihigh), ~ ./ ncells))
R882H_rate = lm(mle ~ age, single_cell_rate)

# barplot for the individual mutation consequences
DNMT3A_counts_consequence = CH_bDM[gene_name == "DNMT3A",.N, by = c("gene_name", "mut_type", "consequence")]
expected_DNMT3A_muts = expected_rates |>
  left_join(DNMT3A_counts_consequence, relationship = "many-to-many") |>
  inner_join(ratios, relationship = "many-to-many") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N * ncells * ratio)) |>
  group_by(consequence, category, sampleID) |>
  summarize(across(c(mle, cilow, cihigh), sum)) |>
  group_by(consequence, category) |>
  summarize(across(c(mle, cilow, cihigh), mean))

DNMT3A_counts_boostdm_ch = CH_bDM[gene_name == "DNMT3A", .N, by = c("gene_name", "mut_type", "driver")]
expected_DNMT3A_muts = expected_rates |>
  left_join(DNMT3A_counts_boostdm_ch, relationship = "many-to-many") |>
  inner_join(ratios, relationship = "many-to-many") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N * ncells * ratio)) |>
  group_by(driver, category, sampleID) |>
  summarize(across(c(mle, cilow, cihigh), sum)) |>
  group_by(driver, category) |>
  summarize(across(c(mle, cilow, cihigh), mean))

# Expected number of cells with double mutations:
DNMT3A_double_muts = expected_rates |>
  left_join(ratios |> filter(gene_name == "DNMT3A")) |>
  left_join(metadata) |>
  inner_join(DNMT3A_counts_boostdm_ch |> filter(driver),
             by = "mut_type",  relationship = "many-to-many") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio)) |>
  group_by(category, donor, age, mut_type) |>
  summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop_last") |>
  summarize(across(c(mle, cilow, cihigh), sum)) |>
  mutate(across(c(mle, cilow, cihigh), ~ ((.^2) / 2) * ncells))

#### profiles:
DNMT3A_blood_count = genie_blood[gene_name == "DNMT3A", .N, by = "position"]
genie_DNMT3A = ggplot(genie_blood[gene_name == "DNMT3A", ], aes(x = position, y = 1)) +
  geom_col(mapping = aes(fill = type)) +
  geom_point(data = DNMT3A_blood_count, aes(x = position, y = N)) +
  scale_fill_manual(values = COLORS6) +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = comma) +
  labs(y = "number of mutations", title = "GENIE DNMT3A colorectal cancer data", x = NULL)

# profile of the boostdm_ch rates:
DNMT3A_muts = expected_rates |>
  group_by(category, mut_type) |>
  summarize(across(c(mle, cilow, cihigh), mean)) |>
  left_join(ratios[gene_name == "DNMT3A"]) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ncells * ratio)) |>
  left_join(triplet_match_substmodel) |>
  left_join(CH_bDM[gene_name == "DNMT3A", ], by = c("mut_type", "gene_name"), relationship = "many-to-many") |>
  mutate(boostdm_ch = ifelse(driver, "driver", "non-driver")) |>
  setDT()

DNMT3A_count = DNMT3A_muts[, .(mle = sum(mle)), by = c("position", "category", "boostdm_ch")]
normal_DNMT3A = DNMT3A_muts[category == "normal", ]
normal_DNMT3A_count = DNMT3A_count[category == "normal", ]

expected_DNMT3A = ggplot(normal_DNMT3A, aes(x = position, y = mle)) +
  geom_col(aes(fill = type)) +
  geom_point(data = normal_DNMT3A_count) +
  scale_fill_manual(values = COLORS6) +
  labs(y = "numer of mutated cells",title = "expected number of HSCs mutated: DNMT3A", fill = NULL) +
  facet_grid(boostdm_ch ~ .) +
  theme_cowplot() +
  panel_border() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = comma) +
  theme(legend.position = "inside", legend.position.inside = c(0.05, 0.75), legend.key.size = unit(4, "mm"))

normal_DNMT3A_count_sum = normal_DNMT3A_count |>
  group_by(position, category) |>
  summarize(mle = sum(mle))

# get the mutation profiles for all the different cases:
DNMT3A_count_total = DNMT3A_count |>
  filter(!is.na(position)) |>
  group_by(category, position) |>
  summarise(mle = sum(mle))
#####
# get the level of mutations for DNMT3A drivers, TP53 drivers, all CH drivers.
####


####
# DNMT3A driver muts
####
site_freqs_DNMT3A_drivers = CH_bDM[gene_name == "DNMT3A" & driver == TRUE, .N, by = c("mut_type", "gene_name")]
mrate_DNMT3A_drivers = get_gene_rate( exp_rates = expected_rates, metadata = metadata,
                                      site_freqs = site_freqs_DNMT3A_drivers, ratios = ratios, ncells = ncells)
plot_mrate(mrate_DNMT3A_drivers, title = "number of cells with DNMT3A driver mutations", colors = blood_colors)
ggsave("plots/blood/ncells_DNMT3A_drivers.png", width = 5.5, height = 4.2, bg = "white")

# numbers for the manuscript:
mrate_DNMT3A_drivers |> filter(age > 0) |> pull(mle) |> min()
mrate_DNMT3A_drivers |> filter(age > 0) |> pull(mle) |> max()


mrate_DNMT3A_drivers_plot = plot_mrate(mrate_DNMT3A_drivers |> filter(category == "normal"), title = "number of cells with DNMT3A driver mutations",
                                       colors = blood_colors)
ggsave("plots/blood/ncells_DNMT3A_drivers_normal.png", width = 5.5, height = 4.2, bg = "white")

# taking the actual difference between individual clones
mrate_DNMT3A_drivers1 = get_gene_rate_fraction(exp_rates = expected_rates,
                                               metadata = metadata,
                                               site_freqs = site_freqs_DNMT3A_drivers,
                                               ratios = ratios,  ncells = ncells)
plot_mrate(mrate_DNMT3A_drivers1 |> filter(category == "normal"),
           title = "number of cells with DNMT3A driver mutations - normal tissue",
           colors = blood_colors)
ggsave("plots/blood/ncells_DNMT3A_drivers_normal_points.png", width = 8, height = 5,  bg = "white")

# get the driver mutations for a list of genes, get the expected rates for a set of different values:
genes_of_i = c("ASXL1", "CHEK2", "TP53", "PPM1D", "DNMT3A", "TET2")
total_list = list()
for (gene_oi in genes_of_i) {

  site_freqs_drivers = CH_bDM[gene_name %in% gene_oi & driver == TRUE, .N, by = c("mut_type", "gene_name")]
  mrate_drivers_100k = get_gene_rate(exp_rates = expected_rates,metadata = metadata,
                                          site_freqs = site_freqs_drivers, ratios = ratios, ncells = ncells) |>
    mutate(ncells = 1e5, cell_level = "estimate")
  mrate_drivers_25 = get_gene_rate(exp_rates = expected_rates,metadata = metadata,
                                        site_freqs = site_freqs_drivers, ratios = ratios, ncells = 2.5e3) |>
    mutate(ncells = 2.5e3,  cell_level = "low")
  mrate_drivers_1.3m = get_gene_rate(exp_rates = expected_rates,metadata = metadata,
                                          site_freqs = site_freqs_drivers, ratios = ratios, ncells = 1.3e6) |>
    mutate(ncells = 1.3e6, cell_level = "high")
  total_cells = rbind(mrate_drivers_25, mrate_drivers_100k, mrate_drivers_1.3m)
  total_list[[gene_oi]] = total_cells
}

# sets to compare with the plot of Masha
gene_rates = rbindlist(total_list, idcol = "gene")
label = gene_rates |> group_by(gene_name, cell_level) |>
  filter(mle == max(mle))
ggplot(gene_rates, aes(x = age, y = mle, group = cell_level)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE, color = 'darkred') +
  geom_text(data = label, aes(label = format(mle, digits = 1, scientific = FALSE ), hjust = -0.1)) +
  facet_wrap(gene_name ~ . , scales = "free") +
  theme_cowplot() +
  scale_x_continuous(limits = c(0, 85)) +
  labs(x = "Age (years)", y = "number of expected driver mutations")

###
# signature-specific modeling of DNMT3A mutation rates
###

# Load signature specific mutation rates
signature_rates_files = "processed_data/blood/blood_sig_donor_rates.tsv.gz"
names(signature_rates_files) = gsub("_sig_patient_rates.tsv.gz", "",basename(signature_rates_files))
signature_rates = fread(signature_rates_files) |>
  mutate(category = "normal") |>
  dplyr::rename(mut_type = name, mle = value) |>
  left_join(metadata |> select(donor, age, sampleID) |> distinct())

mrate_sigs_DNMT3A_driver = get_gene_rate_sig(exp_rates = signature_rates, metadata,
                                             site_freqs = site_freqs_DNMT3A_drivers,ratios = ratios)
mrate_sigs_DNMT3A_driver |>
  mutate(mle = mle * ncells) |>
  group_by(age, signature, category) |>
  summarize(mle = mean(mle)) |>
  ggplot(aes(x = age, y = mle, fill = signature)) +
  geom_col() +
  facet_grid(. ~ category) +
  ggsci::scale_fill_igv() +
  theme_cowplot() +
  panel_border() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = comma) +
  theme(legend.position = "inside", legend.position.inside = c(0.05, 0.65),
        legend.key.size = unit(4, "mm"), legend.title = element_blank()  ) +
  labs(y = "number of mutated cells", x = "Age (years)", title = "DNMT3A driver mutatoins by signature")
ggsave("plots/blood/number_of_cells_mutated_by_sig.png", width = 10, height = 4.2, bg = "white")

# get hotspot mutation levels for blood
DNMT3A_R882H_hotspot_sig_rates = DNMT3A_R882H_hotspot |>
  left_join(signature_rates, by = "mut_type") |>
  left_join(metadata) |>
  group_by(age, signature, category) |>
  summarize(mle = mean(mle))

ggplot(DNMT3A_R882H_hotspot_sig_rates |> filter(category == "normal"),
       aes(x = age, y = mle * ncells, fill = signature)) +
  geom_col() +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ggsci::scale_fill_igv() +
  labs(subtitle = "bloodic crypts with R882H* DNMT3A hotspot mutations",
       y = "Number of cells", x = "Age (years)")
ggsave("plots/blood/R882H_stop_rates_normal.png",width = 5,height = 4,bg = "white")

ggplot(DNMT3A_R882H_hotspot_sig_rates,   aes(x = age, y = mle * ncells, fill = signature)) +
  geom_col() +
  facet_grid(. ~ category) +
  theme_cowplot() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ggsci::scale_fill_igv() +
  labs(subtitle = "bloodic crypts with R882H* DNMT3A hotspot mutations",
       y = "Number of cells",  x = "Age (years)")
ggsave("plots/blood/R882H_stop_rates_all.png", width = 10, height = 4, bg = "white")

####
# mutation rates for individual clones by donor
####
DNMT3A_total = CH_bDM[gene_name == "DNMT3A", .N, by = c("mut_type", "gene_name", "consequence")]
consequence_rates = DNMT3A_total |>
  left_join(signature_rates |> filter(category == "normal"), relationship = "many-to-many",  by = "mut_type") |>
  left_join(ratios) |>
  mutate(mle = mle * N * ncells * ratio) |>
  group_by(consequence, sampleID, signature) |>
  summarize(mle = sum(mle))

# DNMT3A mut types drivers
DNMT3A_total = CH_bDM[gene_name == "DNMT3A",.N, by = c("mut_type", "gene_name", "driver")]
driver_rates = DNMT3A_total |>
  left_join(signature_rates, relationship = "many-to-many", by = "mut_type") |>
  left_join(ratios, ) |>
  mutate(mle = mle * N * ncells * ratio) |>
  group_by(driver, sampleID, signature) |>
  summarize(mle = sum(mle)) |>
  arrange(driver, sampleID, signature)

driver_muts = driver_rates |>
  left_join(metadata) |>
  group_by(category, age, signature, driver) |>
  summarize(mle = mean(mle))

DNMT3A_drivers_by_sig = driver_muts |>
  mutate(boostdm_ch = ifelse(driver, "driver", "no driver")) |>
  ggplot(aes(x = age, y = mle, fill = signature)) +
  geom_col(position = "stack") +
  facet_grid(category ~ boostdm_ch, scales = "free_y") +
  ggsci::scale_fill_igv() +
  labs(y = "number of mutated cells", subtitle = "DNMT3A mutations by position by signature",
       x = "Age (years)", fill = NULL) +
  theme_cowplot() + panel_border() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.05, 0.7)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = comma)
ggsave("plots/blood/driver_muts_by_driver_state.png", DNMT3A_drivers_by_sig,  width = 7, height = 4)

# Read the UKbiobank data:
# make figures for CH, DNMT3A and DNMT3A R882H mutation
DNMT3A_age_frequencies = fread("raw_data/UKBiobank/UKB_age_frequencies_DNMT3A.tsv")
CH_age_frequencies = fread("raw_data/UKBiobank/UKB_age_frequencies.tsv")

CH_frequencies = CH_age_frequencies |>
  mutate(CH = rowSums(select(CH_age_frequencies, c(-Age, -Individuals)))) |>
  mutate(fraction = CH / Individuals) |>
  filter(Individuals > 100)

# arrange plots using patchwork
# example: Age = 0.
df = data.frame(age = 1:100, ncells_mut = c(1, rep(NA, 99)))

# exponentially increase cells:
for (i in 2:nrow(df)) {
  df[i, 2] = df[i - 1, 2] * 1.18
  print(paste(age = 1, df[i, 2]))
}

ggplot(df) +
  geom_line(aes(x = age, y = ncells_mut)) +
  geom_line(aes(x = age, y = ncells_mut)) +
  geom_hline(yintercept = 1e5 * 0.02, linetype = "dotted") +
  scale_y_continuous(limits = c(0, 1e5)) +
  theme_cowplot()

# calculate that this occurs for every individual
age_mut = expected_rate_DNMT3A_R882H |>
  select(age, mle) |>
  mutate(mle = mle, type = "mut_induction")
age_ch = age_mut |>
  mutate(age = age + 46, type = "ch-detection")

ukbiobank_ch = DNMT3A_age_frequencies |>
  mutate(mle = `R/H` / Individuals) |>
  select(Age, mle) |>
  dplyr::rename(age = Age) |>
  mutate(type = "R882H CH UKBiobank")

# plot comparison vs mutation induction and CH onset
comparison_muts_CH = rbind(age_mut, age_ch, ukbiobank_ch)

ggplot(comparison_muts_CH, aes(color = type)) +
  geom_point(aes(x = age, y = mle)) +
  geom_smooth(aes(x = age, y = mle), method = "lm", fullrange = TRUE, se = FALSE) +
  theme_cowplot() +
  scale_y_continuous( labels = percent) +
#  scale_x_continuous(limits = c(0, 80)) +
  labs(y = "percent of individuals with R882H mutation", x = "Age (years") +
  ggtitle("estimations of CH rate - occurrence delay")
comparison_muts_CH

ggsave("plots/blood/comparison_muts_CH.png", width = 6, height = 4.5, bg = "white")

# plot mutation rates in all genes
site_freqs_ch_drivers = CH_bDM[driver == TRUE, .N, by = c("mut_type", "gene_name")]
mrate_CH_drivers = get_gene_rate(exp_rates = expected_rates, metadata = metadata,
                                 site_freqs = site_freqs_ch_drivers,  ratios = ratios,  ncells = ncells)

# Get the last point of each group
df_last <- mrate_CH_drivers |>
  group_by(gene_name) |>
  filter(mle == max(mle)) |>
  arrange(mle) |>
  tail(15)

CH_driver_rates = ggplot(mrate_CH_drivers |> filter(gene_name %in% df_last$gene_name),
                         aes(x = age, y = mle, color = gene_name)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggrepel::geom_text_repel(data = df_last, aes(label = gene_name), nudge_x = 7,
                           size = 4, na.rm = TRUE, max.overlaps = Inf) +
  expand_limits(x = max(df_last$age) + 10) + # Add space for labels
  ggsci::scale_color_igv() +
  theme_cowplot() +
  theme(legend.position = "none") +
  labs(y = "number of CH driver\nmutations", x = "Age (years)")
ggsave("plots/blood/comparison_muts_CH_rates.png", CH_driver_rates,  width = 6, height = 5, bg = "white")

#  TODO Update the blood script with the new estimates for high and low numbers of mutation rate
# min number of stem cells in the blood
min_ncells = 2.5e4
max_ncells = 1.3e6

# DNMT3A mutations all
site_freqs_ch_drivers = CH_bDM[driver == TRUE & gene_name == "DNMT3A", .N, by = c("mut_type", "gene_name")]
estimates_area = get_mut_est_conf(site_freqs = site_freqs_ch_drivers, exp_rates = expected_rates,
                                  metadata = metadata, ratios = ratios, ncells = ncells,
                                  min_ncells = min_ncells, max_ncells = max_ncells)

plot_muts_area = function(estimates) {
  estimates |>
    filter(gene_name == "DNMT3A") |>
    ggplot(aes(x = age, y = ncells)) +
    geom_ribbon(aes(ymin = model_min, ymax = model_max ), alpha = 0.1) +
    geom_line(aes(y = model_ncells), linewidth = 1) +
    geom_line(aes(y = model_min), color = "grey30", linewidth = 1, linetype = "dashed") +
    geom_line(aes(y = model_max), color = "grey30", linewidth = 1, linetype = "dashed") +
    theme_cowplot()  +
    labs(y = "mutated HSCs with DNMT3A driver", x = "Age (years)")
}


label_min = paste0("25,000 HSCs\nEst. muts 80y: ", round(estimates_area[["model_min"]][8], 2))
label_mid = paste0("100,000 HSCs\nEst. muts 80y: ", round(estimates_area[["model_ncells"]][8], 2))
label_max = paste0("1.3M HSCs\nEst. muts: 80y: ", round(estimates_area[["model_max"]][8], 2))

# plot mutation rate estimates for different initial estimates - plot in the confidence interval
plot_confi = plot_muts_area(estimates_area)  +
  labs(y = "Number of cells with DNMT3A\ driver mutation") +
  annotate("text", x = 80, y = estimates_area[["model_ncells"]][8], label = label_mid, hjust = 0, vjust = 0, size = 3.5) +
  annotate("text", x = 80, y = estimates_area[["model_min"]][8], label = label_min, hjust = 0, vjust = 0.6, color = "grey30", size = 3.5) +
  annotate("text", x = 80, y = estimates_area[["model_max"]][8], label = label_max, hjust = 0, color = "grey30", size = 3.5) +
  theme(plot.margin = margin(5,32,5,5, unit = "mm")) +
  coord_cartesian(clip = "off")
plot_confi

site_freqs_ch_drivers = CH_bDM[driver == TRUE & gene_name == "DNMT3A" & aachange == "R882H",.N,by = c("mut_type", "gene_name")]
estimates_area = get_mut_est_area(site_freqs = site_freqs_ch_drivers, exp_rates = expected_rates,
                                  metadata = metadata, ratios = ratios, ncells = ncells,
                                  min_ncells = min_ncells,max_ncells = max_ncells)



# check the ages needed for
delay100k <- log(2000) / log(1.18)
delay25k <- log(500) / log(1.18) # delay is much closer

library(ggforce)
estimates_area_ch = estimates_area |>
  mutate(age = age + delay100k)
estimates_area_ch_25k = estimates_area |>
  mutate(age = age + delay25k)

plot_muts_area(estimates_area) +
  geom_pointpath(data = estimates_area_ch, aes(y = model_min), color = blood_colors ) +
  geom_point(data = ukbiobank_ch, aes(y = mle)) +
  facet_zoom2(ylim = c(0, 0.1)) +
  panel_border() +
  scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.1))) +
  scale_x_continuous(limits = c(0, 100))


# Plot for Figure 4
plot_muts_area_DNMT3A = ggplot(estimates_area, aes(x = age)) +
  geom_line(aes(y = model_ncells)) +
  geom_line(aes(y = model_min), color = "grey30", linetype = "dashed") +
  geom_line(data = estimates_area_ch, aes(y = model_ncells), color = blood_colors) +
  geom_line(data = estimates_area_ch_25k, aes(y = model_min), color = blood_colors, linetype = "dashed") +
  geom_point(data = ukbiobank_ch, aes(y = mle)) +
  scale_y_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.1))) +
  scale_x_continuous(limits = c(0, 110), breaks = c(0, 50, 100)) +
  coord_cartesian(clip = 'off') +
  theme_cowplot() +
  labs(x = "Age (years)", y = "Number of cells with\nDNMT3A R882H mutation/CH")
plot_muts_area_DNMT3A = plot_muts_area_DNMT3A  +
  annotate(geom = "textsegment", y = estimates_area$model_ncells[5], yend = estimates_area$model_ncells[5], x = 48, xend = 94, label = "46-year expansion time", size = 2.7,
           arrow = arrow(type = "closed" ,ends = "both", length = unit(0.075, "inches")))
plot_muts_area_DNMT3A

# wrap plots into a patchwork of plots
design =
  "ABC
 DDD
 EFG"
empty_plot = ggplot() + theme_void()

# make a list of the plots and add them together:
f4p = list(figure_2a = figure_2a, DNMT3A_drivers_by_sig = DNMT3A_drivers_by_sig, expected_DNMT3A = expected_DNMT3A,
           plot_confi = plot_confi, plot_muts_area_DNMT3A = plot_muts_area_DNMT3A, CH_driver_rates = CH_driver_rates)
f4p = lapply(f4p, \(x) x + theme(plot.margin = margin(5,5,5,5, "mm")))
f4p$plot_confi = f4p$plot_confi + theme(plot.margin = margin(6,37,6,7, "mm"))

figure_4 = wrap_plots(empty_plot , f4p$figure_2a , f4pa$DNMT3A_drivers_by_sig,
                      f4p$expected_DNMT3A,
                      f4p$plot_confi, f4p$plot_muts_area_DNMT3A, f4p$CH_driver_rates, nrow = 3, design = design, heights = c(1,1,1.2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave("plots/manuscript/main_figures/figure_4.svg", figure_4, width = 14, height = 12, bg = "white")