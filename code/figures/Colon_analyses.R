# Colon Figure 2 script
library(GenomicRanges)
library(rtracklayerr)
library(cowplot)
library(UpSetR)
library(mutrisk)
source("code/functions/analysis_variables.R")

tissue = "colon"
ncells = tissue_ncells_ci$mid_estimate[1]
colon_colors = c(normal = "#3ca951", IBD = "#6cc5b0", POLE = "#145220", POLD1 = "#222e24")

# load colon metadata
metadata = fread(paste0("processed_data/", tissue, "/", tissue, "_metadata.tsv")) |>   distinct() |>
  mutate(category = factor(category, levels = c("normal", "IBD", "POLD1", "POLE"))) |>
  select(-sensitivity, -coverage) |> distinct()

# load the mutation rates
expected_rates = fread(paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))
ratios = fread(paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))

# load the mutation data:
cancer_bDM = fread("processed_data/boostdm/boostdm_genie_cosmic/colon_boostDM_cancer.txt.gz")

# Get the actual values for the graphic to be correct
APC_1450_hotspot = cancer_bDM[gene_name == "APC" & aachange == "R1450*" , c("gene_name", "mut_type", "aachange", "position", "driver")]
expected_rate_APC_1450 = expected_rates |>
  inner_join(APC_1450_hotspot, by = "mut_type") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ratio * ncells)) |>
  group_by(donor, category) |>
  summarise(across(c(mle, cilow, cihigh, age), mean)) |>
  setDT()

# Values for in the manuscript: Individuals above 35 - APC 1450* mutation rate
expected_rate_APC_1450[category == "normal" & age > 35 , mle] |> mean()
setDT(expected_rate_APC_1450)[category == "normal" & age > 35, mle] |> min()
setDT(expected_rate_APC_1450)[category == "normal" & age > 35, mle] |> max()

figure_2a = ggplot(expected_rate_APC_1450, aes(x = age, y = mle, fill = category)) +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh, color = category), width = 0, show.legend = FALSE) +
  geom_point(shape = 21, size = 2.4, color = "white", stroke = 0.3) +
  scale_fill_manual(values = colon_colors) +
  scale_color_manual(values = colon_colors) +
  theme_cowplot() +
  labs(subtitle = "N cells with APC R1450*", y = "number of cells",
       x = "Age (years)", fill = NULL) +
  theme(legend.position = "inside", legend.position.inside = c(0.05,0.8 ))
figure_2a
ggsave("plots/colon/Figure_2A.png", width = 5, height = 4.5, bg = "white")


# expected rate KRAS driver
KRAS_G12_hotspot = cancer_bDM[gene_name == "KRAS" & driver == TRUE & position == 12 ,
                              c("gene_name", "mut_type", "aachange", "position", "driver")]
expected_rate_KRAS_G12 = expected_rates |>
  inner_join(KRAS_G12_hotspot, by = "mut_type") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(mle = mle * ratio * ncells) |>
  group_by(donor, category) |>
  summarise(mle = mean(mle)) |>
  setDT()

apc_counts_boostdm = cancer_bDM[gene_name == "APC", .N, by = c("gene_name", "mut_type",  "driver")]
expected_apc_muts = expected_rates |>
  left_join(apc_counts_boostdm, relationship = "many-to-many") |>
  inner_join(ratios, relationship = "many-to-many") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N *  ncells * ratio)) |>
  group_by(driver, category, sampleID) |>
  summarize(across(c(mle, cilow, cihigh), sum))  |>
  group_by(driver, category) |>
  summarize(across(c(mle, cilow, cihigh), mean))

# manuscript: expected number of APC mutations in cohort
expected_apc_muts |> filter(driver == TRUE & category == "normal")
metadata |> select(age, donor, category) |> group_by(category) |> summarize(age = mean(age))

# Expected number of cells with double mutations:
apc_double_driver = expected_rates |>
  left_join(ratios |> filter(gene_name == "APC")) |>
  left_join(metadata) |>
  inner_join(apc_counts_boostdm |> filter(driver), by = "mut_type", relationship = "many-to-many") |>
  mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio)) |>
  group_by(category, donor, age,  mut_type) |>
  summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop_last") |>
  summarize(across(c(mle, cilow, cihigh), sum)) |>
  mutate(across(c(mle, cilow, cihigh), ~ ((.^2) / 2) * ncells ))

apc_double_driver_mean = apc_double_driver |>
  group_by(category) |>
  summarize(max = max(mle)*1.1,
            mle = mean(mle))

plot_double_all = apc_double_driver |>
  ggplot(aes(x = age, y = mle, color = category, fill = category)) +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh), linewidth = 1, width = 0) +
  geom_point(shape = 21, color = "white", size = 3, stroke = 0.5) +
  theme_cowplot() +
  scale_color_manual(values = colon_colors) +
  scale_fill_manual(values = colon_colors) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  theme(legend.position = "inside", legend.position.inside = c(0.03, 0.7), legend.key.size = unit(5, "mm"),
        legend.title = element_blank()) +
  panel_border() +
  labs(x = "Age (years)", y = "number of cells with\ndouble mutation") +
  ggforce::facet_zoom(ylim = c(0, 15), zoom.size = 1)
plot_double_all
ggsave("plots/colon/number_of_double_mutated_cells_zoom.png", plot_double_all, width = 8, height = 4.5, bg = "white")

#### profiles:
# load genie data:
genie_colorectal = fread("processed_data/GENIE_17/GENIE_17_processed.txt.gz") |>
  filter(ONCOTREE_CODE %in%  c("COAD", "READ", "COADREAD"))
label_df_genie = data.frame(label = "GENIE APC:\ncolorectal cancer data")

apc_colon_count = genie_colorectal[gene_name == "APC", .N, by = "position"]

# profile of the driver rates for APC:
apc_muts = expected_rates |>
  group_by(category, mut_type) |>
  summarize(across(c(mle, cilow, cihigh), mean)) |>
  left_join(ratios[gene_name == "APC"]) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ncells * ratio)) |>
  left_join(triplet_match_substmodel) |>
  left_join(cancer_bDM[gene_name == "APC", ], by = c("mut_type", "gene_name"), relationship = "many-to-many") |>
  mutate(driver_status = ifelse(driver, "driver", "non-driver"))  |> setDT()

apc_count = apc_muts[, .(mle = sum(mle)) , by = c("position", "category", "driver_status") ]
normal_apc = apc_muts[category == "normal",]
normal_apc_count = apc_count[category == "normal",]

label_df = data.frame(label = "expected number of crypt stem cells mutated",
                      driver_status = "driver")

expected_apc = ggplot(normal_apc, aes(x = position, y = mle)) +
  geom_col(aes(fill = type)) +
  geom_point(data = normal_apc_count) +
  ggpp::geom_text_npc(data = label_df,  aes(label = label),  npcx = 0.05, npcy = 0.96) +
  scale_fill_manual(values = COLORS6) +
  labs(y = "numer of mutated cells") +
  facet_grid(driver_status ~ .) +
  theme_cowplot() + panel_border() +
  scale_y_continuous(expand=expansion(mult=c(0,0.1)), labels = label_comma())

APC_observed_expected = genie_apc /  expected_apc + plot_layout(guides = "collect", heights = c(1, 1.5))
ggsave("plots/colon/APC_observed_expected.png", APC_observed_expected,  width = 15, height = 6)
ggsave("plots/colon/APC_expected.png", expected_apc ,  width = 15, height = 5.5, bg = "white")


#### TP53 mutations for colon:
# profile of the driver rates for colon:
TP53_muts = expected_rates |>
  group_by(category, mut_type) |>
  summarize(across(c(mle, cilow, cihigh), mean)) |>
  left_join(ratios[gene_name == "TP53"]) |>
  mutate(across(c(mle, cilow, cihigh), ~ . * ncells * ratio)) |>
  left_join(triplet_match_substmodel) |>
  left_join(cancer_bDM[gene_name == "TP53", ], by = c("mut_type", "gene_name"), relationship = "many-to-many") |>
  mutate(driver_status = ifelse(driver, "driver", "non-driver"))  |> setDT()

TP53_count = TP53_muts[, .(mle = sum(mle)) , by = c("position", "category", "driver_status") ]
normal_TP53 = TP53_muts[category == "normal",]
normal_TP53_count = TP53_count[category == "normal",]

label_df = data.frame(label = "TP53: expected number of crypt stem cells mutated",
                      driver_status = "driver")

expected_TP53 = ggplot(normal_TP53, aes(x = position, y = mle)) +
  geom_col(aes(fill = type)) +
  geom_point(data = normal_TP53_count) +
  ggpp::geom_text_npc(data = label_df,  aes(label = label),  npcx = 0.05, npcy = 0.96) +
  scale_fill_manual(values = COLORS6) +
  labs(y = "numer of mutated cells") +
  facet_grid(driver_status ~ .) +
  theme_cowplot() + panel_border() +
  scale_y_continuous(expand=expansion(mult=c(0,0.1)), labels = label_comma())
ggsave("plots/colon/TP53_expected.png", expected_TP53 ,  width = 15, height = 5.5, bg = "white")

#####
# get the percentage of mutated cells for APC
####
site_freqs_APC_drivers = cancer_bDM[gene_name == "APC" & driver == TRUE,  .N, by = c("mut_type", "gene_name")]

# percentage all drivers
all_drivers = cancer_bDM[driver == TRUE,  .N, by = c("mut_type", "gene_name")]
driver_rates  = get_gene_rate(exp_rates = expected_rates, metadata = metadata, site_freqs = all_drivers, ratios = ratios, ncells = 1)
mean_mutated = driver_rates |> group_by(category, gene_name) |>
  group_by(category, gene_name) |>
  filter(age > 18) |>
  summarize(across(c(mle, cilow, cihigh), mean)) |>
  arrange(mle) |>
  mutate(gene_name = factor(gene_name, unique(gene_name)))

ggplot(mean_mutated, aes(x = fct_reorder(gene_name, mle, .desc = TRUE), y = mle)) +
  geom_col() +
  facet_grid(category ~ ., scales = "free_y") +
  labs(x = NULL, y = "percent of cells carrying a mutation") +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.1))) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# adress comment ferran - inclusion exclusion:
probs = expected_rates |>
  left_join(distinct(metadata), by = c("category", "sampleID")) |>
  inner_join(all_drivers, by = "mut_type", relationship = "many-to-many") |>
  left_join(ratios, by = c("category", "gene_name")) |>
  group_by(category, mut_type, donor, age, ratio, gene_name, N) |>
  summarize(across(c(mle, cilow, cihigh),  ~ mean(. * ratio))) |>
  filter(!is.na(donor)) |>
  as.data.table()

Rcpp::sourceCpp("code/functions/inclusion_exclusion2.cpp")

fraction_mut_test = tibble(donor = unique(probs$donor), `inclusion exclusion` = NA, sum = NA)
for (select_donor in unique(probs$donor)) {
  print(select_donor)
  probs_donor = probs[donor %in% select_donor, ]
  seq = rep(probs_donor$mle, probs_donor$N)
  fraction_mut_test[fraction_mut_test$donor == select_donor, 2 ]  = inclusion_exclusion2(seq)
  fraction_mut_test[fraction_mut_test$donor == select_donor, 3 ] =  sum(seq)
}

fraction_mut_test = fraction_mut_test |>
  mutate(`fraction multiple muts` = (sum - `inclusion exclusion`) / sum) |>
  left_join(metadata |> select(donor, category) |> distinct())

barplot_ie_vs_sum = fraction_mut_test |> pivot_longer(c(`inclusion exclusion`, sum), values_to = "percent mutated", names_to = "method") |>
  ggplot(aes(x = donor, y = `percent mutated`, alpha = method,  fill = category)) +
  geom_col(position = "dodge", color = "black") +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = colon_colors) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.1))) +
  theme_cowplot()  +
  theme(axis.text.x = element_blank(), legend.position = "inside", legend.position.inside = c(0.05, 0.75))  +
  labs()

percent_double_mut = ggplot(fraction_mut_test, aes(x = donor, y = `fraction multiple muts`, fill = category)) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = colon_colors) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.1))) +
  theme_cowplot()  +
  theme(axis.text.x = element_blank(), legend.position = "none")

percent_non_double_mut =  ggplot(fraction_mut_test, aes(x = donor, y = 1 - `fraction multiple muts`, fill = category)) +
  geom_col() +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = colon_colors) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.1))) +
  theme_cowplot()  +
  theme(axis.text.x = element_blank(), legend.position = "none")

mg = 7
supplementary_note_II_figure = barplot_ie_vs_sum / (percent_double_mut | percent_non_double_mut )  +
  plot_annotation(tag_levels  = "A") &
  theme(plot.margin = margin(mg, mg, mg, mg, unit = "mm"))
supplementary_note_II_figure
ggsave("plots/manuscript/supplementary_notes/Supplementary_Note_II_incl_excl.png", supplementary_note_II_figure, width = 9, height = 9)

# average across clones
site_freqs_APC_drivers = cancer_bDM[gene_name == "APC" & driver == TRUE,  .N, by = c("mut_type", "gene_name")]
mrate_APC_drivers = get_gene_rate(exp_rates = expected_rates, metadata = metadata,
                                  site_freqs = site_freqs_APC_drivers, ratios = ratios, ncells = ncells)
plot_APC_drivers_all = plot_mrate(mrate_APC_drivers, title = "number of cells with APC driver mutations", colors = colon_colors)
ggsave("plots/colon/ncells_APC_drivers.png", plot_APC_drivers_all, width = 5.5, height = 4.2, bg = "white")

# double check - the differences between the mutation rates are so small that this indicates rounding errors
#### APC double driver mutations
mrate_APC_double_drivers = get_double_gene_rate(exp_rates = expected_rates, metadata = metadata,
                                                site_freqs = site_freqs_APC_drivers, ratios = ratios, ncells = ncells)
plot_apc_double_drivers_all = plot_mrate(mrate_APC_double_drivers, title = "number of cells with double APC driver mutations - standard method", colors = colon_colors) +
  ggforce::facet_zoom(ylim = c(0, 10), zoom.size = 1)
ggsave("plots/colon/ncells_APC_double_drivers_all_points.png", plot_apc_double_drivers_all,  width = 10, height = 5, bg = "white")

plot_apc_double_drivers_normal = plot_mrate(mrate_APC_double_drivers |> filter(category == "normal"), title = "number of cells with APC driver mutations", colors = colon_colors) +
  labs(subtitle = "normal tissue - using all cells averaged for each donor")
ggsave("plots/colon/ncells_APC_double_drivers_normal_points.png", plot_apc_double_drivers_normal,  width = 6, height = 5, bg = "white")

# Manuscript: make general overview table of the different conditions. This can be a supplementary or main figure table
mrate_APC_double_drivers_fraction |>
  filter(category == "normal" & age > 35) |>
  pull(mle) |> unique() |> summary()

double_drivers = mrate_APC_double_drivers_fraction |>
  filter(age > 35) |>
  group_by(category) |>
  summarize(`mean exp. \ndouble APC driver` = paste0(round(mean(mle),1), " (", round(min(mle), 1), "-", round(max(mle), 1), ")"),
            `number of donors` = dplyr::n())

single_drivers = mrate_APC_drivers |>  filter(age > 35) |>
  group_by(category) |>
  summarize(`mean age\nin years\n(min-max)` = paste0(round(mean(age),1), " (", round(min(age), 1), "-", round(max(age), 1), ")"),
            `mean cells exp. \nAPC driver (min-max)` = paste0(
              format(round(mean(mle)), big.mark = ","), " (",
              format(round(min(mle)), big.mark = ","), "-",
              format(round(max(mle)), big.mark = ","), ")"))

metadata_counts = metadata |> group_by(category) |>
  summarize(ncrypts = dplyr::n())

table_2 = left_join(single_drivers, double_drivers) |>
  left_join(metadata_counts)
writexl::write_xlsx(table_2, "manuscript/Table_2.xlsx")

# Supplementary Figure 5 manuscript : get the confidence interval of the mean:
mrate_APC_drivers |>  filter(age > 35) |>
  group_by(category) |>
  filter(category == "normal") |>
  summarize(mean = mean(mle),
            low = mean(cilow),
            high = mean(cihigh)) |>
  as.data.frame()

rbindlist(list(`diversity clones` = mrate_APC_double_drivers_fraction |> filter(category == "normal", ),
               `mean across donor` = mrate_APC_double_drivers |> filter(category == "normal", )), idcol = "name") |>
  ggplot(aes(x = fct_reorder(donor, mle), y = mle, color = name)) +
  geom_point(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = cilow, ymax = cihigh), position = position_dodge(width = 0.9)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "inside", legend.position.inside = c(0.02, 0.9)) +
  labs(y = "predicted number of cells with double mutation", x = NULL, color = "computing method",
       subtitle = "Comparison calculation methods # cells with two APC driver mutations")
ggsave("plots/colon/Supplementary_figure_5.png", width = 7, height = 6, bg = 'white')

# all TP53 driver muts:
site_freqs_TP53_drivers = cancer_bDM[gene_name == "TP53" & driver == TRUE,  .N, by = c("mut_type", "gene_name")]
mrate_TP53_drivers = get_gene_rate(exp_rates = expected_rates, metadata = metadata,
                                   site_freqs = site_freqs_TP53_drivers, ratios = ratios, ncells = ncells)
plot_TP53_drivers = plot_mrate(mrate_TP53_drivers, title = "number of cells with TP53 driver mutations", colors = colon_colors) +
  scale_y_continuous(sec.axis = sec_axis(name = "% of cells with mutation", ~ . / ncells, labels = label_percent()),
                     labels = label_comma())
ggsave("plots/colon/ncells_TP53_drivers.png", plot_TP53_drivers, width = 5.5, height = 4.2, bg = "white")

plot_mrate(mrate_TP53_drivers |> filter(category == "normal"), title = "number of cells with TP53 driver mutations", colors = colon_colors) +
  labs(subtitle = "normal tissue") +
  scale_y_continuous(sec.axis = sec_axis(name = "% of cells with mutation", ~ . / ncells, labels = label_percent()),
                     labels = label_comma())
ggsave("plots/colon/ncells_TP53_drivers_normal.png", width = 5.5, height = 4.2, bg = "white")

###
# signature-specific modeling of APC mutation rates
###

# make function to mask all other signature contributions:
mask_sigs = function(mrate_sigs, select_sigs = c("SBS1", "SBS5", "SBS18", "SBS88"), value_columns = c("mle", "cilow", "cihigh")) {
  mrate_sigs |>
    mutate(signature = case_when(!signature %in% c(select_sigs) ~ "other", .default = signature)) |>
    group_by(across(-any_of(value_columns))) |>
    summarize(across(any_of(value_columns), sum), .groups = "drop")
}

# Load signature specific mutation rates
signature_rates = fread("processed_data/colon/colon_sig_patient_rates.tsv.gz") |>
  dplyr::rename(mut_type = name, mle = value) |>
  left_join(metadata, relationship = "many-to-many") |>
  mutate(category = factor(category, levels = levels(metadata$category))) |>
  mask_sigs(select_sigs = c("SBS1", "SBS10a", "SBS10b", "SBS10c", "SBS10d", "SBS13", "SBS2", "SBS5", "SBS88", "SBS89", "SBS18"))

# mrate APC all muts:
site_freqs_APC = cancer_bDM[gene_name == "APC" & driver == TRUE, .N, by = c("mut_type", "gene_name")]
mrate_sigs_APC_driver = get_gene_rate_sig(exp_rates = signature_rates, metadata,
                                          site_freqs = site_freqs_APC, ratios = ratios)

plot_signature_rates_APC_all_muts = mrate_sigs_APC_driver |>
  mutate(mle = mle * ncells) |>
  group_by(age, signature, category) |>
  summarize(mle = mean(mle)) |>
  ggplot(aes(x = age, y = mle, fill = signature)) +
  geom_col()  +
  facet_grid(. ~ category) +
  scale_fill_manual(values = sig_colors) +
  theme_cowplot() +
  panel_border()  +
  scale_y_continuous(expand=expansion(mult=c(0,0.1)), labels = label_comma(),
                     sec.axis = sec_axis(name = "% of cells with mutation", ~ . / ncells, labels = scales::label_percent())) +
  theme(legend.position = "inside", legend.position.inside = c(0.05, 0.65), legend.key.size = unit(4, "mm"), legend.title = element_blank())+
  labs(y = "number of mutated cells", x = "Age (years)", subtitle = "APC mutated cells - all muts")
plot_signature_rates_APC_all_muts
ggsave("plots/colon/number_of_cells_mutated_by_sig_APC.png", plot_signature_rates_APC_all_muts, width = 10, height = 4.2, bg = "white")

# mrate APC driver muts
mrate_sigs_APC_driver = get_gene_rate_sig(exp_rates = signature_rates, metadata, site_freqs = site_freqs_APC_drivers, ratios = ratios)
plot_signature_rates_APC_drivers = mrate_sigs_APC_driver |>
  mutate(mle = mle * ncells) |>
  group_by(age, signature, category) |>
  summarize(mle = mean(mle)) |>
  ggplot(aes(x = age, y = mle, fill = signature)) +
  geom_col()  +
  facet_grid(. ~ category) +
  scale_fill_manual(values = sig_colors) +
  theme_cowplot() +
  panel_border()  +
  scale_y_continuous(expand=expansion(mult = c(0,0.1)), labels = label_comma()) +
  theme(legend.position = "inside", legend.position.inside = c(0.05, 0.65), legend.key.size = unit(4, "mm"), legend.title = element_blank()) +
  labs(y = "number of mutated cells", x = "Age (years)", subtitle = "Cells with APC driver muatations")
ggsave("plots/colon/number_of_cells_mutated_by_sig_driver.png", plot_signature_rates_APC_drivers,  width = 10, height = 4.2, bg = "white")

# get hotspot mutation levels for colon
APC_1450_hotspot_sig_rates = APC_1450_hotspot |>
  left_join(signature_rates, by = "mut_type") |>
  left_join(metadata) |>
  group_by(age, signature, category) |>
  summarize(mle = mean(mle))

normal_APC_1450_rates = APC_1450_hotspot_sig_rates |>
  filter(category == "normal") |>
  mask_sigs()

plot_APC_1450_hotspot = ggplot(normal_APC_1450_rates, aes(x = age, y = mle * ncells, fill = signature)) +
  geom_col() +
  theme_cowplot() +
  theme(legend.position = "inside", legend.position.inside = c(0.05, 0.8)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),
                     sec.axis = sec_axis(name = "% of cells with mutation", ~ . / ncells, labels = scales::label_percent()))  +
  scale_fill_manual(values = sig_colors) +
  labs(subtitle = "Colonic crypts with R1450* APC hotspot mutations", y = "Number of cells", x = "Age (years)")
ggsave("plots/colon/R1450_stop_rates_normal.png", plot_APC_1450_hotspot, width = 5, height = 4, bg = "white")

plot_APC_1450_hotspot_all = ggplot(APC_1450_hotspot_sig_rates, aes(x = age, y = mle * ncells, fill = signature)) +
  geom_col() +
  facet_grid(. ~ category) +
  theme_cowplot() +
  scale_y_continuous(expand=expansion(mult=c(0,0.1)))  +
  scale_fill_manual(values = sig_colors) +
  labs(subtitle = "Colonic crypts with R1450* APC hotspot mutations", y = "Number of cells", x = "Age (years)")
ggsave("plots/colon/R1450_stop_rates_all.png", plot_APC_1450_hotspot_all,  width = 10, height = 4, bg = "white")

####
# mutation rates for individual clones by donor
####
apc_total = cancer_bDM[gene_name == "APC",.N, by = c("mut_type", "gene_name", "consequence")]
consequence_rates = apc_total |>
  left_join(signature_rates |>
              filter(category == "normal"), relationship = "many-to-many", by = "mut_type") |>
  left_join(ratios)  |>
  mutate(mle = mle * N * ncells * ratio) |>
  group_by(consequence, sampleID, signature) |>
  summarize(mle = sum(mle)) |>
  left_join(metadata)

# for manuscript: average number of consequence type in normal tissues
mean_consequence_rates = consequence_rates |>
  left_join(metadata) |>
  filter(age > 35) |>
  group_by(consequence, age, sampleID) |>
  summarize(mle = sum(mle), .groups = "drop_last") |>
  summarize(mle = mean(mle), .groups = "drop_last") |>
  summarize(mean_mle = mean(mle), min = min(mle), max = max(mle))

# manuscript: the outlier with high SBS88 profile
consequence_rates |>
  left_join(metadata) |>
  filter(age > 34) |>
  filter(category == "normal" & consequence == "splicing") |>
  group_by(donor, age,  sampleID) |>
  summarize(mle = sum(mle), .groups = "drop_last") |>
  summarize(mle = mean(mle))  |> arrange(desc(mle))

consequence_rates |>
  left_join(metadata) |>
  filter(age > 34) |>
  filter(category == "normal" & consequence == "splicing")  |>
  ungroup() |>
  select(donor, age) |> distinct() |>
  filter(age == 36)

metadata |> select(donor, age, category) |> distinct() |>
  filter(age == 36)


# apc mut types drivers
apc_total_drivers = cancer_bDM[gene_name == "APC" & driver == TRUE, .N, by = c("mut_type", "gene_name", "driver")]
driver_rates = apc_total_drivers |>
  left_join(signature_rates, relationship = "many-to-many", by = "mut_type") |>
  left_join(ratios[gene_name == "APC",])  |>
  mutate(mle = mle * N * ncells * ratio) |>
  group_by(driver, sampleID, signature) |>
  summarize(mle = sum(mle)) |>
  arrange(driver, sampleID, signature)

driver_muts = driver_rates |> left_join(metadata) |>
  group_by(category, age, signature, driver) |>
  summarize(mle = mean(mle))

# calculate the mutation rate for all drivers
# get the double rates, in a fraction-based approach
# first, all signatures
APC_sig_rates = get_gene_rate_sig(exp_rates = signature_rates, metadata = metadata,
                                  site_freqs = site_freqs_APC_drivers, ratios = ratios)
SBS88_percentage = APC_sig_rates |>
  group_by(donor, signature) |>
  summarize(mle = mean(mle)) |>
  mutate(mle, mle /sum(mle))  |>
  filter(signature == "SBS88")

expected_rates = signature_rates |>
  group_by(category, sampleID, mut_type) |>
  summarise(mle = sum(mle)) |>
  arrange(category, mut_type, sampleID)

drivers2_all = get_double_gene_rate_fraction(exp_rates = expected_rates, metadata = metadata,
                                             site_freqs = site_freqs_APC_drivers, ratios = ratios, ncells = ncells)

expected_rates_nocb = signature_rates |>
  filter(signature != "SBS88") |>
  group_by(category, sampleID, mut_type) |>
  summarise(mle = sum(mle)) |>
  arrange(category, mut_type, sampleID)
drivers2_nocb = get_double_gene_rate_fraction(exp_rates = expected_rates_nocb, metadata = metadata,
                                              site_freqs = site_freqs_APC_drivers, ratios = ratios, ncells = ncells)

expected_cb_muts = list(all_sigs = drivers2_all, no_colibactin = drivers2_nocb) |>
  rbindlist(idcol = "type") |>
  filter(category == "normal") |>
  mutate(type = fct_relevel(type,  "no_colibactin", "all_sigs")) |>
  arrange(type)

# mutation rates for other driver genes
site_freqs_drivers_all = cancer_bDM[driver == TRUE, .N, by = c("mut_type", "gene_name")]
mrate_sigs_all_driver = get_gene_rate_sig(exp_rates = signature_rates, metadata,
                                          site_freqs = site_freqs_drivers_all, ratios = ratios)

plot_all_drivers = mrate_sigs_all_driver |>
  group_by(category, age, gene_name, donor) |>
  summarize(mle = sum(mle) * ncells) |>
  ggplot(aes(x = age, y = mle, color = gene_name)) +
  geom_line() +
  geom_point() +
  facet_grid(. ~ category) +
  ggsci::scale_color_igv() +
  scale_y_continuous(sec.axis = sec_axis(name = "% of cells with mutation", ~ . / ncells, labels = scales::label_percent()),
                     labels = label_comma()) +
  theme_cowplot() + panel_border() +
  labs(y = "number of mutated cells for gene" )
plot_all_drivers
ggsave(paste0("plots/", tissue, "/all_driver_genes.png"), width = 10, height = 4, bg = "white")

# get the average rates for the individual cells
gene_muts = cancer_bDM[gene_name == "APC" , c("chr", "pos", "gene_name", "mut_type", "aachange", "position")]
expected_rate_pos_sc = expected_rates |>
  inner_join(gene_muts, by = "mut_type", relationship = "many-to-many") |>
  left_join(metadata) |>
  left_join(ratios) |>
  mutate(mle = mle* ratio)

# Make the full plots using patchwork
empty_plot = ggplot() + cowplot::theme_nothing()
plot_list = list(empty_plot = empty_plot,
                 plot_APC_1450_hotspot = plot_APC_1450_hotspot,
                 plot_consequence_APC = plot_consequence_APC,
                 plot_muts_driver_normal = plot_muts_driver_normal,
                 plot_APC_drivers_normal = plot_APC_drivers_normal,
                 plot_apc_double_drivers_normal = plot_apc_double_drivers_normal,
                 rel_difference_plot_colibactin = rel_difference_plot_colibactin,
                 plot_signature_rates_APC_drivers = plot_signature_rates_APC_drivers,
                 plot_double_all = plot_double_all,
                 plot_TP53_drivers = plot_TP53_drivers)

mg = 10
pl = lapply(plot_list, "+",  theme(plot.margin = margin(mg, mg, mg, mg, unit = "mm")))

row1 = wrap_plots(pl$empty_plot, pl$plot_APC_1450_hotspot, pl$plot_consequence_APC, nrow = 1, widths = c(1,1,1.5))
row2 = wrap_plots(pl$plot_muts_driver_normal, pl$plot_APC_drivers_normal, pl$plot_apc_double_drivers_normal, nrow = 1)
row3 = wrap_plots(pl$rel_difference_plot_colibactin, pl$plot_signature_rates_APC_drivers, nrow = 1, widths =  c(1, 2))
row4 = wrap_plots(pl$plot_double_all, pl$plot_TP53_drivers,  nrow = 1, widths = c(2, 1))

figure_2 = wrap_plots(row1, row2, row3, row4, ncol = 1) + plot_annotation(tag_levels = "A")
ggsave("plots/manuscript/main_figures/figure_2.svg", figure_2, width = 22, height = 19, bg = "white", dpi = 600)
ggsave("plots/manuscript/main_figures/figure_2.png", figure_2, width = 19, height = 20, bg = "white", dpi = 600)

figure_2 # plot the final figure

#rmarkdown::render("Figure_2_colon.R") # comment because command returns infinite loop

#### Prepare information for figure 5: Colon specific sequence-of events cancer driverness - risk
# load cancer incidence rates
# bowel_cancer = readxl::read_xlsx("raw_data/incidence_rates/cases_crude_mf_bowel_i19.xlsx", skip = 4)[1:19,] |>
#   mutate(start_age = as.numeric(str_split_i(`Age Range`, " ", 1)),
#          end_age = as.numeric(str_split_i(`Age Range`, " ", 3)),
#          `Age Range` = as.factor(`Age Range`),
#          label = str_pad(`Male Rates`, width = max(str_width(`Male Rates`)), side = "both"))  # Center pad
#
# bowel_cancer[19, "start_age"] = 90
# bowel_cancer[19, "end_age"] = 94
#
# saveRDS(bowel_cancer, "processed_data/colon/bowel_cancer.rds")

# APC rates for a single cell
apc_double_driver_sc = apc_double_driver |>
  mutate(across(c(mle, cilow, cihigh), ~ ./ ncells)) |>
  arrange(category, donor, age)

# make a dataframe with the relative mutation rates for all the different changes:
vogelgram_mut_risks_sc = tibble(
  category = apc_double_driver_sc$category,
  donor = apc_double_driver_sc$donor,
  age = apc_double_driver_sc$age,
  `APC 2x SNV` = apc_double_driver_sc$mle,
  KRAS_driver = expected_rate_KRAS_driver_sc$mle) |>
  mutate(`APC 2x +KRAS` = `APC 2x SNV` * KRAS_driver)

vogelgram_mut_risks = vogelgram_mut_risks_sc |>
  mutate(across(-c(category, donor, age), ~ . * ncells))

vogelgram_mut_risks_long = vogelgram_mut_risks |>
  filter(category == "normal") |>
  pivot_longer(-c(category, donor, age), names_to = "mutation", values_to = "expected_mutated_cells")

saveRDS(vogelgram_mut_risks_long, "processed_data/colon/vogelgram_plot.rds")




