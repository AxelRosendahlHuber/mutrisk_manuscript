# prepare blood
library(data.table)
library(dndscv)
library(tidyverse)
library(MutationalPatterns)
library(cowplot)
library(wintr)
library(mutrisk)
library(R.utils)
library(GenomeInfoDb)
library(plyranges)

# [input]
tissue = "blood"
ncells = 1e5  # estimate based on the number of HSPCs actively contributing to the blood at any given moment

outdir = paste0("processed_data/", tissue, "/")

# get the annotated mutations per donor
# for the blood samples this is not directly straightforward
sample_mut_list = list.files("processed_data/blood/processed_blood_normal/", 'sample_mutations', full.names = TRUE)
names(sample_mut_list) = gsub("annotated_mut_set_|_4_01_|_2_01_|_5_01_|standard_rho01.gz_sample_mutations.txt.gz", "", basename(sample_mut_list))
cell_muts = lapply(sample_mut_list, fread) |>
  rbindlist(idcol = "donor") |>
  mutate(category = "normal")

# Check that samples KX007 and AX001 have the same mutation data
x = cell_muts |> filter(donor == "KX007") |> dplyr::select(sampleID, chr, pos, ref, alt)
y = cell_muts |> filter(donor == "AX001") |> dplyr::select(sampleID, chr, pos, ref, alt)

metadata = fread("raw_data/blood/Summary_cut.csv") |>
  dplyr::rename(sampleID = PDID, donor = donor_id) |>
  dplyr::select(donor, sampleID, age, mean_depth) |>
  mutate(category = "normal")

metadata = metadata |>
  filter(donor != "KX007") |> # Since either KX007 and AX001 have the same mutations
  mutate(vaf_estimate = 0.5, # all samples are coming from clonal cultures,
         sensitivity = get_sensitivity(coverage = mean_depth, vaf = vaf_estimate)) |>
  dplyr::rename(coverage = mean_depth)

metadata = metadata |>
  dplyr::select(sampleID, age, donor, sensitivity, coverage, category) |>
  mutate(category = factor(category, levels = c("normal", "chemotherapy")))

# double check on the metadata of Mitchell et al reveals that there are 361 samples with donor AX001 (match), and 315 samples with donor KX007 (no match)
fread("raw_data/blood/Summary_cut.csv") |>
  dplyr::rename(sampleID = PDID, donor = donor_id)  |>
  dplyr::select(sampleID, donor, age) |>
  distinct() |>
  pull(donor) |>
  table()

# plot the variant allele fraction for all donors and for one donor specific curves
# furhtermore make two plots on the VAF of the donors of this data
VAF_per_sample = cell_muts |>
  group_by(donor, sampleID) |>
  summarize(mean_vaf = mean(vaf)) |>
  ggplot(aes(x = donor, y = mean_vaf)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  ggbeeswarm::geom_quasirandom() +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_y_continuous(limits = c(0, 0.75),
                     breaks = scales::extended_breaks(4)) +
  labs(y = "average VAF per sample")

set.seed(10)
AX001_plot = cell_muts |>
  filter(donor == "AX001") |>
  filter(sampleID %in% base::sample(unique(sampleID), 40)) |>
  mutate(y = as.numeric(as.factor(sampleID))) |>
  group_by(y) |>
  ggplot(aes(x = vaf, y = y, group = y)) +
  ggridges::geom_density_ridges() +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme_cowplot() +
  labs(title = "40 clones from donor SX001", y = "HSC/MPP clones", x = "VAF of individiual mutations across sample")

set.seed(10)
SX001_plot = cell_muts |>
  filter(donor == "SX001") |>
  filter(sampleID %in% base::sample(unique(sampleID), 40)) |>
  mutate(y = as.numeric(as.factor(sampleID))) |>
  group_by(y) |>
  ggplot(aes(x = vaf, y = y, group = y)) +
  ggridges::geom_density_ridges() +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme_cowplot() +
  labs(title = "40 clones from donor SX001", y = "HSC/MPP clones", x = "VAF of individiual mutations across sample")

supplementary_note_plot_blood = VAF_per_sample / (AX001_plot | SX001_plot)
ggsave("manuscript/Supplementary_notes/Supplementary_Note_X/figure_blood.png", supplementary_note_plot_blood, width = 10, height = 7)

# thus, samples from donor KX007 can be filtered out
cell_muts = cell_muts |>
  filter(donor != "KX007")  |>
  dplyr::select(sampleID, chr, pos, ref, alt, category, donor)

# check the effect of the coverage on the mutation rate and the adjustment for sensitivity on it
# normal blood
results = effect_coverage_vaf(cell_muts |> filter(category == "normal"),
                              metadata = metadata |> filter(coverage > 9 & category == "normal")) # at a coverage of 9 we observe no effect of coverage
ggsave(paste0("manuscript/Supplementary_notes/Supplementary_Note_I/", tissue, "_coverage_vaf_correction.png"), results$plot, width = 7, height = 4.5)

# filter data for minimal coverage threshold:
metadata_filtered = metadata |> filter(coverage > 9)
cell_muts_filtered = cell_muts |> filter(sampleID %in% unique(metadata_filtered$sampleID)) |>
  filter(donor %in% unique(metadata_filtered$donor))

# save data
fwrite(cell_muts_filtered, file = paste0("processed_data/", tissue, "/", tissue, "_cell_muts.tsv.gz"))
fwrite(metadata_filtered, paste0(outdir, tissue, "_metadata.tsv"))

# make list for signatures used for re-fitting:


# signatures downloaded from:
#https://github.com/emily-mitchell/chemotherapy/ > 5_Mutational_signature_analysis >
# Current link: https://github.com/emily-mitchell/chemotherapy//blob/main/5_Mutational_signature_analysis/mutational_signatures_analysis/SBS_signatures_profiles.txt
mitchell_2025_sigs = fread("raw_data/blood/mutational_signatures_analysis/SBS_signatures_profiles.txt") |>
  column_to_rownames("Type")

mitchell_2025_sigs = mitchell_2025_sigs[mutrisk:::TRIPLETS_96, ]

input_sig_list = list(
  normal = as.matrix(mitchell_2025_sigs[, c("SBS1+SBS5", "SBSBlood")]),
  chemotherapy = as.matrix(mitchell_2025_sigs))

colnames(input_sig_list$normal) = c("SBS1SBS5", "SBSBlood")
### for now, only do the analyses using the normal blood:
metadata = metadata |> filter(category == "normal")

list_results = list()
for (i in unique(metadata$category)) {

  # test_sampleIDs = metadata_filtered |> filter(category == i) |> pull(sampleID) |> unique()
  # test_sampleIDs = test_sampleIDs[1:10]
  mutrisk_results = mutrisk_pipeline(output_path = outdir,
                                     cell_muts = cell_muts_filtered |> filter(category == i),
                                     metadata = metadata_filtered |> filter(category == i),
                                     name = i,
                                     input_signatures = input_sig_list[[i]],
                                     multiple_refit_methods = FALSE,
                                     sensitivity_correction =  TRUE)
}

# summarize the results from the tissue-specific analysis (see following code from and adapt:)
# summarize the mutation rate
sig_rate_files = list.files(paste0("processed_data/", tissue), pattern = 'sig_rate_per_sample',
                            full.names = TRUE, recursive = TRUE)
rates = lapply(sig_rate_files, fread)
names(rates) = gsub("_sig_rate_per_sample.tsv.gz", "", basename(sig_rate_files))
sig_donor_rates = rbindlist(rates, use.names = TRUE, fill = TRUE) |>
  mutate(category = factor(category, levels = levels(metadata$category)))  |>
  inner_join(metadata) |>
  group_by(donor, signature) |>
  summarize(across(contains(">"), mean)) |>
  pivot_longer(contains(">"))
fwrite(sig_donor_rates, file = paste0("processed_data/", tissue, "/", tissue, "_sig_donor_rates.tsv.gz"))

# load the mutation rates
mrates_files = list.files(paste0("processed_data/", tissue), pattern = "_patient_rates.tsv.gz", full.names = TRUE, recursive = TRUE)
mrates_files = mrates_files[!grepl("sig", mrates_files)] # exclude the signature-specific variants
rates = lapply(mrates_files, fread)
names(rates) = gsub("_patient_rates.tsv.gz", "", basename(mrates_files))
expected_rates = rbindlist(rates) |>
  mutate(category = factor(category, levels = levels(metadata$category)))  |>
  inner_join(metadata |> select(-sensitivity, -coverage, -category)) |>
  mutate(across(c(mle, cilow, cihigh), ~ . /sensitivity)) |> # correct for sensitivity
  dplyr::select(-c(donor, age, sensitivity))
fwrite(expected_rates, file = paste0("processed_data/", tissue, "/", tissue, "_expected_rates.tsv.gz"))

# load the relative mutation ratio for each sample for the APC gene:
dnds_files = list.files(paste0("processed_data/", tissue), pattern = "^.*dnds.rds", full.names = TRUE, recursive = TRUE)
dnds_files = dnds_files[!grepl("unique|exon", dnds_files)]
names(dnds_files) = names(rates)
ratios = lapply(dnds_files, \(x) {
  readRDS(x)$genemuts |>
    mutate(ratio = exp_syn_cv / exp_syn) |>
    dplyr::select(gene_name, ratio)}) |>
  rbindlist(idcol = "category") |>
  mutate(category = factor(category, levels = levels(metadata$category)))
fwrite(ratios, file = paste0("processed_data/", tissue, "/", tissue, "_mut_ratios.tsv.gz"))
