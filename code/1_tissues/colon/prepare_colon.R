# Load metadata, mutations, and signature contributions in a common format
# format mutations: cell_muts = 7 columns: sampleID, chr, pos, ref, alt, category, donor
# format metadata: metadata = sample, age, donor, category
source("code/0_functions/analysis_variables.R")
library(wintr)
library(dndscv)
library(mutrisk)

tissue = "colon"

metadata_IBD_normal = fread("raw_data/colon/normal/cell_11495_mmc2.txt")
metadata_IBD_normal <- metadata_IBD_normal[!duplicated(metadata_IBD_normal$crypt_ID) &
                                             !duplicated(metadata_IBD_normal$crypt_ID, fromLast = TRUE),]
# remove duplicated crypt id
metadata_POL = readxl::read_excel("raw_data/colon/hypermutated/Figures/input_files/Extended_Data_Table2.xlsx")

# load data from normal and IBD colon samples
# NOTE: data from normal and IBD colon samples is branch-specific, instead of colon specific
normal_muts = fread("raw_data/colon/normal/All_control_cohort_mutations_mapped_to_branches.txt") |>
  mutate(category = "normal")
IBD_muts = fread("raw_data/colon/normal/All_IBD_cohort_mutations_mapped_to_branches.txt") |>
  dplyr::select(-Patient_ID) |>
  mutate(category = "IBD")

WT_muts = rbind(normal_muts, IBD_muts) |>
  dplyr::rename(sampleID = SampleID, chr = Chr, pos = Pos, ref = Ref, alt = Alt)

# calculate the sensitivity for the individual clones (taking VAF and depth into account)
patient_files = list.files("raw_data/colon/normal/binary_genotype_matrices/",
                           pattern = "binary_matrix.txt", full.names = TRUE)
names(patient_files) = gsub("_.*", "", basename(patient_files))
total_depth_files = list.files("raw_data/colon/normal/Pile-ups_read_counts/",
                               pattern = "total_depth", full.names = TRUE)
names(total_depth_files) = gsub("_.*", "", basename(total_depth_files))
alt_allele_files = list.files(
  "raw_data/colon/normal/Pile-ups_read_counts/",
  pattern = "alt_allele_reads",
  full.names = TRUE)

names(alt_allele_files) = gsub("_.*", "", basename(alt_allele_files))

excl_list = vector("character")
list_patient_muts = list()
for (name in names(total_depth_files)) {
  mut_calls = read.table(patient_files[[name]], sep = " ") |>
    rownames_to_column("mutID") |>
    pivot_longer(-mutID, names_to = "sampleID", values_to = "mut_call") |>
    distinct()
  total_depth = vroom::vroom(total_depth_files[[name]], show_col_types = FALSE) |>
    pivot_longer(-mutID, names_to = "sampleID", values_to = "depth") |>
    distinct() |>
    mutate(mt = paste0(sampleID, mutID))
  alt_alleles = vroom::vroom(alt_allele_files[[name]], show_col_types = FALSE) |>
    pivot_longer(-mutID, names_to = "sampleID", values_to = "alt_depth") |>
    distinct() |>
    mutate(mt = paste0(sampleID, mutID))


  if (all(total_depth$mutID == alt_alleles$mutID)) {
    total_alt_depth = total_depth |>
      mutate(alt_depth = alt_alleles$alt_depth)
  } else {
    total_alt_depth = total_depth |>
      inner_join(alt_alleles |> select(-mutID, -sampleID), by = "mt")

  }

  total_depth_id = total_alt_depth |>
    mutate(vaf =  alt_depth / depth) |>
    filter(alt_depth >= 3 & depth >= 5 & vaf >= 0.05)

  fold_diff = nrow(total_depth_id) / sum(mut_calls$mut_call == 1)
  print(fold_diff)
  if (fold_diff  > 2) {
    warning(paste0("More than 2x more muts expected mutated than mutation calls for sample: ", name))
    next
  }

  # join dataframes and filter out all sites in which the mutation call is not 1.
  donor_muts = mut_calls |>
    inner_join(total_depth_id, by = c("mutID", "sampleID")) |>
    filter(mut_call == 1)

  # calculate the VAF for each of the patients, and split the mutation ID in the respective columns
  list_patient_muts[[name]] = donor_muts |>
    mutate(VAF = alt_depth / depth) |>
    mutate(chr = str_split_i(mutID, pattern = "_", 1),
           pos = str_split_i(mutID, pattern = "_", 2),
           ref = str_split_i(mutID, pattern = "_", 3),
           alt = str_split_i(mutID, pattern = "_", 4)) |>
    dplyr::select(sampleID, chr, pos, ref, alt, VAF)
}

cell_muts = rbindlist(list_patient_muts)

# get VAFS:
vaf_estimates_normal_IBD = cell_muts |>
  group_by(sampleID) |>
  dplyr::summarize(vaf_estimate = estimate_vaf(VAF))

IBD_ids = metadata_IBD_normal |>
  dplyr::rename(sampleID = crypt_ID) |>
  mutate(category = ifelse(cohort == "control_data", "normal", "IBD")) |>
  dplyr::select(sampleID, category) |>
  distinct()

WT_muts = left_join(cell_muts, IBD_ids, by = "sampleID") |>
  filter(!is.na(category))

# load data from POLE POLD1 study (Robinson et al., Nature Genetics 2021, https://doi.org/10.1038/s41588-021-00930-y)
POL_muts = fread("raw_data/colon/hypermutated/DNAPolymerase_NG_somatic_SBS_ID_combined.txt")
metadata_crypts = metadata_POL |>
  filter(sample_type == "intestinal crypt") |>
  mutate(category = gsub(" .*", "", germline_mutation)) |>
  dplyr::select(sample, category)
POLE_POLD1_muts = inner_join(POL_muts, metadata_crypts, by = "sample") |>
  dplyr::rename(sampleID = sample)

cell_muts = rbind(WT_muts |> select(-VAF), POLE_POLD1_muts)

# save the metadata for the different samples:
# PART 2: Modeling of mutation rates for the individual samples
meta_IBD_normal = metadata_IBD_normal |>
  dplyr::mutate(category = ifelse(cohort == "control_data", "normal", "IBD")) |>
  dplyr::select(crypt_ID, Age, patient_ID, category, Coverage) |>
  dplyr::rename(sampleID = crypt_ID,
                age = Age,
                donor = patient_ID,
                coverage = Coverage)

# Update the metadata information parts:
meta_IBD_normal = meta_IBD_normal |>
  left_join(vaf_estimates_normal_IBD) |>
  mutate(sensitivity = get_sensitivity(coverage, vaf_estimate)) |>
  select(sampleID, category, donor, age, sensitivity, coverage) |>
  filter(!is.na(sensitivity)) # TODO - some patients are left out. Check why this is happening.

meta_POLE_POLD1 = metadata_POL |>
  filter(sample_type == "intestinal crypt") |>
  dplyr::rename(category = germline_mutation,
                sampleID = sample,
                donor = individual,
                vaf = median_vaf) |>
  dplyr::select(sampleID, age, donor, category, vaf, coverage) |>
  mutate(coverage = as.numeric(coverage), vaf = as.numeric(vaf),
         category = case_when(grepl("POLE", category) ~ "POLE",
                              grepl("POLD1", category) ~ "POLD1"))

meta_POLE_POLD1 = meta_POLE_POLD1 |>
  mutate(sensitivity = get_sensitivity(coverage, vaf, minalt = 4)) |>
  select(sampleID, category, donor, age, sensitivity, coverage)

metadata = rbind(meta_POLE_POLD1, meta_IBD_normal) |>
  mutate(coverage = as.numeric(coverage),
         category = factor(category, levels = c("normal", "IBD", "POLD1", "POLE")))

# to make the mutation rate model: merge the total mutations with the donor ID to remove shared ancestral mutations from the list during dndscv. (as these would otherwise be counted double)
cell_muts = cell_muts |>
  left_join(distinct(metadata)) |>
  select(-age)

# check the effect of the coverage on the mutation rate and the adjustment for sensitivity on it
# with a minimal coverage of 11 we can correct for all of the samples
metadata_filtered = metadata |> filter(coverage > 11)

list_results = list()
for (i in unique(metadata$category)) {
  list_results[[i]] = effect_coverage_vaf(cell_muts |> filter(category %in% i),
                      metadata = metadata_filtered |> filter(category %in% i))
  list_results[[i]]$plot = list_results[[i]]$plot + ggtitle(i)
}

list_plots = lapply(list_results, \(x) x[["plot"]])
ggsave(paste0("manuscript/Supplementary_notes/Supplementary_Note_I/", tissue, "_coverage_var_correction.png"),
       patchwork::wrap_plots(list_plots), width = 12, height = 7)

# filter data for minimal coverage threshold:
cell_muts_filtered = cell_muts |>
  filter(sampleID %in% unique(metadata_filtered$sampleID))

# save data
fwrite(cell_muts_filtered, file = paste0("processed_data/", tissue, "/", tissue, "_cell_muts.tsv"))
fwrite(metadata_filtered, paste0("processed_data/", tissue, "/", tissue, "_metadata.tsv"))

# pre-select active signatures present in blood
input_signatures = c("SBS1","SBS5","SBS10a","SBS10b","SBS10c",
                "SBS10d","SBS18","SBS2","SBS13","SBS88","SBS89")

list_results = list()
for (i in unique(metadata$category)) {

  mutrisk_results = mutrisk_pipeline(output_path = "processed_data/colon/",
                                   cell_muts = cell_muts_filtered |> filter(category == i),
                                   metadata = metadata |> filter(category == i),
                                   name = i,
                                   input_signatures = input_signatures,
                                   multiple_refit_methods = TRUE,
                                   sensitivity_correction =  TRUE)
}

# summarize the results from the tissue-specific analysis (see following code from and adapt:)
# summarize the mutation rate
sig_rate_files = list.files(paste0("processed_data/", tissue), pattern = 'sig_rate_per_sample',
                            full.names = TRUE, recursive = TRUE)
rates = lapply(sig_rate_files, fread)
names(rates) = gsub("_sig_rate_per_sample.tsv.gz", "", basename(sig_rate_files))
sig_donor_rates = rbindlist(rates, use.names = TRUE, fill = TRUE) |>
  mutate(category = factor(category, levels = levels(metadata$category))) |>
  group_by(donor, signature) |>
  summarize(across(contains(">"), mean), .groups = "drop") |>
  pivot_longer(contains(">"), names_to = "mut_type", values_to = "mle")
fwrite(sig_donor_rates, file = paste0("processed_data/", tissue, "/", tissue, "_sig_donor_rates.tsv.gz"))

# load the mutation rates with CI
mrates_files = list.files(paste0("processed_data/", tissue), pattern = "_rate_per_sample", full.names = TRUE, recursive = TRUE)
mrates_files = mrates_files[!grepl("sig", mrates_files)] # exclude the signature-specific variants
rates = lapply(mrates_files, fread)
names(rates) = gsub("_rate_per_sample.tsv.gz", "", basename(mrates_files))
expected_rates = rbindlist(rates) |>
  mutate(category = factor(category, levels = levels(metadata$category)))  |>
  inner_join(metadata) |>
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

donors = metadata$donor |> unique()
donors3 =  rbindlist(rates, use.names = TRUE, fill = TRUE) |> pull(donor) |> unique()
donors2  = sig_donor_rates$donor |> unique()
setdiff(donors, donors3)
setdiff(donors, donors2)