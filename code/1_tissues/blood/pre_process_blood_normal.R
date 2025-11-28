# pre-process normal blood mutations
# min # min tissue run
# parts of the code taken from: https://github.com/emily-mitchell/normal_haematopoiesis
library(data.table)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)

# make output directories:
output_blood = "processed_data/blood/"
output_blood_preprocessed = "processed_data/blood/processed_blood_normal/"
if (!dir.exists(output_blood)) { dir.create(output_blood) }
if (!dir.exists(output_blood_preprocessed)) { dir.create(output_blood_preprocessed) }

# load files"
clone_files = list.files("raw_data/blood/", pattern  = "annotated_mut_set", recursive = TRUE, full.names = TRUE)
clone_files = clone_files[!grepl("textClipping", clone_files)]

for (file in clone_files) {

  load(file)

  name = gsub("filtering_output_", "", strsplit(file, "/")[[1]][5])
  print(name)

  mut_ids = filtered_muts$COMB_mats.tree.build$Genotype_bin == 1

  muts_per_sample = colSums(filtered_muts$COMB_mats.tree.build$Genotype_bin == 1)
  mean(muts_per_sample); sd(muts_per_sample); range(muts_per_sample)

  sample = names(muts_per_sample)[1]
  muts = as.data.table(filtered_muts$COMB_mats.tree.build$mat)
  vafs = filtered_muts$COMB_mats.tree.build$NV /(filtered_muts$COMB_mats.tree.build$NR)

  sample_mut_list = list()
  for (sample in names(muts_per_sample)) {
    print(sample)
    mut_refs = rownames(mut_ids)[mut_ids[,sample]]
    sample_muts = muts[mut_ref %in% mut_refs, 1:6]
    sample_muts$vaf = vafs[mut_refs, sample]
    sample_mut_list[[sample]] = sample_muts
  }

  mutlist_sample = rbindlist(sample_mut_list, idcol = "clone_id") |>
    select(clone_id, Chrom, Pos, Ref, Alt, mut_ref, vaf) |>
    dplyr::rename(sampleID = "clone_id", chr = "Chrom", pos = "Pos", ref = "Ref", alt = "Alt")

  fwrite(mutlist_sample, paste0("processed_data/blood/processed_blood_normal/", name, "_sample_mutations.txt.gz"))
  rm(filtered_muts) ; rm(mut_ids)
  gc()
}
