# BoostDM RefCDS
# convert the boostdm matrices into refCDS sites. Enables easy extrapolation of the number of driver mutations
# makes it more easy to determine the actual number of expected drivers for each gene
library(tidyverse)
library(data.table)
library(wintr)
library(BSgenome.Hsapiens.UCSC.hg19)
source("code/functions/analysis_variables.R")
# generate 2 outcomes:
# 1. RefCS for all genes (not neccesary driver specific)
# 2. RefCDS for driver mutations. This list is cohort-specific

gene_strands = gr_transcripts |>
  as.data.table() |>
  select(gene_name, strand) |>
  dplyr::rename(gstrand = strand) |>
  mutate(gene_name = case_when(grepl("p16INK4a|p14arf", gene_name) ~ "CDKN2A", .default = gene_name)) |>
  distinct()

# 2. Load the BoostDM files
boostdm_files = list.files("raw_data/boostdm_ch/", full.names = TRUE)
boostdm_output = data.frame(gene = str_split_i(basename(boostdm_files), "\\.", 1),
                            cohort = str_split_i(basename(boostdm_files), "\\.", 2),
                            file = boostdm_files)

boostdm_list = split(boostdm_output, boostdm_output$cohort)

boostdm_cohort_list = list()
cohort = names(boostdm_list)[1]
i = 1
for (cohort in names(boostdm_list)) {
  print(cohort)
  cohort_list = list()
  for (i in 1:nrow(boostdm_list[[cohort]])) {
    print(i)
    bdm = fread(boostdm_list[[cohort]]$file[i]) |>
      select(gene, chr, pos, alt, aachange, starts_with("csqn"), boostDM_class)
    bdm_results = bdm |>  pivot_longer(
      cols = starts_with("csqn"),
      names_to = "category",
      values_to = "value") |>
      filter(value == 1) |>
      mutate(category = sub("csqn_type_", "", category)) |>
      dplyr::select(-value)
    cohort_list[[i]] = bdm_results
  }
  boostdm_cohort_list[[cohort]] = rbindlist(cohort_list)
}

total_boostdm = rbindlist(boostdm_cohort_list, idcol = 'cohort') |>
  dplyr::rename(gene_name = gene) |>
  left_join(gene_strands)

# first we need to get the original trinucleotide to determine which
gr_total_boostdm = GRanges(seqnames = paste0("chr", total_boostdm$chr),
                        ranges = IRanges(start = total_boostdm$pos, end = total_boostdm$pos),
                        strand = total_boostdm$gstrand)

nuc = c("A", "C", "G", "T")
nuc_match = setNames(nuc, rev(nuc))

total_boostdm = total_boostdm |>
  mutate(context = Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, gr_total_boostdm + 1, as.character = TRUE),
         alt_strand = case_when(gstrand == "-" ~ as.character(nuc_match[alt]), .default = alt))

total_boostdm = total_boostdm |>
  mutate(mut_type = paste0(context, ">", substr(context, 1,1), alt_strand, substr(context, 3,3)))

total_boostdm_muts = total_boostdm |>
  distinct()

# BoostDM and liftover to hg19
gr_boostdm_hg38 = GRanges(seqnames = total_boostdm_muts$chr, ranges = IRanges(start = total_boostdm_muts$pos, width = 1))
seqlevelsStyle(gr_boostdm_hg38) = "UCSC"
chain = import.chain("raw_data/resources/hg38ToHg19.over.chain")
gr_boostdm_hg19 = liftOver(gr_boostdm_hg38, chain) |> unlist()
boostdm_hg19 = total_boostdm_muts |>
  mutate(chr = seqnames(gr_boostdm_hg19) |> as.character(),
         pos = start(gr_boostdm_hg19),
         position = parse_number(aachange)) |>
  dplyr::rename(consequence = category)

fwrite(boostdm_hg19, "processed_data/boostdm_ch/boostdm_ch_hg19.txt.gz")

