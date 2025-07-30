# BoostDM RefCDS
# convert the boostdm matrices into refCDS sites. Enables easy extrapolation of the number of driver mutations
# makes it more easy to determine the actual number of expected drivers for each gene
library(tidyverse)
library(data.table)
library(wintr)
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
boostdm_files = list.files("raw_data/boostdm/", full.names = TRUE)
boostdm_output = data.frame(gene = str_split_i(basename(boostdm_files), "\\.", 1),
           cohort = str_split_i(basename(boostdm_files), "\\.", 3),
           file = boostdm_files)

boostdm_list = split(boostdm_output, boostdm_output$cohort)

boostdm_cohort_list = list()
cohort = names(boostdm_list)[1]
i = 1
for (cohort in names(boostdm_list)) {
  print(cohort)
  cohort_list = list()
  for (i in 1:nrow(boostdm_list[[cohort]])) {
    print(paste(i, ":", boostdm_list[[cohort]]$gene[i]))
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

boostdm_list = rbindlist(boostdm_cohort_list, idcol = 'cohort') |>
  dplyr::rename(gene_name = gene) |>
  filter(gene_name != "CNOT9") |>  # filter mismatch of the gene name
  mutate(gene_name = case_match(gene_name, "H3-3A" ~ "H3F3A", .default = gene_name)) |>
  left_join(gene_strands) |>
  mutate(gstrand = as.character(gstrand))

# first we need to get the original trinucleotide to determine which
gr_boostdm_list = GRanges(seqnames = paste0("chr", boostdm_list$chr),
                          ranges = IRanges(start = boostdm_list$pos, end = boostdm_list$pos),
                        strand = boostdm_list$gstrand)

nuc = c("A", "C", "G", "T")
nuc_match = setNames(nuc, rev(nuc))

boostdm_list = boostdm_list |>
  mutate(context = getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, gr_boostdm_list + 1, as.character = TRUE),
         alt_strand = case_when(gstrand == "-" ~ as.character(nuc_match[alt]), .default = alt))

boostdm_list = boostdm_list |>
  mutate(mut_type = paste0(context, ">", substr(context, 1,1), alt_strand, substr(context, 3,3)))

total_boostdm_muts = boostdm_list |>
  distinct()

# print the position of the protein;
total_boostdm_muts |>
  mutate(distance = pos - lead(pos)) |>
  pull(distance) |> table()

fwrite(total_boostdm_muts, "processed_data/boostdm/BoostDM_context.tsv")

# BoostDM and liftover to hg19
total_boostdm_muts
gr_boostdm_hg38 = GRanges(seqnames = total_boostdm_muts$chr, ranges = IRanges(start = total_boostdm_muts$pos, width = 1))
seqlevelsStyle(gr_boostdm_hg38) = "UCSC"
chain = import.chain("raw_data/resources/hg38ToHg19.over.chain")
gr_boostdm_hg19 = liftOver(gr_boostdm_hg38, chain) |> unlist()
boostdm_hg19 = total_boostdm_muts |>
  mutate(chr = seqnames(gr_boostdm_hg19) |> as.character(),
         pos = start(gr_boostdm_hg19),
         position = parse_number(aachange)) |>
  dplyr::rename(consequence = category)

fwrite(boostdm_hg19, "processed_data/boostdm/boostdm_hg19.txt.gz")

# Perform checks - no output needed
# check if no same-same mutations are happening:
sum(substr(boostdm_list$context, 2,2) == boostdm_list$alt_strand) # sum should be 0 - no mutations should be the same as their reference bases

boostdm_list |>
  group_by(gene_name, cohort, category, mut_type) |>
  count()


gene_values = boostdm_list[,.N, by = c("gene_name", "cohort", "category", "mut_type")] |>
  pivot_wider(names_from = category, values_from = N, values_fill = 0)


names = sapply(RefCDS_WGS, \(x) x$gene_name)
RefCDS_WGS[[which(names == "KRAS")]]$L


dnds_intron_KRAS = RefCDS_WGS[[which(names == "KRAS")]]$L


gene_values |>
  filter(gene_name == "KRAS" & cohort == "COREAD") |>
  select(synonymous, missense, splicing, nonsense) |>
  colSums()



