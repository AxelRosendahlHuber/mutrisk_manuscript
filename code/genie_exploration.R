# GENIE analyses
# GENIE is probably the only dataset where we can find the different alterations for each of the CRC cancer genes
library(data.table)
library(tidyverse)
nucs = c("A", "C", "G", "T")

genie_data = fread("raw_data/GENIE_17/data_mutations_extended.txt")
genie_data = genie_data |>
  filter(Reference_Allele %in% nucs & Tumor_Seq_Allele2 %in% nucs) |>
  filter(Consequence %in% c("missense_variant", "synonymous_variant", "stop_gained", "start_lost", "stop_lost")) |>
  mutate(aachange = sub("p.", "", HGVSp_Short),
         type = paste0(Reference_Allele, ">", Tumor_Seq_Allele2),
         type = case_when(type == "A>T" ~ "T>A", type == "A>C" ~ "T>G",
                          type == "A>G" ~ "T>C", type == "G>C" ~ "C>G",
                          type == "G>A" ~ "C>T", type == "G>T" ~ "C>A",
                          .default = type)) |>
  select(Hugo_Symbol, Protein_position, aachange, Tumor_Sample_Barcode, type) |>
  dplyr::rename(gene_name = Hugo_Symbol, position = Protein_position, SAMPLE_ID = Tumor_Sample_Barcode)

# merge the data to obtain sample specificity
genie_sample_info = fread("raw_data/GENIE_17/GENIE_release_17.0/data_clinical_sample.txt", skip = 4) |>
  select(SAMPLE_ID, CANCER_TYPE, ONCOTREE_CODE)
genie_cancer_type = inner_join(genie_data, genie_sample_info) |>
  arrange(ONCOTREE_CODE, gene_name) |>
  mutate(genie_present = TRUE)
fwrite(genie_cancer_type, "processed_data/GENIE_17/GENIE_17_processed.txt.gz")


# get for all of the cancer types the mutations:
genie_cancer_type$CANCER_TYPE |> table() |> sort()
genie_crc = genie_cancer_type[CANCER_TYPE == "Colorectal Cancer"]

# number of CRC samples =
unique(genie_crc$SAMPLE_ID) |> length()
table(genie_crc$SAMPLE_ID) |> hist(400)

genie_crc[gene_name == "APC", .N, by = "SAMPLE_ID"] |>
  arrange(N)
