# idea: check the expected rate of specific mutations to occur in the ukbiobank population
# take specific variants and check the relative rate
library(readxl)
GE_CRC = readxl::read_excel("raw_data/genomicsengland/41586_2024_7747_MOESM4_ESM.xlsx", sheet = 7, skip = 6)

GE_CRC |>
  filter(Gene == "KRAS") |>
  select(Gene, Mutation, Carriers...8, `Non-carriers...9`) |>
  dplyr::rename(mutated = Carriers...8, non_mutated = `Non-carriers...9`) |>
  mutate(fraction_mutated = mutated / non_mutated)
