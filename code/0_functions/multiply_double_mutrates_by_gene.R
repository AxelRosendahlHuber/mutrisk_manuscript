# Compare double mutation rates:
# goal of the functions: Compare the mutation rates for the individual rates
# performed the half-squared approach to calculate the expected number of cells having double cancer driver mutations
# the question is if it is more realistic to allow for a bit of variation in to the different mutational rates

# function: calculate total number of mutations expected:
get_double_gene_rate = function(exp_rates, metadata, site_freqs, ratios, ncells) {
  colnames(site_freqs)[3] = "N" # set count table to unify between data.table and dplyr output
  x = exp_rates |>
    left_join(distinct(metadata), by = c("category", "sampleID")) |>
    inner_join(site_freqs, by = "mut_type", relationship = "many-to-many") |>
    left_join(ratios, by = c("category", "gene_name")) |>
    group_by(category, age, donor, sampleID) |>
    mutate(across(any_of(c("mle", "cilow", "cihigh")), ~ . * N * ratio)) |>
    summarize(across(any_of(c("mle", "cilow", "cihigh")), sum)) |>
    summarize(across(any_of(c("mle", "cilow", "cihigh")), mean), .groups = "drop_last") |>
    mutate(across(any_of(c("mle", "cilow", "cihigh")), ~ (. * (./4)) * ncells)) # calculate the fraction of cells which are expected to be mutated double
}

# Idea to separate the number of cells in the tissue by the number of cells sequenced
# In that case, we can calculate the relative 'fractions' of mutated cells in a given tissue according the observed mutation rates
get_double_gene_rate_fraction = function(exp_rates, metadata, site_freqs, ratios, ncells) {
  colnames(site_freqs)[3] = "N" # set count table to unify between data.table and dplyr output
  y = exp_rates |>
    left_join(distinct(metadata), by = c("category", "sampleID")) |>
    inner_join(site_freqs, by = "mut_type", relationship = "many-to-many") |>
    left_join(ratios, by = c("category", "gene_name")) |>
    group_by(donor) |>
    mutate(across(any_of(c("mle", "cilow", "cihigh")), ~ . * N * ratio)) |>
    group_by(category, age, donor, sampleID) |>
    summarize(across(any_of(c("mle", "cilow", "cihigh")), sum), .groups = "drop_last") |>  # sum all mutation types
    mutate(count = length(unique(sampleID)),  # count number of samples in the sample
           ncells_part = ncells / count) |>            # divide the number of cells by this sample
    mutate(across(any_of(c("mle", "cilow", "cihigh")), ~ (. * (./4)) * ncells_part)) |>  # multiply the number of cells by these rates
    summarize(across(any_of(c("mle", "cilow", "cihigh")), sum), .groups = "drop")
}