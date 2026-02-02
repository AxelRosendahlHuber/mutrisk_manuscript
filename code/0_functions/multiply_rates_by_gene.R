# function: calculate total number of mutations expected:
get_gene_rate = function(exp_rates, metadata, site_freqs, ratios, ncells) {
  colnames(site_freqs)[3] = "N" # set count table to unify between data.table and dplyr output
  exp_rates |>
    left_join(distinct(metadata), by = c("category", "sampleID")) |>
    inner_join(site_freqs, by = "mut_type", relationship = "many-to-many") |>
    left_join(ratios, by = c("category", "gene_name")) |>
    group_by(category, age, gene_name, donor, sampleID) |>
    mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio * ncells)) |>
    summarize(across(c(mle, cilow, cihigh), sum), .groups = "drop_last") |>
    summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop")
}

# Idea to better represent the diversity within samples
# in this case, we multiply as if the mutation rates ar equal to the observed number of cells in a patient
get_gene_rate_clone = function(exp_rates, metadata, site_freqs, ratios, ncells) {
  colnames(site_freqs)[3] = "N" # set count table to unify between data.table and dplyr output
  exp_rates |>
    left_join(distinct(metadata), by = c("category", "sampleID")) |>
    inner_join(site_freqs, by = "mut_type", relationship = "many-to-many") |>
    left_join(ratios, by = c("category", "gene_name")) |>
    mutate(across(c(mle, cilow, cihigh), ~ . * N * ratio * ncells)) |>  # multiply the number of cells by these rates
    group_by(category, age, donor, sampleID) |>
    summarize(across(c(mle, cilow, cihigh), sum))  # sum all mutation types
}

# Idea to separate the number of cells in the tissue by the number of cells sequenced
# In that case, we can calculate the relative 'fractions' of mutated cells in a given tissue according the observed mutation rates
get_gene_rate_fraction = function(exp_rates, metadata, site_freqs, ratios, ncells) {
  colnames(site_freqs)[3] = "N" # set count table to unify between data.table and dplyr output
  exp_rates |>
    left_join(distinct(metadata), by = c("category", "sampleID")) |>
    inner_join(site_freqs, by = "mut_type", relationship = "many-to-many") |>
    left_join(ratios, by = c("category", "gene_name")) |>
    group_by(donor) |>
    mutate(count = length(unique(sampleID)),  # count number of samples in the sample
           ncells = ncells / count,           # divide the number of cells by this sample
           across(c(mle, cilow, cihigh), ~ . * N * ratio * ncells)) |>  # multiply the number of cells by these rates
    group_by(category, age, donor, sampleID) |>
    summarize(across(c(mle, cilow, cihigh), sum), .groups = "drop_last") |>  # sum all mutation types
    summarize(across(c(mle, cilow, cihigh), sum), .groups = "drop")
}


# same function, but now also group by signature
get_gene_rate_sig = function(exp_rates, metadata, site_freqs, ratios) {

  colnames(site_freqs)[3] = "N" # set count table to unify between data.table and dplyr output
  # correct the trinucleotide rates for
  exp_rates_trinuc = exp_rates |>
    left_join(distinct(metadata))

  exp_rates_trinuc |>
    inner_join(site_freqs, by = "mut_type") |>
    left_join(ratios, by = c("category", "gene_name")) |>
    group_by(category, age, gene_name, donor, signature, sampleID) |>
    mutate(mle = mle * N * ratio) |>
    summarize(mle = sum(mle), .groups = "drop_last") |>
    summarize(mle = mean(mle), .groups = "drop")
}

