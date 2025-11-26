# functions to calculate the fraction of a given gene expected to be mutated

# calculate the probability for a position to be mutated given there are n cells
get_prob_mutated = function(risk, ncells) {
  1  - ((1 - risk)^ncells)
}

# adjust the mutation rates based on the relative mutagenicity
get_adjusted_rates = function(expected_rates_sample, gene_counts, ratios_cat) {

  trinuc_vec = setNames(expected_rates_sample$mle, expected_rates_sample$mut_type)
  ratio_vec = setNames(ratios_cat$ratio, ratios_cat$gene_name)

  setDT(gene_counts)
  gene_counts = gene_counts[count > 0,]

  ratios_gene = ratio_vec[gene_counts$gene_name]
  trinucs_gene = trinuc_vec[gene_counts$mut_type]
  adjusted_rates = ratios_gene * trinucs_gene

  return(list(adjusted_rates = adjusted_rates, count = gene_counts$count))
}


# calculate the mutation probablity across a range of values
calc_prob_mut_range = function(adjusted_rates, range) {

  prob_results = data.frame(ncells = range, prob = NA)

  for (i in 1:length(range)) {
    ncells = round(range[i])

    rates = get_prob_mutated(adjusted_rates$adjusted_rates, ncells)
    prob_results[i,2] = weighted.mean(rates, adjusted_rates$count)
  }

  return(prob_results)
}

# get the probability across a range. Input the expected rates for a single individual
get_prob_mutated_range = function(expected_rates_sample, gene_counts, ratios_cat, range, step = 0.1) {

  if (missing(range)) {
    range = 10^(seq(4,10, step)) # if a range is not provided, create an own set of breaks, separated on a logarithmic scale
  }

  adjusted_rates = get_adjusted_rates(expected_rates_sample = expected_rates_sample,
                                      gene_counts = gene_counts, ratios_cat = ratios_cat)

  calc_prob_mut_range(adjusted_rates = adjusted_rates, range = range)
}

### EXPERIMENTAL - plotting figures
# consider not adding this in the figure, and make it possible for every indivudal to plot separately

# function to plot the results (experimental and needs to be further optimized)
plot_prob_range = function(prob_results) {
  prob_results = prob_results |>
    as.data.frame()

  prob_results |>
    ggplot() +
    geom_line(aes(x = ncells, y = prob,), alpha = 0.5, linewidth = 0.5) +
    scale_x_log10() +
    cowplot::theme_cowplot()
}
