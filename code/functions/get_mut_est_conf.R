# Functions to model and plot mutation rates with a confidence interval (as used in the interval model)
# is used in the Blood_Analyses.R script
get_mut_est_conf = function(site_freqs, exp_rates, metadata,
                            ratios, ncells, min_ncells, max_ncells) {

  rates = setNames(c(ncells, min_ncells, max_ncells), c("ncells", "min_ncells", "max_ncells"))

  list_rates = list()
  for (name in names(rates)) {
    print(name)
    n = rates[name]
    list_rates[[name]] = get_gene_rate(exp_rates = exp_rates, metadata = metadata,
                                       site_freqs = site_freqs, ratios = ratios, ncells = n)

  }

  mut_estimates = rbindlist(list_rates, idcol = "cell_estimate")
  estimates_ncells = mut_estimates |>
    filter(gene_name == "DNMT3A") |>
    select(-c(cilow, cihigh)) |>
    pivot_wider(names_from = cell_estimate, values_from = mle)

  model_min = lm(min_ncells ~ age, estimates_ncells)
  model_ncells = lm(ncells ~ age, estimates_ncells)
  model_max = lm(max_ncells ~ age, estimates_ncells)

  estimates = estimates_ncells |>
    mutate(model_min = predict(model_min, estimates_ncells),
           model_ncells = predict(model_ncells, estimates_ncells),
           model_max = predict(model_max, estimates_ncells))
}


