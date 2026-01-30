merge_mutrisk_drivers = function(boostdm, ratios, expected_rates,  gene_of_interest, tissue_select = "colon", tissue_name,
                                 category_select = "normal", cell_probabilities = FALSE,
                                 individual = FALSE, filter_age = TRUE) {

  older_individuals = metadata |> filter(tissue == tissue_select,
                                         category %in% category_select)

  if (filter_age == TRUE) {older_individuals = older_individuals |> filter(age > 30)}
  ratio_gene_tissue = ratios |> filter(gene_name == gene_of_interest &
                                         category %in% category_select,
                                       tissue == tissue_select) |> pull(ratio)
  expected_rates_select = expected_rates[category %in% category_select &
                                           tissue == tissue_select, ] |>
    left_join(older_individuals, by = c("tissue", "sampleID", "category")) |>
    filter(donor %in% older_individuals$donor)

  ncells_select = tissue_ncells_ci[tissue_ncells_ci$tissue == tissue_select, "mid_estimate"]

  if (cell_probabilities == TRUE) {
    ncells_select = 1
  }

  # group by donor individually
  mutated_rates = expected_rates_select |>
    group_by(donor, mut_type, tissue) |>
    summarize(across(c(mle, cilow, cihigh), mean), .groups = "drop") |>
    mutate(across(c(mle, cilow, cihigh), ~ . * ratio_gene_tissue * ncells_select))

  # modify for specific individual, or for the entiriety
  if (individual == FALSE) {
    print("taking mean of the mutation rates")

    mutated_rates_select = mutated_rates |>
      group_by(mut_type, tissue) |>
      summarize(mle = mean(mle))

    individuals = older_individuals |>
      select(donor, age) |> distinct()

    label = paste(tissue_name, "- average age:", format(mean(individuals$age), digits = 3))


  } else if (individual %in% unique(mutated_rates$donor)) {
    mutated_rates_select = mutated_rates |>
      filter(donor == individual) |>
      select(mut_type, tissue, mle)

    label = paste(category_select, "donor", individual, " age:", older_individuals[donor == individual] |> pull(age))

  }  else if (individual == "all") {
    mutated_rates_select = mutated_rates
    label = "no_label"
  }   else {print("parameter 'individual' must either be FALSE, or donor-id")}

  boostdM_goi = boostdm[[tissue_select]][gene_name == gene_of_interest, c("mut_type", "position", "driver")]
  expected_gene_muts = boostdM_goi |>
    full_join(mutated_rates_select, relationship = "many-to-many", by = "mut_type")

  return(list(expected_gene_muts = expected_gene_muts, label = label))
}


