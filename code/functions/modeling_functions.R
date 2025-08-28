# simulate blood age of onset: # model blood rates:
# this function is deprecated for a cpp function enabling a >10x speedup
model_clones = function(B, s, max_time, mrate, start_rate) {

  VAF_threshold = B * 0.04
  A = rpois(1, lambda = B * start_rate) # start with intercept rate
  B_model = B
  for (i in 1:max_time) {

    cells_mutated <- rpois(1, lambda = B_model * mrate)
    cells_expanded <- rpois(1, lambda = s * A)
    delta_A <- cells_mutated + cells_expanded

    A = A + delta_A
    B_model = B_model - cells_mutated
    if (B_model <= 0) {B_model = 0}

    if (A > VAF_threshold)   {break}
  }


  if (A > VAF_threshold)   { return(i)
  } else {return(NULL)}
}

model_multi_clones = function(nreps, B, s, max_time, mrate, start_rate) {

  n_cores <- detectCores() - 1
  results <- mclapply(1:nreps, function(i) model_clones_rcpp(B, s, max_time, mrate, start_rate), mc.cores = n_cores)
  unlist(results)
}


full_modeling = function(nreps, s, max_time, B, mrate, intercept) {

  output_sims = model_multi_clones(nreps, B, s, max_time, mrate = mrate, start_rate = intercept)

  total_sum = numeric(max_time)
  if (!is.null(output_sims)) {
    total_counts = tabulate(output_sims)
    total_sum[1:length(total_counts)] = total_counts
  }

  simulation_df = data.frame(Freq = cumsum(total_sum),
                              age = 1:max_time) |>
    mutate(fraction = Freq / nreps,
           type = "simulation CH") |>
    select(-Freq)

  return(simulation_df)
}



library(Rcpp)

cppFunction('
int model_clones_rcpp(int B, double s, int max_time, double mrate, double start_rate) {
    double VAF_threshold = B * 0.04;
    int A = R::rpois(B * start_rate); // start with intercept rate
    int B_model = B;

    for (int i = 1; i <= max_time; i++) {
        int cells_mutated = R::rpois(B_model * mrate);
        int cells_expanded = R::rpois(s * A);
        int delta_A = cells_mutated + cells_expanded;

        A += delta_A;
        B_model -= cells_mutated;
        if (B_model <= 0) B_model = 0;

        if (A > VAF_threshold) return i;
    }

    return NA_INTEGER; // same as returning NULL in R
}
')

