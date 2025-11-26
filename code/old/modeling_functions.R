# simulate blood age of onset: # model blood rates:


library(Rcpp)
# use a fast Rcpp function for maximum modeling speed
sourceCpp("model_clones.cpp")

# this function is deprecated for a cpp function enabling a >10x speedup
model_clones = function(B, s, max_time, mrate, start_rate) {

  A_time = numeric(max_time)

  A = rpois(1, lambda = B * start_rate) # start with
  for (i in 1:max_time) {

    cells_mutated <- rpois(1, lambda = B * mrate)
    cells_expanded <- rpois(1, lambda = s * A)
    delta_A <- cells_mutated + cells_expanded

    A = A + delta_A
    B = B - cells_mutated

    A_time[i] = A # save the progression of A over time
  }

  if (max(A_time > 2000)) {
    return(which(A_time > 2000))
  }
}


model_multi_clones = function(nreps, B, s, max_time, mrate, start_rate) {

  n_cores <- detectCores() - 1
  results <- mclapply(1:nreps, function(i) model_clones_rcpp(B, s, max_time, mrate, start_rate), mc.cores = n_cores)
  unlist(results)
}

nreps = 1e6


full_modeling = function(nreps, s, max_time, B, mutation) {

  mrate = mutrates |> filter(aachange %in% mutation) |> pull(slope)
  intercept = mutrates |> filter(aachange %in% mutation) |> pull(intecept)
  intercept_individual = intercept * B
  mrate_individual = mrate * B

  output_sims = model_multi_clones(nreps, B, s, max_time, mrate = mrate, start_rate = intercept)

  simulation_df = data.frame(Freq = tabulate(output_sims),
                              age = 1:max_time) |>
    mutate(fraction = Freq / nreps,
           type = "simulation CH") |>
    select(-Freq)

  return(simulation_df)
}