# exploratory script to model the individual CH rates
# model blood rates:
library(tidyverse)
library(data.table)
# read mutation rates:
mutrates = read.delim("~/Desktop/DNMT3A_HSC_rates.tsv")
B = 1e5
mrate = mutrates$slope[2]
A = 0
s = 0.15
max_time = 70

model_clones = function(B, s, max_time) {

  A_time = numeric(max_time)

  for (i in 1:max_time) {

    cells_mutated <- rpois(1, lambda = B * mrate)
    cells_expanded <- rpois(1, lambda = s * A)
    delta_A <- cells_mutated + cells_expanded

    A = A + delta_A
    B = B - cells_mutated

    A_time[i] = A # save the progression of A over time
  }
  return(A_time)
}







model_multi_clones = function(nreps, B, s, max_time) {

  list_trajectories = list()
  for (i in 1:nreps){
    list_trajectories[[paste0("sim_",i)]] = model_clones_rcpp(B, s, max_time)
  }
  data.frame(list_trajectories)
}

nsims = 100000
multi_clones = model_multi_clones(nsims, B, s, max_time) |>
  mutate(years = 1:max_time) |>
  pivot_longer(-years, names_to = "simulation", values_to = "ncells")

# for tomorrow: Make a new plot where we simulate the rates for a large set of samples
# check if it is possible to make the function faster. For instance, only the year of ncells > 2000 could be saved for instance..

# get the year of appearance:

age_onset = log(4000)/log(1.15)

multi_clones |>
  group_by(simulation) |>
  filter(ncells > 4000)  |>
  ungroup() |>
  dplyr::count(years) |>
  mutate(n = n/nsims) |>
  ggplot(aes(x = years, y = n)) +
  geom_line() +
  geom_abline(slope = mrate * B, intercept = - age_onset * mrate * B, linetype = "dashed") +
  theme_bw() +
  scale_y_continuous(expand = c(0, NA)) +
  scale_x_continuous(limits = c(0, NA))



