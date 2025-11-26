# Supplementary Note 3 Figure A

# model a population of cells which can grow and divide
# function to model expansions (grow a cell over a specific time interval)
grow_pop = function(pop, division_rate, dt, s = 1 ) {
  pop_division_rate = pop * division_rate * dt
  n_div <- rpois(1, pop_division_rate)
  d  = 1 + 1 + s # denominator
  probs = c(s/d, 1/d, 1/d)
  if (n_div > 0) {
    # For each division event, choose an outcome:
    # "expansion" gives +1 cell, "no_expansion" gives 0, "death" gives -1.
    outcomes <- sample(c(1, 0, -1),
                       size = n_div,
                       replace = TRUE,
                       prob = probs)

    # Update population: no_expansion events do not change the count.
    pop <- pop + sum(outcomes)
  }
  pop[pop < 0] = 0 # reset population to 0 every time there is a 'negative event'
  pop
}


set.seed(12356)

ncells_sim = 1e5
cells = rep(1, ncells_sim)

cells_100_years = vector("numeric", ncells_sim)
sum_cells = vector("numeric",100)
ncells = vector("numeric", 100)
for (year in 1:1001) {
  print(year)
  ncells[year] = length(cells)
  sum_cells[year] = sum(cells)
  cells = sapply(cells, grow_pop, division_rate = 1.3, dt = 0.1)
  cells = cells[cells > 0]
}

hist(cells)

plot(sum_cells)
plot(ncells)
hsc_simulation_df = tibble(age = seq(0, 100, 0.1),
                               `cells in population` = sum_cells,
                               `surviving HSC lineages` = ncells)

hsc_simulation_df |>
  pivot_longer(-age) |>
  ggplot(aes(x = age)) +
  geom_line(aes(y = value, color = name)) +
  annotate("text", label = paste0("HSCs age 100:\n", ncells[1001]), x = 100, y = ncells[1001], hjust = 1, vjust = -.2) +
  theme_cowplot() +
  theme(legend.position = "inside", legend.position.inside = c(0.5, 0.5)) +
  labs(x = "age (years)", y = "number of cells", subtitle = "simulation 100,000 HSCs\nDivision rate of 1.3/year",
       color = "number of:")
ggsave("manuscript/Supplementary_notes/Supplementary_Note_III/SNIII_Figure_1.png", width = 5, height = 4, bg = "white")