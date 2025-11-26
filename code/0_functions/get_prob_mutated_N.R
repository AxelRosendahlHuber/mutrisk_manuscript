get_prob_mutated_N = function(risk, ncells, N = 1) {
  1 - pbinom(N - 1, size = ncells, prob = risk)
}