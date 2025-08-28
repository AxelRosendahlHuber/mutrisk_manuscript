#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector model_clones_rcpp(double B, double s, int max_time, double mrate, double start_rate) {
  NumericVector A_time(max_time);

  // Initial A
  double A = R::rpois(B * start_rate);

  for (int i = 0; i < max_time; i++) {
    int cells_mutated = R::rpois(B * mrate);
    int cells_expanded = R::rpois(s * A);
    double delta_A = cells_mutated + cells_expanded;

    A += delta_A;
    B -= cells_mutated;

    A_time[i] = A;
  }

  // Check if any A exceeds 2000
  if (max(A_time) > 2000) {
    // Return 1-based indices where A_time > 2000
    std::vector<int> indices;
    for (int i = 0; i < max_time; i++) {
      if (A_time[i] > 2000) {
        indices.push_back(i + 1); // R is 1-indexed
      }
    }
    return wrap(indices);
  }

  return IntegerVector::create(); // return empty vector if none exceed 2000
}
