#include <Rcpp.h>
#include <vector>

// Helper function to handle recursion efficiently
double inclusion_exclusion_recursive(const std::vector<double>& plist, int start, int end) {
  int n = end - start;
  
  if (n > 2) {
    int mid = start + n / 2;
    
    // Recursive calls on two halves
    double res1 = inclusion_exclusion_recursive(plist, start, mid);
    double res2 = inclusion_exclusion_recursive(plist, mid, end);
    
    // Combine results
    return res1 + res2 - res1 * res2;
  }
  
  if (n == 2) {
    // Base case for two elements
    return plist[start] + plist[start + 1] - plist[start] * plist[start + 1];
  }
  
  if (n == 1) {
    // Base case for one element
    return plist[start];
  }
  
  // Default return value (should not reach here in correct cases)
  return 0.0;
}

// [[Rcpp::export]]
double inclusion_exclusion2(std::vector<double> plist) {
  return inclusion_exclusion_recursive(plist, 0, plist.size());
}
