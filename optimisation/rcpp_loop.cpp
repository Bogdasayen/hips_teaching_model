#include <Rcpp.h>
using namespace Rcpp;

// You can source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
// Arguments are the transition matrices in a column of n_cycles stacked matrices 
// and starting cohort vector used to find initial distribution
NumericMatrix rcpp_loop(NumericMatrix cohort_vectors_in, NumericMatrix transition_matrices) {
  // Determine number of cycles
  int n_cycles = cohort_vectors_in.nrow();
  int n_states = cohort_vectors_in.ncol();
  NumericMatrix cohort_vectors_out(n_cycles, n_states);
  // Initialise the output cohort vector
  cohort_vectors_out(0, _) = cohort_vectors_in(0, _);
  NumericVector rm1, cm2;
  NumericMatrix transition_matrices_cycle(n_states, n_states);
  for (int i_cycle = 1; i_cycle < n_cycles; i_cycle++) {
    rm1 = cohort_vectors_out(i_cycle - 1,_);
    transition_matrices_cycle = transition_matrices(Range((i_cycle - 1) *n_states, i_cycle * n_states - 1), Range(0, n_states - 1));
    for (int j_state = 0; j_state < n_states; ++j_state) {
      //cm2 = transition_matrices(_,j_state);
      // Use the correct rows for this cycle
      cm2 = transition_matrices_cycle(_, j_state);
      cohort_vectors_out(i_cycle, j_state) = std::inner_product(rm1.begin(), rm1.end(), cm2.begin(), 0.);              
    }
  }
  return(cohort_vectors_out);
}

