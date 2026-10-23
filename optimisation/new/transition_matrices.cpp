// Qian Xin October 2026
//
// Fills the 5-D transition array in C++. R reads the inputs and life table;
// this function does the element-wise filling. Column-major indexing as in
// markov_loop.cpp: element [c, i, s, from, to] of an array with dims
// (C, I, S, N, N) is at c + C*(i + I*(s + S*(from + N*to))).

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fill_transition_matrices_cpp(
    NumericMatrix prob_1st_revision,   // n_samples x n_implants
    NumericVector prob_2nd_revision,   // n_samples
    NumericVector prob_death,          // n_cycles
    int n_cycles, int n_implants, int n_samples, int n_states,
    int post_thr, int post_1st_rev, int post_2nd_rev, int dead) {  // zero-based state indices
  
  NumericVector tm(n_cycles * n_implants * n_samples * n_states * n_states, 0.0);
  auto idx = [&](int c, int i, int s, int from, int to) {
    return c + n_cycles * (i + n_implants * (s + n_samples * (from + n_states * to)));
  };
  
  for (int c = 0; c < n_cycles; ++c) {
    for (int i = 0; i < n_implants; ++i) {
      for (int s = 0; s < n_samples; ++s) {
        tm[idx(c, i, s, post_thr,     post_1st_rev)] = prob_1st_revision(s, i);
        tm[idx(c, i, s, post_1st_rev, post_2nd_rev)] = prob_2nd_revision[s];
        for (int st = 0; st < n_states; ++st) {
          if (st != dead) tm[idx(c, i, s, st, dead)] = prob_death[c];
        }
        // Diagonal: remaining probability, so each row sums to one
        for (int from = 0; from < n_states; ++from) {
          double off = 0.0;
          for (int to = 0; to < n_states; ++to) {
            if (to != from) off += tm[idx(c, i, s, from, to)];
          }
          tm[idx(c, i, s, from, from)] = 1.0 - off;
        }
      }
    }
  }
  tm.attr("dim") = IntegerVector::create(n_cycles, n_implants, n_samples, n_states, n_states);
  return tm;
}
