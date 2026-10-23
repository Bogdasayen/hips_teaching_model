// HIPS Teaching model based on published data in
// Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
// Howard Thom February 2022 (original rcpp_loop_full.cpp)
// Rewritten by Qian Xin October 2026 to take the model's arrays directly,
//   so the R wrapper has the same signature as generate_net_benefit()
//
// Runs the whole Markov simulation for every implant and sample: the cohort
// update over cycles, the per-cycle costs and QALYs, and the discounted
// totals. Only the net benefit calculation is left to R.
//
// The helper functions below do that arithmetic so the loops read like R.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List markov_loop_cpp(NumericVector transition_matrices,  // n_cycles x n_implants x n_samples x n_states x n_states
                     NumericVector state_costs,          // n_implants x n_samples x n_states
                     NumericVector state_qalys,          // n_implants x n_samples x n_states
                     NumericMatrix implant_costs,        // n_implants x n_samples
                     NumericVector discount_factors,     // n_cycles
                     int n_cycles, int n_implants, int n_samples, int n_states,
                     int start_state) {                  // zero-based index of the starting state
  
  NumericMatrix total_costs(n_implants, n_samples);
  NumericMatrix total_qalys(n_implants, n_samples);
  
  // Index into the 5-D transition array
  auto tm_at = [&](int c, int i, int s, int from, int to) {
    return transition_matrices[c + n_cycles * (i + n_implants * (s + n_samples * (from + n_states * to)))];
  };
  // Index into the 3-D state value arrays
  auto sv_at = [&](const NumericVector& x, int i, int s, int st) {
    return x[i + n_implants * (s + n_samples * st)];
  };
  
  std::vector<double> cohort(n_states), cohort_next(n_states);
  
  for (int i = 0; i < n_implants; ++i) {
    for (int s = 0; s < n_samples; ++s) {
      
      // Everyone starts in the starting state
      std::fill(cohort.begin(), cohort.end(), 0.0);
      cohort[start_state] = 1.0;
      
      double cost_sum = 0.0, qaly_sum = 0.0;
      
      for (int c = 0; c < n_cycles; ++c) {
        // Costs and QALYs accrued in this cycle: occupancy x state values, discounted
        double cycle_cost = 0.0, cycle_qaly = 0.0;
        for (int st = 0; st < n_states; ++st) {
          cycle_cost += cohort[st] * sv_at(state_costs, i, s, st);
          cycle_qaly += cohort[st] * sv_at(state_qalys, i, s, st);
        }
        cost_sum += cycle_cost * discount_factors[c];
        qaly_sum += cycle_qaly * discount_factors[c];
        
        // Markov update for the next cycle: cohort_next = cohort %*% P[c]
        if (c < n_cycles - 1) {
          for (int to = 0; to < n_states; ++to) {
            double acc = 0.0;
            for (int from = 0; from < n_states; ++from)
              acc += cohort[from] * tm_at(c, i, s, from, to);
            cohort_next[to] = acc;
          }
          cohort.swap(cohort_next);
        }
      }
      
      total_costs(i, s) = implant_costs(i, s) + cost_sum;
      total_qalys(i, s) = qaly_sum;
    }
  }
  
  return List::create(Named("total_costs") = total_costs,
                      Named("total_qalys") = total_qalys);
}

