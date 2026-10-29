# Helper for the test suite. testthat sources helper-*.R files before any test
# file, so everything defined here is available to every test.
#

library(readxl)

package_root <- testthat::test_path("..", "..")

# A small model is enough to test behaviour and keeps the suite fast
model_settings <- list(
  n_samples      = 50,
  n_cycles       = 30,
  cycle_length   = 1,
  initial_age    = 65,
  discount_rate  = 0.035,
  state_names    = c("post_thr", "post_1st_rev", "post_2nd_rev", "dead"),
  implant_names  = c("cemented", "uncemented", "hybrid", "reverse_hybrid"),
  data_directory = file.path(package_root, "data")
)
n_samples  <- model_settings$n_samples
n_cycles   <- model_settings$n_cycles
n_states   <- length(model_settings$state_names)
n_implants <- length(model_settings$implant_names)

# Fixed inputs shared by the tests
set.seed(1)
input_parameters    <- generate_input_parameters(model_settings)
transition_matrices <- generate_transition_matrices(input_parameters, model_settings)
state_costs         <- generate_state_costs(input_parameters, model_settings)
state_qalys         <- generate_state_qalys(input_parameters, model_settings)

# Helper so each test only names the argument it changes
run_model <- function(...) {
  generate_net_benefit(transition_matrices, state_costs, state_qalys,
                       input_parameters, model_settings, ...)
}

# Independent reference implementation of the Markov simulation: the original
# nested loops of the monolithic script, written out in full here
# so that it shares no code with the package. Used to check every
# implementation of generate_net_benefit, in R/ and in optimisation/new/.
reference_totals <- function(transition_matrices, state_costs, state_qalys,
                             input_parameters, model_settings) {
  discount_factors <- (1 / (1 + model_settings$discount_rate))^
    (model_settings$cycle_length * (0:(n_cycles - 1)))
  implant_costs <- t(input_parameters[, paste0("implant_cost_",
                                               model_settings$implant_names)])
  ref_costs <- ref_qalys <- matrix(NA, n_implants, n_samples)
  for (i_implant in 1:n_implants) {
    for (i_sample in 1:n_samples) {
      cohort <- matrix(0, n_cycles, n_states)
      cohort[1, 1] <- 1
      for (i_cycle in 2:n_cycles) {
        cohort[i_cycle, ] <- cohort[i_cycle - 1, ] %*%
          transition_matrices[i_cycle - 1, i_implant, i_sample, , ]
      }
      ref_costs[i_implant, i_sample] <- implant_costs[i_implant, i_sample] +
        sum((cohort %*% state_costs[i_implant, i_sample, ]) * discount_factors)
      ref_qalys[i_implant, i_sample] <-
        sum((cohort %*% state_qalys[i_implant, i_sample, ]) * discount_factors)
    }
  }
  list(total_costs = ref_costs, total_qalys = ref_qalys)
}
