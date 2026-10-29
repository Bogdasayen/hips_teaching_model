# Helper for the test suite. testthat sources helper-*.R files before any test
# file, so everything defined here is available to every test.
#

library(readxl)

package_root <- testthat::test_path("..", "..")
for (function_file in list.files(file.path(package_root, "paper", "R"),
                                 pattern = "\\.R$", full.names = TRUE)) {
  source(function_file)
}

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
