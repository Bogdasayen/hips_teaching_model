# Regression tests for the optimised implementations in optimisation/new/.
# They are not part of the package, so they are sourced here. Each must
# reproduce the package functions (and the independent reference) exactly.

optimisation_dir <- file.path(package_root, "optimisation", "new")
skip_if_not(dir.exists(optimisation_dir),
            "optimisation/new/ not available outside the development repository")
source(file.path(optimisation_dir, "generate_transition_matrices_vectorised.R"))
source(file.path(optimisation_dir, "generate_net_benefit_vectorised.R"))

test_that("vectorised transition matrices match the package function", {
  expect_equal(generate_transition_matrices_vectorised(input_parameters, model_settings),
               transition_matrices)
})

test_that("vectorised net benefit matches the package function and the reference", {
  output <- generate_net_benefit_vectorised(transition_matrices, state_costs, state_qalys,
                                            input_parameters, model_settings)
  expect_equal(output, run_model(), tolerance = 1e-10)
  reference <- reference_totals(transition_matrices, state_costs, state_qalys,
                                input_parameters, model_settings)
  expect_equal(unname(output$total_costs), reference$total_costs, tolerance = 1e-10)
  expect_equal(unname(output$total_qalys), reference$total_qalys, tolerance = 1e-10)
})

test_that("vectorised net benefit keeps the input checks", {
  expect_error(generate_net_benefit_vectorised(transition_matrices, state_costs[, , 1:3],
                                               state_qalys, input_parameters, model_settings),
               "unexpected dimensions")
})

test_that("C++ implementations match the package functions", {
  skip_if_not_installed("Rcpp")
  skip_if(Sys.which("g++") == "" && Sys.which("gcc") == "" && .Platform$OS.type != "windows",
          "No C++ compiler available")
  source(file.path(optimisation_dir, "generate_transition_matrices_cpp.R"))
  source(file.path(optimisation_dir, "generate_net_benefit_cpp.R"))
  Rcpp::sourceCpp(file.path(optimisation_dir, "transition_matrices.cpp"))
  Rcpp::sourceCpp(file.path(optimisation_dir, "markov_loop.cpp"))
  
  expect_equal(generate_transition_matrices_cpp(input_parameters, model_settings),
               transition_matrices)
  output <- generate_net_benefit_cpp(transition_matrices, state_costs, state_qalys,
                                     input_parameters, model_settings)
  expect_equal(output, run_model(), tolerance = 1e-10)
})
