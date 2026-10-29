test_that("generate_net_benefit reproduces the independent nested-loop implementation", {
  reference <- reference_totals(transition_matrices, state_costs, state_qalys,
                                input_parameters, model_settings)
  output <- run_model()
  expect_equal(unname(output$total_costs), reference$total_costs, tolerance = 1e-10)
  expect_equal(unname(output$total_qalys), reference$total_qalys, tolerance = 1e-10)
})

test_that("results match the saved reference output", {
  # Regression test: the totals for these fixed inputs were saved when the
  # model was first verified. A change here means the results have changed.
  reference_file <- test_path("reference_output.rds")
  output <- run_model()
  if (!file.exists(reference_file)) {
    saveRDS(output[c("total_costs", "total_qalys")], reference_file)
    skip("Reference output created; rerun the tests to compare against it")
  }
  reference <- readRDS(reference_file)
  expect_equal(output$total_costs, reference$total_costs, tolerance = 1e-10)
  expect_equal(output$total_qalys, reference$total_qalys, tolerance = 1e-10)
})

test_that("input checks reject malformed arguments", {
  expect_error(generate_net_benefit(transition_matrices, state_costs, state_qalys,
                                    as.matrix(input_parameters), model_settings),
               "must be a data frame")
  expect_error(generate_net_benefit(transition_matrices, state_costs, state_qalys,
                                    input_parameters[1:10, ], model_settings),
               "n_samples rows")
  expect_error(generate_net_benefit(transition_matrices[1:5, , , , ], state_costs,
                                    state_qalys, input_parameters, model_settings),
               "unexpected dimensions")
  expect_error(generate_net_benefit(transition_matrices, state_costs[, , 1:3],
                                    state_qalys, input_parameters, model_settings),
               "unexpected dimensions")
  expect_error(generate_net_benefit(transition_matrices, state_costs,
                                    state_qalys[1:2, , ], input_parameters, model_settings),
               "unexpected dimensions")
  expect_error(run_model(lambda = -1), "non-negative")
  expect_error(run_model(lambda = c(20000, 30000)), "non-negative")
})

test_that("output has the expected structure", {
  output <- run_model()
  expect_type(output, "list")
  expect_named(output, c("total_costs", "total_qalys", "net_benefit",
                         "incremental_net_benefit"))
  for (element in output) {
    expect_equal(dim(element), c(n_implants, n_samples))
    expect_equal(rownames(element), model_settings$implant_names)
  }
  expect_true(all(output$total_costs >= 0))
  expect_true(all(output$total_qalys >= 0))
})

test_that("net benefit follows from costs, QALYs and lambda", {
  lambda <- 30000
  output <- run_model(lambda = lambda)
  expect_equal(output$net_benefit, output$total_qalys * lambda - output$total_costs)
  expect_true(all(output$incremental_net_benefit[1, ] == 0))
  expect_equal(output$incremental_net_benefit,
               sweep(output$net_benefit, 2, output$net_benefit[1, ]))
})

test_that("total costs and QALYs are plausible", {
  output <- run_model()
  # Lifetime QALYs cannot exceed 30 years at full health, and with 3.5%
  # discounting and mortality from age 65 they should be well below that
  expect_true(all(output$total_qalys > 5 & output$total_qalys < 20))
  # Costs at least cover the primary implant
  implant_costs <- t(input_parameters[, paste0("implant_cost_", model_settings$implant_names)])
  expect_true(all(output$total_costs >= implant_costs))
})
