test_that("state_qalys has dimensions implant x sample x state", {
  expect_equal(dim(state_qalys), c(n_implants, n_samples, n_states))
  expect_equal(dimnames(state_qalys)[[1]], model_settings$implant_names)
  expect_equal(dimnames(state_qalys)[[3]], model_settings$state_names)
})

test_that("state QALYs equal utility times cycle length and are the same for every implant", {
  utilities <- as.matrix(input_parameters[, paste0("state_utility_", model_settings$state_names)])
  for (i_implant in 1:n_implants) {
    expect_equal(unname(state_qalys[i_implant, , ]),
                 unname(utilities) * model_settings$cycle_length)
  }
  expect_true(all(state_qalys[, , "dead"] == 0))
})

test_that("cycle length scales the QALYs", {
  half_year <- model_settings
  half_year$cycle_length <- 0.5
  expect_equal(generate_state_qalys(input_parameters, half_year), state_qalys / 2)
})
