test_that("state_costs has dimensions implant x sample x state", {
  expect_equal(dim(state_costs), c(n_implants, n_samples, n_states))
  expect_equal(dimnames(state_costs)[[1]], model_settings$implant_names)
  expect_equal(dimnames(state_costs)[[3]], model_settings$state_names)
})

test_that("state costs are non-negative and zero in the dead state", {
  expect_true(all(state_costs >= 0))
  expect_true(all(state_costs[, , "dead"] == 0))
})

test_that("only the post-THR cost differs between implants", {
  for (state_name in c("post_1st_rev", "post_2nd_rev")) {
    for (i_implant in 2:n_implants) {
      expect_equal(state_costs[i_implant, , state_name], state_costs[1, , state_name])
    }
  }
  expect_false(isTRUE(all.equal(state_costs[1, , "post_thr"], state_costs[2, , "post_thr"])))
})

test_that("state cost equals probability of revision times cost of revision", {
  prob_1st <- 1 - exp(-exp(input_parameters$log_rate_1st_revision_cemented))
  prob_2nd <- 1 - exp(-exp(input_parameters$log_rate_2nd_revision))
  prob_higher <- 1 - exp(-exp(input_parameters$log_rate_higher_revision))
  expect_equal(unname(state_costs["cemented", , "post_thr"]),
               input_parameters$cost_revision * prob_1st)
  expect_equal(unname(state_costs[1, , "post_1st_rev"]),
               input_parameters$cost_revision * prob_2nd)
  expect_equal(unname(state_costs[1, , "post_2nd_rev"]),
               input_parameters$cost_revision * prob_higher)
})
