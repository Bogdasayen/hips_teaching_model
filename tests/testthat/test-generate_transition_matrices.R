test_that("transition_matrices has dimensions cycle x implant x sample x state x state", {
  expect_equal(dim(transition_matrices),
               c(n_cycles, n_implants, n_samples, n_states, n_states))
  expect_equal(dimnames(transition_matrices)[[2]], model_settings$implant_names)
  expect_equal(dimnames(transition_matrices)[[4]], model_settings$state_names)
  expect_equal(dimnames(transition_matrices)[[5]], model_settings$state_names)
})

test_that("every row of every transition matrix sums to one", {
  # The interim test in the function checks one matrix; here check all of them
  row_sums <- apply(transition_matrices, c(1, 2, 3, 4), sum)
  expect_true(all(abs(row_sums - 1) < 1e-8))
})

test_that("all transition probabilities lie between zero and one", {
  expect_true(all(transition_matrices >= 0 & transition_matrices <= 1))
})

test_that("dead is absorbing and revisions move only forward", {
  expect_true(all(transition_matrices[, , , "dead", "dead"] == 1))
  expect_true(all(transition_matrices[, , , "dead", -n_states] == 0))
  # No transitions back to an earlier revision state, and no skipping a revision
  expect_true(all(transition_matrices[, , , "post_1st_rev", "post_thr"] == 0))
  expect_true(all(transition_matrices[, , , "post_2nd_rev", c("post_thr", "post_1st_rev")] == 0))
  expect_true(all(transition_matrices[, , , "post_thr", "post_2nd_rev"] == 0))
})

test_that("probability of death matches the life table and is the same for every living state", {
  lifetable <- read_excel(file.path(model_settings$data_directory, "hips_input_data.xlsx"),
                          sheet = "uk_lifetables")
  ages <- model_settings$initial_age + (1:n_cycles) - 1
  expected <- 1 - exp(-lifetable$Males[match(ages, lifetable$Age)])
  expect_equal(unname(transition_matrices[, 1, 1, "post_thr", "dead"]), expected)
  expect_true(all(diff(expected) >= 0))
  expect_equal(transition_matrices[, , , "post_1st_rev", "dead"],
               transition_matrices[, , , "post_thr", "dead"])
  expect_equal(transition_matrices[, , , "post_2nd_rev", "dead"],
               transition_matrices[, , , "post_thr", "dead"])
})

test_that("revision probabilities follow the sampled log rates", {
  for (implant_name in model_settings$implant_names) {
    expected <- 1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_", implant_name)]))
    expect_equal(unname(transition_matrices[1, implant_name, , "post_thr", "post_1st_rev"]),
                 expected)
  }
  expected_2nd <- 1 - exp(-exp(input_parameters$log_rate_2nd_revision))
  expect_equal(unname(transition_matrices[1, 1, , "post_1st_rev", "post_2nd_rev"]), expected_2nd)
})

test_that("the function fails clearly if the state names do not include dead", {
  bad_settings <- model_settings
  bad_settings$state_names <- c("post_thr", "post_1st_rev", "post_2nd_rev", "alive")
  expect_error(generate_transition_matrices(input_parameters, bad_settings))
})
