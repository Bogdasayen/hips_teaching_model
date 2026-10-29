test_that("input_parameters has one row per sample and the expected columns", {
  expect_s3_class(input_parameters, "data.frame")
  expect_equal(nrow(input_parameters), n_samples)
  expect_equal(ncol(input_parameters), n_implants + 2 + 1 + n_states + n_implants)
  expect_true(all(c(paste0("log_rate_1st_revision_", model_settings$implant_names),
                    "log_rate_2nd_revision", "log_rate_higher_revision",
                    "cost_revision",
                    paste0("state_utility_", model_settings$state_names),
                    paste0("implant_cost_", model_settings$implant_names)) %in%
                    colnames(input_parameters)))
  expect_false(anyNA(input_parameters))
})

test_that("fixed parameters are constant across samples", {
  expect_equal(length(unique(input_parameters$cost_revision)), 1)
  for (implant_name in model_settings$implant_names) {
    expect_equal(length(unique(input_parameters[, paste0("implant_cost_", implant_name)])), 1)
  }
  expect_true(all(input_parameters$state_utility_dead == 0))
})

test_that("sampled parameters come from the workbook distributions", {
  workbook <- file.path(model_settings$data_directory, "hips_input_data.xlsx")
  log_rates <- read_excel(workbook, sheet = "log_rates_1st_revision")
  utilities <- read_excel(workbook, sheet = "state_utilities")
  # Draw a larger sample so the means can be checked against the workbook
  set.seed(2)
  large_settings <- model_settings
  large_settings$n_samples <- 5000
  large_sample <- generate_input_parameters(large_settings)
  for (implant_name in model_settings$implant_names) {
    expect_equal(mean(large_sample[, paste0("log_rate_1st_revision_", implant_name)]),
                 log_rates$Mean[log_rates$Implant == implant_name], tolerance = 0.01)
  }
  expect_equal(mean(large_sample$state_utility_post_thr),
               utilities$Mean[utilities$State == "post_thr"], tolerance = 0.01)
  # Utilities of the living states should lie in (0, 1]
  living <- grepl("state_utility_post", colnames(large_sample))
  expect_true(all(large_sample[, living] > 0 & large_sample[, living] <= 1))
})

test_that("the same seed reproduces the same samples", {
  set.seed(123)
  first <- generate_input_parameters(model_settings)
  set.seed(123)
  second <- generate_input_parameters(model_settings)
  expect_identical(first, second)
})
