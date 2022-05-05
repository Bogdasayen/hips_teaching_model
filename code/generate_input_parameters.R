# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom February 2022

# Function to generate random parameters to the Markov model
# These are based on 'parameters' defined in hips_input_data.xlsx

require(readxl)

generate_input_parameters <- function(n_samples) {
  # 15 random parameters in the HIPS teaching model
  input_parameters <- as.data.frame(matrix(nrow = n_samples, ncol = 15))
  colnames(input_parameters) <- c(paste0("log_rate_1st_revision_", implant_names),
                                  "log_rate_2nd_revision", "log_rate_higher_revision",
                                  "cost_revision",
                                  paste0("state_utility_", state_names),
                                  paste0("implant_cost_", implant_names)
  )
  
  # Log rates of first revision for each implant
  log_rates_1st_revision_raw <- read_excel(paste0(data_directory, "/hips_input_data.xlsx"), sheet = "log_rates_1st_revision")
  # Log rates of second and higher revision
  other_log_rates_raw <- read_excel(paste0(data_directory, "/hips_input_data.xlsx"), sheet = "other_log_rates")
  # Implant and revision costs (All fixed)
  costs_raw <- read_excel(paste0(data_directory, "/hips_input_data.xlsx"), sheet = "costs")
  # Utilities for the health states
  state_utilities_raw <- read_excel(paste0(data_directory, "/hips_input_data.xlsx"), sheet = "state_utilities")
  

  
  # Log rates of 1st revision are Normally distributed
  input_parameters$log_rate_1st_revision_cemented <- with(log_rates_1st_revision_raw, rnorm(n_samples, mean = Mean[Implant == "cemented"],
                                                                                            sd = SE[Implant == "cemented"]))
  input_parameters$log_rate_1st_revision_uncemented <- with(log_rates_1st_revision_raw, rnorm(n_samples, mean = Mean[Implant == "uncemented"],
                                                                                              sd = SE[Implant == "uncemented"]))
  input_parameters$log_rate_1st_revision_hybrid <- with(log_rates_1st_revision_raw, rnorm(n_samples, mean = Mean[Implant == "hybrid"],
                                                                                          sd = SE[Implant == "hybrid"]))
  input_parameters$log_rate_1st_revision_reverse_hybrid <- with(log_rates_1st_revision_raw, rnorm(n_samples, mean = Mean[Implant == "reverse_hybrid"],
                                                                                               sd = SE[Implant == "reverse_hybrid"]))

  # Log rates of 2nd and higher revision are Normally distributed
  input_parameters$log_rate_2nd_revision <- with(other_log_rates_raw, rnorm(n_samples, mean =  Value[Parameter == "log_hazard_2nd_revision"],
                                                  sd = Value[Parameter == "log_hazard_2nd_revision_se"]))
  input_parameters$log_rate_higher_revision <- with(other_log_rates_raw, rnorm(n_samples, mean =  Value[Parameter == "log_hazard_higher_revision"],
                                                   sd = Value[Parameter == "log_hazard_higher_revision_se"]))

  # Implant costs are fixed
  input_parameters$implant_cost_cemented <- rep(with(costs_raw, Cost[Resource == "cost_cemented"]), n_samples)
  input_parameters$implant_cost_uncemented <- rep(with(costs_raw, Cost[Resource == "cost_uncemented"]), n_samples)
  input_parameters$implant_cost_hybrid <- rep(with(costs_raw, Cost[Resource == "cost_hybrid"]), n_samples)
  input_parameters$implant_cost_reverse_hybrid <- rep(with(costs_raw, Cost[Resource == "cost_reverse_hybrid"]), n_samples)

  # Cost of revision is fixed
  input_parameters$cost_revision <- rep(with(costs_raw, Cost[Resource == "cost_revision"]), n_samples)
  
  
  # State utilities are Normally distributed
  input_parameters$state_utility_post_thr <- with(state_utilities_raw, rnorm(n_samples, mean = Mean[State == "post_thr"],
                                                                             sd = SE[State == "post_thr"]))
  input_parameters$state_utility_post_1st_rev <- with(state_utilities_raw, rnorm(n_samples, mean = Mean[State == "post_1st_rev"],
                                                                             sd = SE[State == "post_1st_rev"]))
  input_parameters$state_utility_post_2nd_rev <- with(state_utilities_raw, rnorm(n_samples, mean = Mean[State == "post_2nd_rev"],
                                                                                 sd = SE[State == "post_2nd_rev"]))
  # Zero utility in dead state
  input_parameters$state_utility_dead <- 0
  
  return(input_parameters)
}



