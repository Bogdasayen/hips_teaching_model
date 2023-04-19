# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom April 2023

#' Function to generate state costs for each implant, sample, and state. This is
#' used internally by the generate_net_benefit() function
#' @param input_parameters Matrix with row for each sample and column 
#' for each parameter, with values samples for the model input parameters
#' @return Array of n_implants x n_samples x n_states, with values equal to costs 
#' for each state for each implant and sample
#' @examples 
#' # First sample the input parameters
#' input_parameters <- generate_input_parameters(n_samples)
#' 
#' # Generate the state costs
#' state_costs <- generate_state_costs(input_parameters)
#' 
#' # Sampled costs for first implant
#' state_costs[1, , ]
#' @export 
generate_state_costs <- function(input_parameters) {

  # Define an array to store state costs for each treatment, sample and state
  state_costs <- array(NA, dim = c(n_implants, n_samples, n_states),
                       dimnames = list(implant_names, NULL, state_names))
  # Cost of 1 year in post total hip replacement depends on annual probability of 1st revision
  # This is implant dependent
  for(implant_name in implant_names) {
    state_costs[implant_name, , "post_thr"] <- input_parameters$cost_revision * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_", implant_name)])))
  }
  # Cost of 1 year in post 1st revision depends on annual probability of 2nd revision
  # It is the same for each implant to repeat each sample several times
  state_costs[, , "post_1st_rev"] <- rep(input_parameters$cost_revision * 
                                           (1 - exp(-exp(input_parameters$log_rate_2nd_revision))), each = n_implants)
  # Cost of 1 year in post 2nd revision depends on annual probability of higher revision
  # It is the same for each implant to repeat each sample several times
  state_costs[, , "post_2nd_rev"] <- rep(input_parameters$cost_revision * 
                                           (1 - exp(-exp(input_parameters$log_rate_higher_revision))), each = n_implants)
  
  # Zero cost of being in the dead state
  state_costs[, , "dead"] <- 0
  
  return(state_costs)
}