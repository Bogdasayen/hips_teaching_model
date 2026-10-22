# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom April 2023
# Modified by Qian Xin October 2026

#' Generate transition matrices
#' 
#' Builds one transition matrix for every cycle, implant, and probabilistic
#' sample. Revision transitions are derived from the sampled log rates in
#' \code{input_parameters}, converted to annual probabilities as
#' \code{1 - exp(-exp(log_rate))}. Death transitions are taken from the UK
#' life table in \code{hips_input_data.xlsx}, which is read here rather than
#' stored in \code{input_parameters} because it is fixed and age-specific.
#' Each diagonal is set so that every row sums to one. The transition from
#' post second revision to a further revision is not modelled as a change of
#' state; its cost is handled in \code{\link{generate_state_costs}}.
#'
#' @param input_parameters A data frame with one row per probabilistic sample
#'   and one column per parameter, as returned by
#'   \code{\link{generate_input_parameters}}.
#' @param model_settings A list of model settings. This function uses
#'   \code{n_samples}, \code{n_cycles}, \code{initial_age}, \code{state_names},
#'   \code{implant_names}, and \code{data_directory}.
#'
#' @return A numeric array with dimensions n_cycles x n_implants x n_samples x
#' n_states x n_states, holding the probability of moving from each state to each state 
#' in each cycle, for each implant and sample. The implant and state dimensions are named.
#'
#' @examples
#' model_settings <- list(
#'   n_samples = 10, n_cycles = 30, initial_age = 65,
#'   state_names    = c("post_thr", "post_1st_rev", "post_2nd_rev", "dead"),
#'   implant_names  = c("cemented", "uncemented", "hybrid", "reverse_hybrid"),
#'   data_directory = "data")
#' input_parameters    <- generate_input_parameters(model_settings)
#' transition_matrices <- generate_transition_matrices(input_parameters, model_settings)
#'
#' # The transition matrix for the first cycle, first implant and first sample
#' transition_matrices[1, 1, 1, , ]
#'
#' # The same for the first cycle, third implant and second sample
#' transition_matrices[1, 3, 2, , ]
#'
#' # Transitions out of post_thr for all samples, first cycle, cemented implant
#' transition_matrices[1, "cemented", , "post_thr", ]
#'
#' @importFrom readxl read_excel
#' @export
generate_transition_matrices <- function(input_parameters, model_settings) {
  
  # Extract the settings this function needs
  n_samples      <- model_settings$n_samples
  n_cycles       <- model_settings$n_cycles
  initial_age    <- model_settings$initial_age
  state_names    <- model_settings$state_names
  implant_names  <- model_settings$implant_names
  data_directory <- model_settings$data_directory
  
  # Derive dimensions 
  n_states   <- length(state_names)
  n_implants <- length(implant_names)
  
  # One n_states x n_states transition matrix for each cycle, implant and sample
  transition_matrices <- array(0,
                               dim = c(n_cycles, n_implants, n_samples, n_states, n_states),
                               dimnames = list(NULL, implant_names, NULL, state_names, state_names))
  
  # The UK life table is fixed, so it is read here rather than stored in
  # input_parameters, which would repeat it for every sample
  uk_lifetables <- read_excel(
    paste0(data_directory, "/hips_input_data.xlsx"), sheet = "uk_lifetables")
  
  # Living states are all those except dead
  living_states <- state_names[state_names != "dead"]
  
  for (i_cycle in 1:n_cycles) {
    
    # Age of the cohort in this cycle, and the matching annual death rate
    # (looked up by the Age column, not by row position)
    age_this_cycle <- initial_age + i_cycle - 1
    prob_death <- 1 - exp(-uk_lifetables$Males[uk_lifetables$Age == age_this_cycle])
    
    for (implant_name in implant_names) {
      
      # First revision: rate differs by implant
      transition_matrices[i_cycle, implant_name, , "post_thr", "post_1st_rev"] <-
        1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_", implant_name)]))
      
      # Second revision: same rate for all implants
      transition_matrices[i_cycle, implant_name, , "post_1st_rev", "post_2nd_rev"] <-
        1 - exp(-exp(input_parameters$log_rate_2nd_revision))
      
      # Death: same probability from every living state
      transition_matrices[i_cycle, implant_name, , living_states, "dead"] <- prob_death
      
      # Remaining patients stay in their state, so that each row sums to one
      for (i_state in 1:n_states) {
        transition_matrices[i_cycle, implant_name, , i_state, i_state] <-
          1 - rowSums(transition_matrices[i_cycle, implant_name, , i_state, -i_state])
      }
    }
  }
  
  # Interim test: rows of a randomly chosen matrix must sum to one
  row_sums <- rowSums(transition_matrices[sample(n_cycles, 1),
                                          sample(implant_names, 1),
                                          sample(n_samples, 1), , ])
  if (any(abs(row_sums - 1) > 1e-8)) {
    stop("Rows of the transition matrix must sum to one")
  }
  
  return(transition_matrices)
}
