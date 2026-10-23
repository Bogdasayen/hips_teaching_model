# HIPS Teaching model based on published data in
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Qian Xin October 2026

#' Generate transition matrices for the Markov model (vectorised)
#'
#' A faster implementation of \code{\link{generate_transition_matrices}} that
#' gives an identical array. The original loops over cycles and implants and
#' assigns one slice at a time. Only the death probability varies by cycle, and
#' only the first-revision probability varies by implant, so each transition
#' can be filled across all cycles, implants and samples in one assignment.
#' The life table is looked up for all cycles at once with \code{match()}.
#'
#' @inheritParams generate_transition_matrices
#' @inherit generate_transition_matrices return
#'
#' @importFrom readxl read_excel
#' @export
generate_transition_matrices_vectorised <- function(input_parameters, model_settings) {
  
  n_samples      <- model_settings$n_samples
  n_cycles       <- model_settings$n_cycles
  initial_age    <- model_settings$initial_age
  state_names    <- model_settings$state_names
  implant_names  <- model_settings$implant_names
  data_directory <- model_settings$data_directory
  n_states   <- length(state_names)
  n_implants <- length(implant_names)
  
  transition_matrices <- array(0,
                               dim = c(n_cycles, n_implants, n_samples, n_states, n_states),
                               dimnames = list(NULL, implant_names, NULL, state_names, state_names))
  
  uk_lifetables <- read_excel(
    paste0(data_directory, "/hips_input_data.xlsx"), sheet = "uk_lifetables")
  
  # Death: one probability per cycle, looked up for all cycles at once.
  # Assigning a length-n_cycles vector to a [cycle, implant, sample, state]
  # slice recycles it along the first (cycle) dimension
  ages <- initial_age + (1:n_cycles) - 1
  prob_death <- 1 - exp(-uk_lifetables$Males[match(ages, uk_lifetables$Age)])
  living_states <- state_names[state_names != "dead"]
  transition_matrices[, , , living_states, "dead"] <- prob_death
  
  # First revision: one probability per sample and implant, the same in every
  # cycle. rep(each = n_cycles) lays the samples out along the cycle dimension
  for (implant_name in implant_names) {
    prob_1st_revision <-
      1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_", implant_name)]))
    transition_matrices[, implant_name, , "post_thr", "post_1st_rev"] <-
      rep(prob_1st_revision, each = n_cycles)
  }
  
  # Second revision: one probability per sample, the same for every cycle
  # and implant
  prob_2nd_revision <- 1 - exp(-exp(input_parameters$log_rate_2nd_revision))
  transition_matrices[, , , "post_1st_rev", "post_2nd_rev"] <-
    rep(prob_2nd_revision, each = n_cycles * n_implants)
  
  # Remaining patients stay in their state: diagonal = 1 - sum of the rest,
  # summed over the to-state dimension for all cycles, implants and samples
  for (i_state in 1:n_states) {
    transition_matrices[, , , i_state, i_state] <-
      1 - rowSums(transition_matrices[, , , i_state, -i_state], dims = 3)
  }
  
  # Interim test, as in the original
  row_sums <- rowSums(transition_matrices[sample(n_cycles, 1),
                                          sample(implant_names, 1),
                                          sample(n_samples, 1), , ])
  if (any(abs(row_sums - 1) > 1e-8))
    stop("Rows of the transition matrix must sum to one")
  
  return(transition_matrices)
}
