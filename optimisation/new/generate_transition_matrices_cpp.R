# HIPS Teaching model based on published data in
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Qian Xin October 2026

#' Generate transition matrices for the Markov model (C++)
#'
#' A compiled implementation of \code{\link{generate_transition_matrices}}
#' that gives an identical array. R converts the sampled log rates and the
#' life table to probabilities; the element-wise filling of the five-
#' dimensional array is done in C++ (\code{transition_matrices.cpp}).
#'
#' @inheritParams generate_transition_matrices
#' @inherit generate_transition_matrices return
#'
#' @importFrom readxl read_excel
#' @export
generate_transition_matrices_cpp <- function(input_parameters, model_settings) {
  
  n_samples      <- model_settings$n_samples
  n_cycles       <- model_settings$n_cycles
  initial_age    <- model_settings$initial_age
  state_names    <- model_settings$state_names
  implant_names  <- model_settings$implant_names
  data_directory <- model_settings$data_directory
  n_states   <- length(state_names)
  n_implants <- length(implant_names)
  
  uk_lifetables <- read_excel(
    paste0(data_directory, "/hips_input_data.xlsx"), sheet = "uk_lifetables")
  ages <- initial_age + (1:n_cycles) - 1
  prob_death <- 1 - exp(-uk_lifetables$Males[match(ages, uk_lifetables$Age)])
  
  prob_1st_revision <- as.matrix(
    1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_", implant_names)])))
  prob_2nd_revision <- 1 - exp(-exp(input_parameters$log_rate_2nd_revision))
  
  state_index <- function(name) which(state_names == name) - 1   # zero-based for C++
  transition_matrices <- fill_transition_matrices_cpp(
    prob_1st_revision, prob_2nd_revision, prob_death,
    n_cycles, n_implants, n_samples, n_states,
    state_index("post_thr"), state_index("post_1st_rev"),
    state_index("post_2nd_rev"), state_index("dead"))
  dimnames(transition_matrices) <-
    list(NULL, implant_names, NULL, state_names, state_names)
  
  # Interim test, as in the original
  row_sums <- rowSums(transition_matrices[sample(n_cycles, 1),
                                          sample(implant_names, 1),
                                          sample(n_samples, 1), , ])
  if (any(abs(row_sums - 1) > 1e-8))
    stop("Rows of the transition matrix must sum to one")
  
  return(transition_matrices)
}
