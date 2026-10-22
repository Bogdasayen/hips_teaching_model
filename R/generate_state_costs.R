# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom April 2023
# Modified by Qian Xin October 2026

#' Generate state costs
#' 
#' Computes the expected cost of spending one cycle in each health state, for
#' each implant and probabilistic sample. In each living state the cost is the
#' annual probability of the next revision multiplied by the cost of revision
#' surgery, with the probability converted from the sampled log rate as
#' \code{1 - exp(-exp(log_rate))}. The post-THR cost differs by implant because
#' the first-revision rate is implant-specific; the post-revision costs are the
#' same for every implant. The post second revision state carries the cost of
#' third and higher revisions, which are modelled as a recurring cost rather
#' than a change of state. The dead state has zero cost. Implant purchase costs
#' are not included here; they are added once, at the start of the model, in
#' \code{\link{generate_net_benefit}}.
#'
#' @param input_parameters A data frame with one row per probabilistic sample
#'   and one column per parameter, as returned by
#'   \code{\link{generate_input_parameters}}.
#' @param model_settings A list of model settings. This function uses
#'   \code{n_samples}, \code{state_names}, and \code{implant_names}.
#'
#' @return A numeric array with dimensions n_implants x n_samples x n_states,
#'   holding the cost (in pounds) of one cycle in each state. The implant and
#'   state dimensions are named.
#'
#' @examples
#' model_settings <- list(
#'   n_samples = 10,
#'   state_names    = c("post_thr", "post_1st_rev", "post_2nd_rev", "dead"),
#'   implant_names  = c("cemented", "uncemented", "hybrid", "reverse_hybrid"),
#'   data_directory = "data")
#' input_parameters <- generate_input_parameters(model_settings)
#' state_costs      <- generate_state_costs(input_parameters, model_settings)
#'
#' # Sampled state costs for the cemented implant
#' state_costs["cemented", , ]
#'
#' # Mean cost of a year in each state, by implant
#' apply(state_costs, c(1, 3), mean)
#'
#' @export
generate_state_costs <- function(input_parameters, model_settings) {
  
  # Extract the settings this function needs
  n_samples     <- model_settings$n_samples
  state_names   <- model_settings$state_names
  implant_names <- model_settings$implant_names
  
  # Derive dimensions
  n_states   <- length(state_names)
  n_implants <- length(implant_names)
  
  # One cost for each implant, sample and state
  state_costs <- array(NA, dim = c(n_implants, n_samples, n_states),
                       dimnames = list(implant_names, NULL, state_names))
  
  # Post THR: cost of revision x annual probability of first revision,
  # which depends on the implant
  for (implant_name in implant_names) {
    prob_1st_revision <-
      1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_", implant_name)]))
    state_costs[implant_name, , "post_thr"] <-
      input_parameters$cost_revision * prob_1st_revision
  }
  
  # Post 1st revision: cost of revision x annual probability of second revision,
  # the same for every implant, so each sample's value is repeated n_implants times
  prob_2nd_revision <- 1 - exp(-exp(input_parameters$log_rate_2nd_revision))
  state_costs[, , "post_1st_rev"] <-
    rep(input_parameters$cost_revision * prob_2nd_revision, each = n_implants)
  
  # Post 2nd revision: cost of revision x annual probability of a higher revision
  prob_higher_revision <- 1 - exp(-exp(input_parameters$log_rate_higher_revision))
  state_costs[, , "post_2nd_rev"] <-
    rep(input_parameters$cost_revision * prob_higher_revision, each = n_implants)
  
  # No cost in the dead state
  state_costs[, , "dead"] <- 0
  
  return(state_costs)
}
