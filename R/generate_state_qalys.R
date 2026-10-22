# HIPS Teaching model based on published data in
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Qian Xin October 2026: new function, separating the per-cycle QALYs from
#   generate_net_benefit() so that costs and QALYs are built the same way

#' Generate state QALYs
#'
#' Computes the quality-adjusted life years accrued by spending one cycle in
#' each health state, for each implant and probabilistic sample. This is the
#' sampled state utility multiplied by the cycle length in years. Utilities do
#' not differ by implant, so the same values are repeated for every implant;
#' the implant dimension is kept so that the array has the same shape as the
#' one returned by \code{\link{generate_state_costs}} and the two can be used
#' identically in \code{\link{generate_net_benefit}}. The dead state has zero
#' utility.
#'
#' @param input_parameters A data frame with one row per probabilistic sample
#'   and one column per parameter, as returned by
#'   \code{\link{generate_input_parameters}}.
#' @param model_settings A list of model settings. This function uses
#'   \code{n_samples}, \code{cycle_length}, \code{state_names}, and
#'   \code{implant_names}.
#'
#' @return A numeric array with dimensions n_implants x n_samples x n_states,
#'   holding the QALYs accrued in one cycle in each state. The implant and
#'   state dimensions are named.
#'
#' @examples
#' model_settings <- list(
#'   n_samples = 10, cycle_length = 1,
#'   state_names    = c("post_thr", "post_1st_rev", "post_2nd_rev", "dead"),
#'   implant_names  = c("cemented", "uncemented", "hybrid", "reverse_hybrid"),
#'   data_directory = "data")
#' input_parameters <- generate_input_parameters(model_settings)
#' state_qalys      <- generate_state_qalys(input_parameters, model_settings)
#'
#' # Sampled QALYs per cycle for the cemented implant
#' state_qalys["cemented", , ]
#'
#' # Mean QALYs per cycle in each state (identical across implants)
#' apply(state_qalys, c(1, 3), mean)
#'
#' @export
generate_state_qalys <- function(input_parameters, model_settings) {
  
  # Extract the settings this function needs
  n_samples     <- model_settings$n_samples
  cycle_length  <- model_settings$cycle_length
  state_names   <- model_settings$state_names
  implant_names <- model_settings$implant_names
  
  # Derive dimensions
  n_states   <- length(state_names)
  n_implants <- length(implant_names)
  
  # One value for each implant, sample and state
  state_qalys <- array(NA, dim = c(n_implants, n_samples, n_states),
                       dimnames = list(implant_names, NULL, state_names))
  
  # QALYs per cycle = state utility x cycle length (in years).
  # Utilities are the same for every implant, so each sample's value is
  # repeated n_implants times, as for the post-revision costs
  for (state_name in state_names) {
    state_qalys[, , state_name] <-
      rep(input_parameters[, paste0("state_utility_", state_name)] * cycle_length,
          each = n_implants)
  }
  
  return(state_qalys)
}
