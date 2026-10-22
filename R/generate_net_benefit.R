# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom April 2023
# Modified by Qian Xin October 2026

#' Generate net benefit
#' 
#' Runs the Markov cohort simulation and returns discounted total costs, total
#' QALYs, net monetary benefit, and incremental net benefit for each implant
#' and probabilistic sample. The cohort starts in the post-THR state and is
#' traced through the states by multiplying the cohort vector by the
#' transition matrix in each cycle. Costs and QALYs accrued in each cycle are
#' the state occupancy multiplied by the per-cycle state costs and QALYs,
#' discounted at the rate in \code{model_settings}, and summed over the time
#' horizon. The implant purchase cost is added once, at the start.
#'
#' The function takes the model's components as arguments rather than
#' building them, so its logic can be tested on known inputs without rebuilding
#' the rest of the model. It checks each argument before the simulation begins.
#'
#' @param transition_matrices Numeric array, n_cycles x n_implants x n_samples
#'   x n_states x n_states, as returned by
#'   \code{\link{generate_transition_matrices}}.
#' @param state_costs Numeric array, n_implants x n_samples x n_states, as
#'   returned by \code{\link{generate_state_costs}}.
#' @param state_qalys Numeric array, n_implants x n_samples x n_states, as
#'   returned by \code{\link{generate_state_qalys}}.
#' @param input_parameters A data frame with one row per probabilistic sample
#'   and one column per parameter, as returned by
#'   \code{\link{generate_input_parameters}}. Used here for the implant costs.
#' @param model_settings A list of model settings. This function uses
#'   \code{n_samples}, \code{n_cycles}, \code{cycle_length},
#'   \code{discount_rate}, \code{state_names}, and \code{implant_names}.
#' @param lambda Willingness-to-pay threshold in pounds per QALY, a single
#'   non-negative number. Default 25000.
#'
#' @return A list of four numeric matrices, each n_implants x n_samples:
#'   \code{total_costs}, \code{total_qalys}, \code{net_benefit}, and
#'   \code{incremental_net_benefit} (relative to the first implant).
#'
#' @examples
#' model_settings <- list(
#'   n_samples = 10, n_cycles = 30, cycle_length = 1, initial_age = 65,
#'   discount_rate = 0.035,
#'   state_names    = c("post_thr", "post_1st_rev", "post_2nd_rev", "dead"),
#'   implant_names  = c("cemented", "uncemented", "hybrid", "reverse_hybrid"),
#'   data_directory = "data")
#' input_parameters    <- generate_input_parameters(model_settings)
#' transition_matrices <- generate_transition_matrices(input_parameters, model_settings)
#' state_costs         <- generate_state_costs(input_parameters, model_settings)
#' state_qalys         <- generate_state_qalys(input_parameters, model_settings)
#' model_output <- generate_net_benefit(transition_matrices, state_costs,
#'                                      state_qalys, input_parameters,
#'                                      model_settings, lambda = 25000)
#'
#' # Average costs, QALYs and incremental net benefit by implant
#' rowMeans(model_output$total_costs)
#' rowMeans(model_output$total_qalys)
#' rowMeans(model_output$incremental_net_benefit)
#'
#' @export
generate_net_benefit <- function(transition_matrices, state_costs,
                                 state_qalys, input_parameters,
                                 model_settings, lambda = 25000) {
  
  # Extract the settings this function needs
  n_samples     <- model_settings$n_samples
  n_cycles      <- model_settings$n_cycles
  cycle_length  <- model_settings$cycle_length
  discount_rate <- model_settings$discount_rate
  state_names   <- model_settings$state_names
  implant_names <- model_settings$implant_names
  
  # Derive dimensions
  n_states   <- length(state_names)
  n_implants <- length(implant_names)
  
  # Tests on inputs
  # input_parameters must be a data frame with one row per sample
  if (!is.data.frame(input_parameters))
    stop("input_parameters must be a data frame")
  if (nrow(input_parameters) != n_samples)
    stop("input_parameters must have n_samples rows")
  
  # transition_matrices must have the expected five dimensions
  if (!all(dim(transition_matrices) ==
           c(n_cycles, n_implants, n_samples, n_states, n_states)))
    stop("transition_matrices has unexpected dimensions")
  
  # state_costs and state_qalys must match implant x sample x state
  if (!all(dim(state_costs) == c(n_implants, n_samples, n_states)))
    stop("state_costs has unexpected dimensions")
  if (!all(dim(state_qalys) == c(n_implants, n_samples, n_states)))
    stop("state_qalys has unexpected dimensions")
  
  # lambda must be a single non-negative number
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0)
    stop("lambda must be a single non-negative number")
  
  # Set up 
  # Implant purchase costs, one per implant and sample
  # (transposed to the n_implants x n_samples convention used throughout)
  implant_costs <- t(as.matrix(
    input_parameters[, paste0("implant_cost_", implant_names)]))
  
  # Cohort vector for each cycle, implant and sample: everyone starts in post_thr
  cohort_vectors <- array(0, dim = c(n_cycles, n_implants, n_samples, n_states),
                          dimnames = list(NULL, implant_names, NULL, state_names))
  cohort_vectors[1, , , "post_thr"] <- 1
  
  # Costs and QALYs accrued in each cycle, and their discounted totals
  cycle_costs <- array(dim = c(n_cycles, n_implants, n_samples),
                       dimnames = list(NULL, implant_names, NULL))
  cycle_qalys <- array(dim = c(n_cycles, n_implants, n_samples),
                       dimnames = list(NULL, implant_names, NULL))
  total_costs <- array(dim = c(n_implants, n_samples),
                       dimnames = list(implant_names, NULL))
  total_qalys <- array(dim = c(n_implants, n_samples),
                       dimnames = list(implant_names, NULL))
  
  # Discount factor for each cycle: 1 in the first cycle, then
  # 1 / (1 + rate)^years elapsed, where years = cycle_length x cycles elapsed
  discount_factors <- (1 + discount_rate)^-(cycle_length * (0:(n_cycles - 1)))
  
  # Markov simulation 
  for (i_implant in 1:n_implants) {
    for (i_sample in 1:n_samples) {
      
      # Markov update: multiply the previous cohort vector by the transition
      # matrix. Cycle 1 is already defined, so update cycles 2 to n_cycles
      for (i_cycle in 2:n_cycles) {
        cohort_vectors[i_cycle, i_implant, i_sample, ] <-
          cohort_vectors[i_cycle - 1, i_implant, i_sample, ] %*%
          transition_matrices[i_cycle - 1, i_implant, i_sample, , ]
      }
      
      # Costs and QALYs in each cycle: state occupancy x per-cycle values
      cycle_costs[, i_implant, i_sample] <-
        cohort_vectors[, i_implant, i_sample, ] %*% state_costs[i_implant, i_sample, ]
      cycle_qalys[, i_implant, i_sample] <-
        cohort_vectors[, i_implant, i_sample, ] %*% state_qalys[i_implant, i_sample, ]
      
      # Discount and sum; implant purchase cost is added once, undiscounted
      total_costs[i_implant, i_sample] <- implant_costs[i_implant, i_sample] +
        sum(cycle_costs[, i_implant, i_sample] * discount_factors)
      total_qalys[i_implant, i_sample] <-
        sum(cycle_qalys[, i_implant, i_sample] * discount_factors)
    }
  }
  
  # Tests on outputs 
  # Totals have one row per implant and one column per sample
  if (!all(dim(total_costs) == c(n_implants, n_samples)))
    stop("total_costs has unexpected dimensions")
  if (!all(dim(total_qalys) == c(n_implants, n_samples)))
    stop("total_qalys has unexpected dimensions")
  
  # Discounted lifetime costs and QALYs should be non-negative, and QALYs
  # cannot exceed the undiscounted time horizon in years
  if (any(total_costs < 0))
    stop("total_costs contains negative values")
  if (any(total_qalys < 0) || any(total_qalys > n_cycles * cycle_length))
    stop("total_qalys outside the plausible range")
  
  # Net benefit
  net_benefit <- total_qalys * lambda - total_costs
  
  # Incremental net benefit relative to the first implant, subtracting the
  # reference row from every row, sample by sample
  incremental_net_benefit <- sweep(net_benefit, 2, net_benefit[1, ])
  
  return(list(total_costs             = total_costs,
              total_qalys             = total_qalys,
              net_benefit             = net_benefit,
              incremental_net_benefit = incremental_net_benefit))
}
