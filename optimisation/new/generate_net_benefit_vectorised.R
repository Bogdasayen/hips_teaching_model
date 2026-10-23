# HIPS Teaching model based on published data in
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Qian Xin October 2026

#' Generate net benefit(vectorised)
#'
#' A faster implementation of \code{\link{generate_net_benefit}} that gives
#' identical results. The Markov update is performed for all implants and
#' samples at once in each cycle, and the per-cycle costs and QALYs are
#' computed with a single array multiplication rather than inside a loop.
#' The remaining loop over cycles is unavoidable, since each cycle depends on
#' the previous one. Input and output tests are the same as in the original.
#'
#' @inheritParams generate_net_benefit
#' @inherit generate_net_benefit return
#'
#' @examples
#' # See generate_net_benefit(); the call is identical
#'
#' @export
generate_net_benefit_vectorised <- function(transition_matrices, state_costs,
                                            state_qalys, input_parameters,
                                            model_settings, lambda = 25000) {
  
  n_samples     <- model_settings$n_samples
  n_cycles      <- model_settings$n_cycles
  cycle_length  <- model_settings$cycle_length
  discount_rate <- model_settings$discount_rate
  state_names   <- model_settings$state_names
  implant_names <- model_settings$implant_names
  n_states   <- length(state_names)
  n_implants <- length(implant_names)
  
  # Tests on inputs
  if (!is.data.frame(input_parameters))
    stop("input_parameters must be a data frame")
  if (nrow(input_parameters) != n_samples)
    stop("input_parameters must have n_samples rows")
  if (!all(dim(transition_matrices) ==
           c(n_cycles, n_implants, n_samples, n_states, n_states)))
    stop("transition_matrices has unexpected dimensions")
  if (!all(dim(state_costs) == c(n_implants, n_samples, n_states)))
    stop("state_costs has unexpected dimensions")
  if (!all(dim(state_qalys) == c(n_implants, n_samples, n_states)))
    stop("state_qalys has unexpected dimensions")
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0)
    stop("lambda must be a single non-negative number")
  
  # Set up (preallocated at full size)
  implant_costs <- t(as.matrix(
    input_parameters[, paste0("implant_cost_", implant_names)]))
  
  cohort_vectors <- array(0, dim = c(n_cycles, n_implants, n_samples, n_states),
                          dimnames = list(NULL, implant_names, NULL, state_names))
  cohort_vectors[1, , , "post_thr"] <- 1
  
  discount_factors <- (1 + discount_rate)^-(cycle_length * (0:(n_cycles - 1)))
  
  # Markov simulation, vectorised across implants and samples
  for (i_cycle in 2:n_cycles) {
    # Extract this cycle's slices once, not once per state
    cohort_prev <- cohort_vectors[i_cycle - 1, , , ]          # implant x sample x state
    tm_cycle    <- transition_matrices[i_cycle - 1, , , , ]   # implant x sample x state x state
    cohort_next <- array(0, dim = dim(cohort_prev))
    for (i_state in 1:n_states) {
      # Arrivals into i_state for all implants and samples at once:
      # sum over source states of (cohort in source) x P(source -> i_state)
      cohort_next[, , i_state] <-
        rowSums(cohort_prev * tm_cycle[, , , i_state], dims = 2)
    }
    cohort_vectors[i_cycle, , , ] <- cohort_next
  }
  
  # Costs and QALYs, vectorised
  # Per-cycle values: multiply occupancy by state values and sum over states.
  # state_costs is implant x sample x state; replicate it along the cycle
  # dimension so it lines up with cohort_vectors (cycle x implant x sample x state)
  expand <- function(x) aperm(array(x, dim = c(dim(x), n_cycles)), c(4, 1, 2, 3))
  cycle_costs <- rowSums(cohort_vectors * expand(state_costs), dims = 3)   # cycle x implant x sample
  cycle_qalys <- rowSums(cohort_vectors * expand(state_qalys), dims = 3)
  
  # Discount over cycles and sum: multiply each cycle's slice by its factor
  total_costs <- implant_costs + apply(cycle_costs * discount_factors, c(2, 3), sum)
  total_qalys <- apply(cycle_qalys * discount_factors, c(2, 3), sum)
  dimnames(total_costs) <- dimnames(total_qalys) <- list(implant_names, NULL)
  
  # Tests on outputs
  if (!all(dim(total_costs) == c(n_implants, n_samples)))
    stop("total_costs has unexpected dimensions")
  if (!all(dim(total_qalys) == c(n_implants, n_samples)))
    stop("total_qalys has unexpected dimensions")
  if (any(total_costs < 0))
    stop("total_costs contains negative values")
  if (any(total_qalys < 0) || any(total_qalys > n_cycles * cycle_length))
    stop("total_qalys outside the plausible range")
  
  # Net benefit
  net_benefit <- total_qalys * lambda - total_costs
  incremental_net_benefit <- sweep(net_benefit, 2, net_benefit[1, ])
  
  return(list(total_costs             = total_costs,
              total_qalys             = total_qalys,
              net_benefit             = net_benefit,
              incremental_net_benefit = incremental_net_benefit))
}
