# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom February 2022
# Modified by Qian Xin October 2026

# Function to use Markov modelling to simulate total costs, QALYs and net benefit

#' A compiled implementation of \code{\link{generate_net_benefit}} that gives
#' identical results. The Markov update, per-cycle costs and QALYs, and
#' discounting are all performed in C++ via Rcpp; R does the argument checks,
#' passes the arrays through, and computes net benefit from the totals.
#'
#' @inheritParams generate_net_benefit
#' @inherit generate_net_benefit return
#'
#' @export
generate_net_benefit_cpp <- function(transition_matrices, state_costs,
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
  
  implant_costs <- t(as.matrix(
    input_parameters[, paste0("implant_cost_", implant_names)]))
  discount_factors <- (1 + discount_rate)^-(cycle_length * (0:(n_cycles - 1)))
  
  # Markov simulation in C++ 
  # The arrays are passed as they are; R's column-major layout is handled
  # inside markov_loop.cpp
  totals <- markov_loop_cpp(
    transition_matrices = transition_matrices,
    state_costs         = state_costs,
    state_qalys         = state_qalys,
    implant_costs       = implant_costs,
    discount_factors    = discount_factors,
    n_cycles = n_cycles, n_implants = n_implants,
    n_samples = n_samples, n_states = n_states,
    start_state = which(state_names == "post_thr") - 1)
  
  total_costs <- totals$total_costs
  total_qalys <- totals$total_qalys
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

# Before calling this function, compile the C++ loop once per session with
#   Rcpp::sourceCpp("optimisation/markov_loop.cpp")