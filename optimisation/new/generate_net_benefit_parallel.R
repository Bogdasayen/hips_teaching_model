# HIPS Teaching model based on published data in
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Qian Xin October 2026

#'
#' Parallelisation does not speed up a single run of the model; it runs
#' several independent, complete evaluations at once on separate cores. This
#' function runs \code{n_runs} full evaluations (sampling, transition
#' matrices, state values, net benefit), each with its own random draws, and
#' returns a list of their outputs.
#'
#' The gain depends on the number of cores. The approach pays off on a high-performance
#' computing cluster, where tens or hundreds of cores are available and
#' the setup cost is negligible against the work.
#' 
#' \code{plan(multisession)} starts separate R processes and works on
#' Windows, macOS and Linux. One core is left free. Workers are given
#' independent, reproducible random streams by \code{future.seed}.
#'
#' @param n_runs Number of complete model evaluations.
#' @param model_settings A list of model settings, as in the main script.
#' @param net_benefit_fun The net benefit implementation to use in each run;
#'   default \code{generate_net_benefit}.
#' @param transition_fun The transition matrix implementation to use;
#'   default \code{generate_transition_matrices}.
#' @param lambda Willingness-to-pay threshold; default 25000.
#' @param parallel If \code{FALSE}, runs sequentially with \code{lapply}
#'   (used for the comparison).
#' @param seed Seed passed to \code{future.seed} (parallel) or
#'   \code{set.seed} (sequential).
#'
#' @return A list of length \code{n_runs}, each element the list returned by
#'   the net benefit function.
#'
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @export
run_model_parallel <- function(n_runs, model_settings,
                               net_benefit_fun = generate_net_benefit,
                               transition_fun = generate_transition_matrices,
                               lambda = 25000, parallel = TRUE, seed = 123456) {
  
  # One complete evaluation of the model
  run_one <- function(i) {
    input_parameters    <- generate_input_parameters(model_settings)
    transition_matrices <- transition_fun(input_parameters, model_settings)
    state_costs         <- generate_state_costs(input_parameters, model_settings)
    state_qalys         <- generate_state_qalys(input_parameters, model_settings)
    net_benefit_fun(transition_matrices, state_costs, state_qalys,
                    input_parameters, model_settings, lambda = lambda)
  }
  
  if (!parallel) {
    set.seed(seed)
    return(lapply(seq_len(n_runs), run_one))
  }
  
  future::plan(future::multisession,
               workers = max(1, parallel::detectCores() - 1))
  on.exit(future::plan(future::sequential), add = TRUE)
  future.apply::future_lapply(seq_len(n_runs), run_one, future.seed = seed)
}
