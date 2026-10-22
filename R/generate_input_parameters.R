# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom April 2023
# Modified by Qian Xin October 2026

#' Generate input parameters for the Markov model
#'
#' Reads the model inputs from the Excel workbook and draws the probabilistic
#' samples used by the rest of the model. The values are defined in the named
#' sheets of \code{hips_input_data.xlsx}: \code{log_rates_1st_revision},
#' \code{other_log_rates}, \code{costs}, and \code{state_utilities}.
#' Revision rates are sampled on the log-rate scale and converted to
#' probabilities later, in \code{\link{generate_transition_matrices}}.
#'
#' @param model_settings A list of model settings. This function uses
#'   \code{n_samples} (number of probabilistic samples), \code{state_names}
#'   and \code{implant_names} (character vectors naming the health states and
#'   implants), and \code{data_directory} (path to the folder containing
#'   \code{hips_input_data.xlsx}).
#'
#' @return A data frame with one row per probabilistic sample and one column
#'   per parameter: the log rates of first revision for each implant, the log
#'   rates of second and higher revision, the (fixed) cost of revision, the
#'   health-state utilities, and the (fixed) cost of each implant.
#'
#' @examples
#' model_settings <- list(
#'   n_samples      = 10,
#'   state_names    = c("post_thr", "post_1st_rev", "post_2nd_rev", "dead"),
#'   implant_names  = c("cemented", "uncemented", "hybrid", "reverse_hybrid"),
#'   data_directory = "data")
#' input_parameters <- generate_input_parameters(model_settings)
#' head(input_parameters)
#'
#' @importFrom readxl read_excel
#' @importFrom stats rnorm
#' @export
#' 
generate_input_parameters <- function(model_settings) {
  
  # Extract the settings this function needs
  n_samples      <- model_settings$n_samples
  state_names    <- model_settings$state_names
  implant_names  <- model_settings$implant_names
  data_directory <- model_settings$data_directory
  
  # Derive dimensions 
  n_states   <- length(state_names)
  n_implants <- length(implant_names)
  
  # One row per sample, one column per parameter
  input_parameters <- as.data.frame(
    matrix(nrow = n_samples, ncol = (n_implants + 2 + 1 + n_states + n_implants)))
  colnames(input_parameters) <- c(
    paste0("log_rate_1st_revision_", implant_names),
    "log_rate_2nd_revision", "log_rate_higher_revision",
    "cost_revision",
    paste0("state_utility_", state_names),
    paste0("implant_cost_", implant_names))
  
  # Read inputs from the named sheets of the workbook
  input_file <- paste0(data_directory, "/hips_input_data.xlsx")
  log_rates_1st_revision <- read_excel(input_file, sheet = "log_rates_1st_revision")
  other_log_rates        <- read_excel(input_file, sheet = "other_log_rates")
  costs                  <- read_excel(input_file, sheet = "costs")
  state_utilities        <- read_excel(input_file, sheet = "state_utilities")
  
  # First-revision log rates are sampled per implant on the log-rate scale
  for (implant_name in implant_names) {
    input_parameters[, paste0("log_rate_1st_revision_", implant_name)] <-
      with(log_rates_1st_revision,
           rnorm(n_samples, mean = Mean[Implant == implant_name],
                 sd   = SE[Implant == implant_name]))
  }
  
  # Second and higher revision log rates
  input_parameters$log_rate_2nd_revision <-
    with(other_log_rates,
         rnorm(n_samples, mean = Value[Parameter == "log_hazard_2nd_revision"],
               sd   = Value[Parameter == "log_hazard_2nd_revision_se"]))
  input_parameters$log_rate_higher_revision <-
    with(other_log_rates,
         rnorm(n_samples, mean = Value[Parameter == "log_hazard_higher_revision"],
               sd   = Value[Parameter == "log_hazard_higher_revision_se"]))
  
  # Costs are fixed: the single value is recycled across all samples
  for (implant_name in implant_names) {
    input_parameters[, paste0("implant_cost_", implant_name)] <-
      with(costs, Cost[Resource == paste0("cost_", implant_name)])
  }
  input_parameters$cost_revision <- with(costs, Cost[Resource == "cost_revision"])
  
  # State utilities are sampled for the living states; dead state has zero utility
  for (state_name in setdiff(state_names, "dead")) {
    input_parameters[, paste0("state_utility_", state_name)] <-
      with(state_utilities,
           rnorm(n_samples, mean = Mean[State == state_name],
                 sd   = SE[State == state_name]))
  }
  input_parameters$state_utility_dead <- 0
  
  return(input_parameters)
}