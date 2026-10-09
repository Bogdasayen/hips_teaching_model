# HIPS Model - Monolithic Script
# Based on Fawsitt et al. 2019

# Set seed for reproducibility
set.seed(123456)

# Load libraries
library(readxl)
library(BCEA)

# Set data directory
data_directory <- "data"

# Define global parameters
n_samples <- 1000
n_states <- 4
n_implants <- 4
initial_age <- 65 # cohort starts at age 65
n_cycles <- 30 # 30 one-year cycles

state_names <- c("post_thr", "post_1st_rev", "post_2nd_rev", "dead")
implant_names <- c("cemented", "uncemented", "hybrid", "reverse_hybrid")

# Generate Input Parameters
# Create the matrix to hold all sampled parameters (15 parameters)
input_parameters <- as.data.frame(matrix(nrow = n_samples, ncol = 15))
colnames(input_parameters) <- c(
  paste0("log_rate_1st_revision_", implant_names),
  "log_rate_2nd_revision", "log_rate_higher_revision",
  "cost_revision",
  paste0("state_utility_", state_names),
  paste0("implant_cost_", implant_names)
)

# Log rates of 1st revision (Normal on the log-rate scale), by implant
input_parameters$log_rate_1st_revision_cemented       <- rnorm(n_samples, -5.36497, 0.03)
input_parameters$log_rate_1st_revision_uncemented     <- rnorm(n_samples, -5.80764, 0.03)
input_parameters$log_rate_1st_revision_hybrid         <- rnorm(n_samples, -5.60664, 0.03)
input_parameters$log_rate_1st_revision_reverse_hybrid <- rnorm(n_samples, -5.51946, 0.03)

# Log rates of 2nd and higher revision (same across implants)
input_parameters$log_rate_2nd_revision    <- rnorm(n_samples, -3.14364, 0.07259)
input_parameters$log_rate_higher_revision <- rnorm(n_samples, -2.84900, 0.05300)

# Revision cost (fixed)
input_parameters$cost_revision <- 9004.95

# State utilities (Normal); dead state fixed at zero
input_parameters$state_utility_post_thr     <- rnorm(n_samples, 0.762, 0.004)
input_parameters$state_utility_post_1st_rev <- rnorm(n_samples, 0.575, 0.009)
input_parameters$state_utility_post_2nd_rev <- rnorm(n_samples, 0.575, 0.009)
input_parameters$state_utility_dead         <- 0

# Implant costs (fixed)
input_parameters$implant_cost_cemented       <- 756
input_parameters$implant_cost_uncemented     <- 2386
input_parameters$implant_cost_hybrid         <- 1597
input_parameters$implant_cost_reverse_hybrid <- 1446

# Convenience objects used later in the script
state_utilities <- as.matrix(
  input_parameters[, paste0("state_utility_", state_names)])
implant_costs <- t(as.matrix(
  input_parameters[, paste0("implant_cost_", implant_names)]))

# Read UK lifetables for mortality
uk_lifetables <- read_excel(
  paste0(data_directory, "/hips_input_data.xlsx"), sheet = "uk_lifetables")

# Generate Transition Matrices
# Array: n_cycles x n_implants x n_samples x n_states x n_states
transition_matrices <- array(0,
                             dim = c(n_cycles, n_implants, n_samples, n_states, n_states),
                             dimnames = list(NULL, implant_names, NULL, state_names, state_names))

# Fill transition matrices
for (i_cycle in 1:n_cycles) {
  for (implant_name in implant_names) {
    # Probability of 1st revision (from log hazard rate)
    transition_matrices[i_cycle, implant_name, , "post_thr", "post_1st_rev"] <-
      1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_", implant_name)]))
    
    # Probability of 2nd revision
    transition_matrices[i_cycle, implant_name, , "post_1st_rev", "post_2nd_rev"] <-
      1 - exp(-exp(input_parameters$log_rate_2nd_revision))
    
    # Age-dependent mortality (same for all living states)
    prob_death <- 1 - exp(-uk_lifetables$Males[initial_age + i_cycle - 1])
    transition_matrices[i_cycle, implant_name, , c(1:3), "dead"] <- prob_death
    
    # Probability of staying in state (ensures rows sum to 1)
    for (i_state in 1:n_states) {
      transition_matrices[i_cycle, implant_name, , i_state, i_state] <-
        1 - rowSums(transition_matrices[i_cycle, implant_name, , i_state, -i_state])
    }
  }
}

# Generate State Costs
# Array to store state costs: n_implants x n_samples x n_states
state_costs <- array(NA, dim = c(n_implants, n_samples, n_states),
                     dimnames = list(implant_names, NULL, state_names))

# Cost in post THR state = P(1st revision) * cost of revision
for (implant_name in implant_names) {
  prob_1st_rev <- 1 - exp(-exp(
    input_parameters[, paste0("log_rate_1st_revision_", implant_name)]))
  state_costs[implant_name, , "post_thr"] <-
    input_parameters$cost_revision * prob_1st_rev
}

# Cost in post-1st-revision state = P(2nd revision) * cost of revision
prob_2nd_rev <- 1 - exp(-exp(input_parameters$log_rate_2nd_revision))
state_costs[, , "post_1st_rev"] <- rep(
  input_parameters$cost_revision * prob_2nd_rev, each = n_implants)

# Cost in post-2nd-revision state = P(3rd or higher revision) * cost of revision
prob_higher_rev <- 1 - exp(-exp(input_parameters$log_rate_higher_revision))
state_costs[, , "post_2nd_rev"] <- rep(
  input_parameters$cost_revision * prob_higher_rev, each = n_implants)

# Zero cost in dead state
state_costs[, , "dead"] <- 0

# Run Markov Model and Calculate Net Benefit
# Initialize cohort vectors (everyone starts in post THR)
cohort_vectors <- array(dim = c(n_cycles, n_implants, n_samples, n_states),
                        dimnames = list(NULL, implant_names, NULL, state_names))
cohort_vectors[1, , , ] <- 0
cohort_vectors[1, , , "post_thr"] <- 1

# Arrays for costs and QALYs
cycle_costs <- array(dim = c(n_cycles, n_implants, n_samples))
cycle_qalys <- array(dim = c(n_cycles, n_implants, n_samples))
total_costs <- array(dim = c(n_implants, n_samples))
total_qalys <- array(dim = c(n_implants, n_samples))

# Discount factors (3.5% per year)
discount_factors <- (1 / 1.035)^(0:(n_cycles - 1))

# Main Markov loop
for (i_implant in 1:n_implants) {
  for (i_sample in 1:n_samples) {
    # Markov updating: multiply cohort vector by transition matrix
    for (i_cycle in 2:n_cycles) {
      cohort_vectors[i_cycle, i_implant, i_sample, ] <-
        cohort_vectors[i_cycle - 1, i_implant, i_sample, ] %*%
        transition_matrices[i_cycle - 1, i_implant, i_sample, , ]
    }
    
    # Calculate cycle costs and QALYs
    cycle_costs[, i_implant, i_sample] <-
      cohort_vectors[, i_implant, i_sample, ] %*% state_costs[i_implant, i_sample, ]
    cycle_qalys[, i_implant, i_sample] <-
      cohort_vectors[, i_implant, i_sample, ] %*% state_utilities[i_sample, ]
    
    # Apply discounting and sum
    total_costs[i_implant, i_sample] <- implant_costs[i_implant, i_sample] +
      sum(cycle_costs[, i_implant, i_sample] * discount_factors)
    total_qalys[i_implant, i_sample] <-
      sum(cycle_qalys[, i_implant, i_sample] * discount_factors)
  }
}

# Calculate net monetary benefit (at threshold = 20,000 per QALY)
lambda <- 20000
net_benefit <- total_qalys * lambda - total_costs

# Incremental net benefit relative to the first implant (cemented)
incremental_net_benefit <- net_benefit - matrix(net_benefit[1, ],
                                                nrow = n_implants, ncol = n_samples, byrow = TRUE)

# Summarise Results
# Mean costs, QALYs and net benefit by implant
results <- data.frame(
  implant     = implant_names,
  mean_costs  = rowMeans(total_costs),
  mean_qalys  = rowMeans(total_qalys),
  mean_nb     = rowMeans(net_benefit),
  mean_inb    = rowMeans(incremental_net_benefit),
  prob_ce     = rowMeans(net_benefit == matrix(apply(net_benefit, 2, max),
                                               nrow = n_implants, ncol = n_samples, byrow = TRUE))
)
print(results)

# Cost-effectiveness analysis with BCEA (expects samples x interventions)
hips_bcea <- bcea(e = t(total_qalys), c = t(total_costs),
                  ref = 1, interventions = implant_names, Kmax = 50000)
summary(hips_bcea)
ceac.plot(hips_bcea)