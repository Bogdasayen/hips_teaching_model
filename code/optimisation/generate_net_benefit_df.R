# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom February 2022

# Function to use Markov modelling to simulate total costs, QALYs and net benefit

# Changed transition matrices and cohort vectors to data.frame
# This allows for larger numbers of cycles, samples, implants and states
# Extremely slow however as no efficient way to matrix multiply data frames in 
# R so for this (non-C/C++) converts each data frame to a matrix and numeric vector

# Idea is that it makes it easier to move core loop to C/C++

# Packages needed to easily manipulate data frames
require(reshape)
require(dplyr)

generate_net_benefit_df <- function(input_parameters, lambda = 20000) {
  # First generate components needed for simulation
  transition_matrices <- generate_transition_matrices_optimised(input_parameters)
  state_costs <- generate_state_costs(input_parameters)
  state_qalys <- input_parameters[, grepl("state_utility", colnames(input_parameters))]
  
  # Convert transition matrices array to a dataframe
  # Must be in cycle, implant, sample, "from state", "to states" order
  # Also remove dimnames so data.frame can index using numbers (important if going to C/C++)
  transition_matrices_temp <- transition_matrices
  dimnames(transition_matrices_temp) <- list(NULL, NULL, NULL, NULL, state_names)
  # List of data frames for transition from each state
  transition_matrices_temp_df <- list()
  for(i_state in 1:n_states) {
    # Each stores the transition probabilities from i_state
    transition_matrices_temp_df[[i_state]] <- transition_matrices_temp[, , , , i_state]
    # Convert the multidimensional array to a data frame
    transition_matrices_temp_df[[i_state]] <- 
      melt(transition_matrices_temp_df[[i_state]], varnames = c("cycle", "implant", "sample", "from"))
    # Name the state to which you're transiting
    colnames(transition_matrices_temp_df[[i_state]])[5] <- state_names[i_state]
  }
  # Combine the data frames
  transition_matrices_df <- do.call(cbind, transition_matrices_temp_df)
  # Only keep unique columns
  transition_matrices_df <- transition_matrices_df[, unique(colnames(transition_matrices_df))]
  # Sort so that from asce
  transition_matrices_df %>% arrange(cycle, implant, sample, from)
  
  # Check if we've set them up correctly
  #transition_matrices_df[with(transition_matrices_df, cycle == 5 & implant == 3 & sample == 101), c(4:8)]
  #transition_matrices[5, 3, 101, , ]

  
  
  # Store the cohort vectors as a data frame with one row for each cycle, implant and sample
  cohort_vectors <- data.frame(cycle = rep(c(1:n_cycles), n_implants * n_samples),
                               implant = rep(c(1:n_implants), each = n_cycles, n_samples),
                               sample = rep(c(1:n_samples), each = n_cycles * n_implants))
  # One column for each states
  for(i_state in 1:n_states) {cohort_vectors <- cbind(cohort_vectors, rep(0, n_cycles * n_implants * n_samples))}
  colnames(cohort_vectors)[4:(3 + n_states)] <- state_names
  
  # Assume everyone starts in the post_thr state
  cohort_vectors[cohort_vectors$cycle == 1, "post_thr"] <- 1
  # All other proportions start at zero by default when setting up data frame
  
  
  # Implant costs (transpose to keep convention of n_implants, n_samples)
  implant_costs <- t(input_parameters[, grepl("implant_cost_", colnames(input_parameters))])
  
  # Build an array to store the costs and QALYs accrued per cycle
  # One for each cycle, implant, and PSA sample
  # These will be filled in below in the main model code
  # Then discounted and summed to contribute to total costs and total QALYs
  cycle_costs <- array(dim = c(n_cycles, n_implants, n_samples), 
                       dimnames = list(NULL, implant_names, NULL))
  cycle_qalys <- array(dim = c(n_cycles, n_implants, n_samples), 
                       dimnames = list(NULL, implant_names, NULL))
  
  # Build arrays to store the total costs and total QALYs
  # There is one for each treatment and each PSA sample
  # These are filled in below using cycle_costs, 
  # implant_costs, and cycle_qalys
  total_costs <- array(dim = c(n_implants, n_samples), 
                       dimnames = list(implant_names, NULL))
  total_qalys <- array(dim = c(n_implants, n_samples), 
                       dimnames = list(implant_names, NULL))
  
  
  # The remainder of the cohort_vectors will be filled in by Markov updating below
  

  # Main model code
  # Loop over the implants
  
  for (i_implant in 1:n_implants) {
    # Loop over the PSA samples
    for (i_sample in 1:n_samples) {
      # Loop over the cycles
      # Cycle 1 is already defined so only need to update cycles 2:n_cycles
      for (i_cycle in 2:n_cycles) {
        # Markov update
        # Multiply previous cycle's cohort vector by transition matrix
        # i_e_ pi_j = pi_(j-1)*P
        
        # Extremely slow due to type conversions
        cohort_vectors[with(cohort_vectors, cycle == i_cycle & implant == i_implant & sample == i_sample), c(4:7)] <-
          as.data.frame(as.numeric(as.vector(cohort_vectors[with(cohort_vectors, cycle == i_cycle - 1 & implant == i_implant & sample == i_sample), c(4:7)])) %*%
                          as.matrix(transition_matrices_df[with(transition_matrices_df, cycle == i_cycle - 1 & implant == i_implant & sample == i_sample), c(5:8)]))
        
      }
    }
  }
    
  # Original loops to calculate total costs and qalys
  for(i_implant in 1:n_implants) {
    for(i_sample in 1:n_samples) {
      # Now use the cohort vectors to calculate the 
      # total costs for each cycle
      cycle_costs[, i_implant, i_sample] <- 
        cohort_vectors[, i_implant, i_sample, ] %*% state_costs[i_implant, i_sample, ]
      # And total QALYs for each cycle"
      cycle_qalys[, i_implant, i_sample] <- 
        cohort_vectors[, i_implant, i_sample, ] %*% t(state_qalys[i_sample, ])
      
      # Combine the cycle_costs and implant_costs to get total costs
      # Apply the discount factor 
      # (1 in first year, 1_035 in second, 1_035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total_costs[i_implant, i_sample] <- implant_costs[i_implant, i_sample] + 
        cycle_costs[, i_implant, i_sample] %*%
        (1 / 1.035)^rep(c(0:(n_cycles-1)), each = 1)
      
      # Combine the cycle_qalys to get total qalys
      # Apply the discount factor 
      # (1 in first year, 1_035 in second, 1_035^2 in third, and so on)
      # Each year acounts for two cycles so need to repeat the discount values
      total_qalys[i_implant, i_sample] <- cycle_qalys[, i_implant, i_sample]%*%
        (1 / 1.035)^rep(c(0:(n_cycles-1)), each = 1)
    }
  }
  
  # Calculate net benefit and incremental net benefit at
  # willingness-to-pay that was supplied
  net_benefit <- total_qalys * lambda - total_costs
  incremental_net_benefit <- net_benefit - net_benefit[1, ]
  
  return(list("total_costs" = total_costs,
              "total_qalys" = total_qalys,
              "net_benefit" = net_benefit,
              "incremental_net_benefit" = incremental_net_benefit))
}