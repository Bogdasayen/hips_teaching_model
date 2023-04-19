# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom February 2022

# Very simple optimization using the rep() function to avoid double loops
# Less readable and only about 10% faster!

generate_transition_matrices_optimised <- function(input_parameters) {
  # One 4x4 transition matrix for each cycle, implant_name and sample.
  transition_matrices <- array(0, dim = c(n_cycles, n_implants, n_samples, n_states, n_states),
                               dimnames = list(NULL, implant_names, NULL, state_names, state_names))
  
  # UK lifetables are fixed so are not part of the input parameters matrix
  # This saves memory of an n_sample*100 (i.e., number of ages) matrix
  uk_lifetables <- read_excel(paste0(data_directory, "/hips_input_data.xlsx"), sheet = "uk_lifetables")
  
  
  for(implant_name in implant_names) {
    # Probabilities are the same for every cycle
    transition_matrices[, implant_name, , "post_thr", "post_1st_rev"] <- 
      rep(1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_", implant_name)])), each = n_cycles)
    transition_matrices[, implant_name, , "post_1st_rev", "post_2nd_rev"] <- 
      rep(1 - exp(-exp(input_parameters$log_rate_2nd_revision)), each = n_cycles)
  }
  
  for(i_cycle in 1:n_cycles) {
    # Annual probability of death the same for all samples, implants and states (except dead itself)
    transition_matrices[i_cycle, , , -n_states, "dead"] <-
      1 - exp(-with(uk_lifetables, Males[Age = initial_age + i_cycle - 1]))
  }

  # Ensure remaining patients stay in state
  # Sum probabilities of transitions to other states
  # This ensures probabilities sum to 1.
  # Same action repeated for each cycle, sample and implant so use an apply statement  
  for(i_state in 1:length(state_names)) {
    transition_matrices[, , , i_state, i_state] <- 1 - 
      apply(transition_matrices[, , , i_state, -i_state], c(1, 2, 3), sum, na.rm=TRUE)
  } # End loop over states
  
  
  # At this point do a test and throw and error if it fails
  if(prod(rowSums(transition_matrices[sample(1:n_cycles, 1),
                              sample(implant_names, 1),
                              sample(1:n_samples, 1), , ]) == 1) != 1) {
    stop("Rows must sum to 1!")
  }
  
  return(transition_matrices)
}