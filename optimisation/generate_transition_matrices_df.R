# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom February 2022

# Rewritten to output a dataframe in required format by cpp_full

generate_transition_matrices_df <- function(input_parameters) {
  transition_matrices_df <-  data.frame(cycle = rep(c(1:n_cycles), each = n_states, n_implants * n_samples),
                                        implant = rep(c(1:n_implants), each = n_states * n_cycles, n_samples),
                                        sample = rep(c(1:n_samples), each = n_states * n_cycles * n_implants),
                                        from = rep(c(1:n_states), n_cycles * n_implants * n_samples))
  for(i_state in 1:n_states) {
    transition_matrices_df[[state_names[i_state]]] <- 0
  }
  
  # UK lifetables are fixed so are not part of the input parameters matrix
  # This saves memory of an n_sample*100 (i.e., number of ages) matrix
  uk_lifetables <- read_excel(paste0(data_directory, "/hips_input_data.xlsx"), sheet = "uk_lifetables")
  
  
  for(i_implant in 1:n_implants) {
    # Probabilities are the same for every cycle
    transition_matrices_df[with(transition_matrices_df, implant == i_implant & from == 1), "post_1st_rev"] <- 
      rep(1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_", implant_names[i_implant])])), each = n_cycles)
    transition_matrices_df[with(transition_matrices_df, implant == i_implant & from == 2), "post_2nd_rev"] <- 
      rep(1 - exp(-exp(input_parameters$log_rate_2nd_revision)), each = n_cycles)
  }

  for(i_cycle in 1:n_cycles) {
    # Annual probability of death the same for all samples, implants and states (except dead itself)
      transition_matrices_df[with(transition_matrices_df, cycle == i_cycle & from != 4), "dead"] <-
      1 - exp(-with(uk_lifetables, Males[Age = initial_age + i_cycle - 1]))
  }

  # Ensure remaining patients stay in state
  # Sum probabilities of transitions to other states
  # This ensures probabilities sum to 1.
  # Same action repeated for each cycle, sample and implant so use an apply statement  
  for(i_state in 1:length(state_names)) {
    transition_matrices_df[with(transition_matrices_df, from == i_state), state_names[i_state]] <- 1 - 
      rowSums(transition_matrices_df[with(transition_matrices_df, from == i_state), state_names[-i_state]])
  } # End loop over states
  
  # Sort the transitions
  transition_matrices_df  <- transition_matrices_df %>% arrange(cycle, implant, sample, from)
  
  
  # At this point do a test and throw and error if it fails
  if(sum(rowSums(transition_matrices_df[, c(5:8)])) != dim(transition_matrices_df)[1]) {
    stop("Rows must sum to 1!")
  }
  
  # Could check if we've set them up correctly
  #transition_matrices_df[with(transition_matrices_df, cycle == 5 & implant == 3 & sample == 101), c(4:8)]
  #transition_matrices[5, 3, 101, , ]
  
  # Can be accessed by indexing the correct row
  #i_cycle <- 5; i_implant <- 3; i_sample <- 101
  #transition_matrices_df[(i_cycle - 1) * (n_implants * n_samples * n_states) +
  #                         (i_implant - 1) * (n_samples * n_states) +
  #                         (i_sample - 1) * n_states + c(1:4), ]
  
  
  return(transition_matrices_df)
}