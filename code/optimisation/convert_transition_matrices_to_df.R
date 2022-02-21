# Function to convert multidimensional array of transition matrices to a dataframe
# Howard Thom 16Feb2022

convert_transition_matrices_to_df <- function(transition_matrices) {
  
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
  # Sort 
  transition_matrices_df  <- transition_matrices_df %>% arrange(cycle, implant, sample, from)
  
  return(transition_matrices_df)
}