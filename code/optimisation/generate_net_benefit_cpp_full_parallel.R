# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom February 2022

# Function to use Markov modelling to simulate total costs, QALYs and net benefit

# Optimised by replacing the loop over cycles as C++
# Loop over treatments and samples remains in R

# Uses the Rcpp package
require(Rcpp)
# Compiles the C++ file for the Markov loop
Rcpp::sourceCpp(file = paste0(optimisation_directory, "/rcpp_loop_full.cpp"))


generate_net_benefit_cpp_full_parallel <- function(input_parameters, lambda = 20000) {
  # First generate components needed for simulation
  transition_matrices_df <- generate_transition_matrices_df(input_parameters)
  state_costs <- generate_state_costs(input_parameters)
  state_qalys <- input_parameters[, grepl("state_utility", colnames(input_parameters))]


  # Store the cohort vectors as a data frame with one row for each cycle, implant and sample
  cohort_vectors <- data.frame(cycle = rep(c(1:n_cycles), n_implants * n_samples),
                               implant = rep(c(1:n_implants), each = n_cycles, n_samples),
                               sample = rep(c(1:n_samples), each = n_cycles * n_implants))
  # One column for each states
  for(i_state in 1:n_states) {cohort_vectors <- cbind(cohort_vectors, rep(0, n_cycles * n_implants * n_samples))}
  colnames(cohort_vectors)[4:(3 + n_states)] <- state_names
  
  # Sort 
  cohort_vectors <- cohort_vectors %>% arrange(cycle, implant, sample)
  
  
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
  
  # Pre-calculate the discount vector to reduce runtime
  discount_vector <- (1 / 1.035)^rep(c(0:(n_cycles-1)), each = 1)
  state_qalys <- as.matrix(state_qalys)
  # Main model code
  # Use lapply instead of loop
  
  # Function currently requires conversion to matrices and returns a matrix, not a dataframe
  cohort_vectors <- 
    (rcpp_loop_full(cohort_vectors_in = as.matrix(cohort_vectors), 
                   transition_matrices = as.matrix(transition_matrices_df),
                   n_cycles = n_cycles, n_implants = n_implants, n_samples = n_samples, n_states = n_states))
  
  
  future_lapply(c(1:n_implants), function(i_implant){
    # Pre-index to reduce runtime
    #cohort_vectors_tr <- cohort_vectors[i_implant ] 
    cycle_costs_tr <- cycle_costs[, i_implant, ]
    cycle_qalys_tr <- cycle_qalys[, i_implant, ]
    treatment_costs_tr <- implant_costs[i_implant, ]
    total_costs_tr <- total_costs[i_implant, ]
    total_qalys_tr <- total_qalys[i_implant, ]
    state_costs_tr <- state_costs[i_implant, , ]
    # In this case state qalys are the same for all treatments/implants 
    # but in general allow this for optimization
    state_qalys_tr <- state_qalys
    
    for(i_sample in 1:n_samples) {
      cohort_vectors_tr_sample <- cohort_vectors[c(0:(n_cycles-1)) * (n_implants * n_samples) +
                                                     (i_implant - 1) * n_samples + i_sample, c(4:(3+n_states))]
      # Use cohort vectors to calculate cycle costs and qalys
      cycle_costs_tr[, i_sample] <- cohort_vectors_tr_sample[,  ] %*% state_costs_tr[i_sample, ]
      cycle_qalys_tr[, i_sample] <- cohort_vectors_tr_sample[, ] %*% state_qalys_tr[i_sample, ]
      # Sum  and discount to get total costs and qalys
      # Add implant costs to total costs
      total_costs_tr[i_sample] <- treatment_costs_tr[i_sample] + cycle_costs_tr[, i_sample] %*% discount_vector
      total_qalys_tr[i_sample] <- cycle_qalys_tr[, i_sample] %*% discount_vector
    }
    
    return(list(total_qalys = total_qalys_tr, total_costs = total_costs_tr))
    
  }) -> output_list # End lapply
  
  names(output_list) <- implant_names 
  
  # Reconvert result to a matrix
  total_costs <- sapply(implant_names, function(implant_name) {total_costs[implant_name, ] <- output_list[[implant_name]]$total_costs})
  total_qalys <- sapply(implant_names, function(implant_name) {total_qalys[implant_name, ] <- output_list[[implant_name]]$total_qalys})
  # sapply inverts the matrices to uninvert
  total_costs <- t(total_costs)
  total_qalys <- t(total_qalys)
  
  # Calculate net benefit and incremental net benefit at
  # willingness-to-pay that was supplied
  net_benefit <- total_qalys * lambda - total_costs
  incremental_net_benefit <- net_benefit - net_benefit[1, ]
  
  return(list("total_costs" = total_costs,
              "total_qalys" = total_qalys,
              "net_benefit" = net_benefit,
              "incremental_net_benefit" = incremental_net_benefit))
}