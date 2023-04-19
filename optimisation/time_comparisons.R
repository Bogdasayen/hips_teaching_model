# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom February 2022

# These lines are only needed if using the future.apply package but this didn't give speed beenfit
#library(future.apply)
#plan(multisession) ## Run in parallel on local computer

# Set seed for random number generation
set.seed(2345295)

# Target population is 65 year old males

# Load necessary libraries
library(readxl)
library(BCEA)
library(ggplot2)


# Define global parameters for the data and source code directories
# Data could be separate location than GitHub project if data were sensitive
data_directory <- "data"
code_directory <- "code"
optimisation_directory <- "code/optimisation"

# Load necessary functions
source(paste0(code_directory, "/generate_input_parameters.R"))
source(paste0(code_directory, "/generate_transition_matrices.R"))
source(paste0(code_directory, "/generate_state_costs.R"))
source(paste0(code_directory, "/generate_net_benefit.R"))

# Run optimisations scripts
source(paste0(optimisation_directory, "/generate_transition_matrices_optimised.R"))
source(paste0(optimisation_directory, "/generate_transition_matrices_df.R"))
source(paste0(optimisation_directory, "/generate_net_benefit_lapply.R"))
source(paste0(optimisation_directory, "/generate_net_benefit_lapply_vectorised.R"))
source(paste0(optimisation_directory, "/generate_net_benefit_lapply_vectorised_parallel.R"))
source(paste0(optimisation_directory, "/generate_net_benefit_cpp_partial.R"))
source(paste0(optimisation_directory, "/convert_transition_matrices_to_df.R"))
source(paste0(optimisation_directory, "/generate_net_benefit_cpp_full.R"))
source(paste0(optimisation_directory, "/generate_net_benefit_cpp_full_parallel.R"))
source(paste0(optimisation_directory, "/generate_net_benefit_df.R"))

# Define global parameters that are accessed by all functions
n_samples <- 1000
n_states <- 4
n_implants <- 4
# Age of cohort is used to access correct lifetable
initial_age <- 65
# 30 one-year cycles for patients initially aged 65.
n_cycles <- 30

state_names <- c("post_thr", "post_1st_rev", "post_2nd_rev", "dead")
implant_names <- c("cemented", "uncemented", "hybrid", "reverse_hybrid")

# Generate a matrix of input parameters
# One row for each sample, one column for each parameter
input_parameters <- generate_input_parameters(n_samples)

# Note that some of the improvement comes from an optimised generate_transition_matrices function
# Runtime reduces by about 9%
system.time({
  transition_matrices_slow <- generate_transition_matrices(input_parameters)
})
system.time({
  transition_matrices_fast <- generate_transition_matrices_optimised(input_parameters)
})
# Further reduces runtime by 65%
system.time({
  transition_matrices_df <- generate_transition_matrices_df(input_parameters)
})

# Generate the net benefit using the above input parameters
system.time({
model_output <- generate_net_benefit(input_parameters)
})
# Replace implant loop with lapply reduces by ~30%
system.time({
  mo_lapply <- generate_net_benefit_lapply(input_parameters)
})
# Reduced by ~53% (fastest with no C/C++ or parallelisation)
system.time({
  mo_lapply_vectorised <- generate_net_benefit_lapply_vectorised(input_parameters)
})

# Run lapply in parallel using future.apply package
# Slower if using 1000 or 10000 samples as set-up for parallelisation takes time
# At 100000 got an error that the exports exceeded the limit for 'future' expressions
#system.time({
#  mo_lapply_vectorised_parallel <- generate_net_benefit_lapply_vectorised_parallel(input_parameters)
#})

# Conversion of loop over cycles to C/C++ by ~45%
system.time({
  mo_cpp_partial <- generate_net_benefit_cpp_partial(input_parameters)
})
# Full conversion of Markov loop to C/C++ and uses transition matrices as data frame
# Reduces runtime by 75% (though mostly from transition matrix as data frame)
system.time({
  mo_cpp_full <- generate_net_benefit_cpp_full(input_parameters)
})

# Run lapply of C++ in parallel, but this is slower due to parallelisation setup time
# Again slower that non-parallel case
#system.time({
#  mo_cpp_full_parallel <- generate_net_benefit_cpp_full_parallel(input_parameters)
#})
# Conversion to data frames is an interim step before C/C++ - don't run as VERY slow
#system.time({
 # mo_df <- generate_net_benefit_df(input_parameters)
#})



# Ensure different methods give the same results
rowMeans(model_output$net_benefit)
rowMeans(mo_lapply_vectorised$net_benefit)
rowMeans(mo_lapply_vectorised_parallel$net_benefit)
rowMeans(mo_cpp_partial$net_benefit)
rowMeans(mo_cpp_full$net_benefit)



