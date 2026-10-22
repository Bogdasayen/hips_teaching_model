# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom April 2023
# Modified by Qian Xin October 2026

# Set seed for random number generation
set.seed(2345295)

# Target population is 65 year old males

# Load necessary libraries
# install.packages("BCEA")
# install.packages("devtools")
# install.packages("roxygen2")
# install.packages("readxl")
# install.packages("ggplot2")
library(readxl)
library(BCEA)

# Load necessary functions from the package
devtools::load_all()
# Update the documentation (if necessary)
# roxygen2::roxygenise()

# Collect all structural settings in a single list
model_settings <- list(
  n_samples      = 1000,
  n_cycles       = 30,        # 30 one-year cycles
  cycle_length   = 1,         # years
  initial_age    = 65,        # cohort starts at age 65
  discount_rate  = 0.035,     # 3.5% per year
  state_names    = c("post_thr", "post_1st_rev", "post_2nd_rev", "dead"),
  implant_names  = c("cemented", "uncemented", "hybrid", "reverse_hybrid"),
  data_directory = "data"    
)

# Run the model: each function takes what it needs as arguments
input_parameters    <- generate_input_parameters(model_settings)
transition_matrices <- generate_transition_matrices(input_parameters, model_settings)
state_costs         <- generate_state_costs(input_parameters, model_settings)
state_qalys         <- generate_state_qalys(input_parameters, model_settings)

model_output <- generate_net_benefit(transition_matrices, state_costs,
                                     state_qalys, input_parameters,
                                     model_settings, lambda = 20000)

