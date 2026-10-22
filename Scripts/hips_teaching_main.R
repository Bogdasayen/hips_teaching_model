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
# Summary results
# Mean costs, QALYs and net benefit by implant, and the probability that each
# implant has the highest net benefit across the probabilistic samples
results_table <- data.frame(
  implant    = model_settings$implant_names,
  mean_costs = rowMeans(model_output$total_costs),
  mean_qalys = rowMeans(model_output$total_qalys),
  mean_nb    = rowMeans(model_output$net_benefit),
  mean_inb   = rowMeans(model_output$incremental_net_benefit),
  prob_ce    = rowMeans(model_output$net_benefit ==
                          matrix(apply(model_output$net_benefit, 2, max),
                                 nrow = length(model_settings$implant_names),
                                 ncol = model_settings$n_samples, byrow = TRUE)))
print(results_table, row.names = FALSE)

# ICERs relative to the cemented implant
incremental_costs <- sweep(model_output$total_costs, 2, model_output$total_costs[1, ])
incremental_qalys <- sweep(model_output$total_qalys, 2, model_output$total_qalys[1, ])
icer <- rowMeans(incremental_costs) / rowMeans(incremental_qalys)
print(round(icer[-1]))

# Analyses with BCEA
# BCEA expects samples in rows and interventions in columns, so transpose
hips_bcea <- bcea(e = t(model_output$total_qalys),
                  c = t(model_output$total_costs),
                  ref = 1, interventions = model_settings$implant_names,
                  Kmax = 50000)

# Summary at the chosen threshold
summary(hips_bcea, wtp = lambda)

# Cost-effectiveness plane
ceplane.plot(hips_bcea, wtp = lambda, graph = "ggplot",
             xlim = c(-0.05, 0.05), ylim = c(-1800, 1800))

# Cost-effectiveness acceptability curves for all four implants
hips_multi_ce <- multi.ce(hips_bcea)
ceac.plot(hips_multi_ce, graph = "ggplot", pos = c(0, 0.5))


