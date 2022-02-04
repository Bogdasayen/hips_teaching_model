# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom February 2022

# Target population is 60 year old males

library(readxl)
library(BCEA)

# Define global parameters for the data and source code directories
# Data could be separate location than GitHub project if data were sensitive
data_directory <- "data"
code_directory <- "code"

# Load necessary functions
source(paste0(code_directory, "/generate_input_parameters.R"))
source(paste0(code_directory, "/generate_transition_matrices.R"))
source(paste0(code_directory, "/generate_state_costs.R"))
source(paste0(code_directory, "/generate_net_benefit.R"))

# Define global parameters that are accessed by all functions
n_samples <- 1000
n_states <- 4
n_implants <- 4
# Age of cohort is used to access correct lifetable
initial_age <- 60
# 30 one-year cycles for patients initially aged 60.
n_cycles <- 30

state_names <- c("post_thr", "post_1st_rev", "post_2nd_rev", "dead")
implant_names <- c("cemented", "uncemented", "hybrid", "reverse_hybrid")

# Generate a matrix of input parameters
# One row for each sample, one column for each parameter
input_parameters <- generate_input_parameters(n_samples)

# Generate the net benefit using the above input parameters
model_output <- generate_net_benefit(input_parameters)


# Now use the BCEA package to analyse the results___
# Note costs and QALYs need to be transposed in this example for BCEA to run
hips_bcea <- bcea(e = t(model_output$total_qalys), c = t(model_output$total_costs), ref = 1, interventions = implant_names) 
#get summary statistics
summary(hips_bcea, wtp = 20000)
#plot the cost-effectiveness plane
ceplane.plot(hips_bcea, wtp = 20000)
#plot a CEAC
hips_multi_ce <- multi.ce(hips_bcea)
mce.plot(hips_multi_ce, pos = c(1,0))

# Conduct EVPPI analysis
evpi_table <- matrix(nrow = 4, ncol = 2)
rownames(evpi_table) <- c("Total", "1st revision probabilities", "2nd and higher revision probabilities",
                          "Utilities")
colnames(evpi_table) <- c("Per person", "Population")

# Assume implant decision remains relevant for at least 10 years
# And that there are 10000 revisions per year
technology_horizon <- 10
discounted_population_size <- sum((1/1.035)^(0:(technology_horizon - 1))) * 10000 

# Total EVPI
evpi_table["Total", c("Per person", "Population")] <-  hips_bcea$evi[201] * c(1, discounted_population_size)


# Log rate of first revision (same EVPPI as if calculating probability)
# Use GP instead of GAM as 4 parameters
evppi_gp_1st_revision <- evppi(parameter =  
                               c("log_rate_1st_revision_cemented",
                                 "log_rate_1st_revision_uncemented",    
                                 "log_rate_1st_revision_hybrid",        
                                 "log_rate_1st_revision_reverse_hybrid"), 
                             input = input_parameters, he = hips_bcea, method = 'gp')
evpi_table["1st revision probabilities", c("Per person", "Population")] <- evppi_gam_1st_revision$evppi[201] * c(1, discounted_population_size)

input_parameters$log_rate_2nd_revision
# Log rate 2nd and higher revision probabilities
evppi_gam_2nd_and_higher_revision <- evppi(parameter =  
                                  c("log_rate_2nd_revision",
                                    "log_rate_higher_revision"), 
                                input = input_parameters, he = hips_bcea, method = 'gam')
evpi_table["2nd and higher revision probabilities", c("Per person", "Population")] <- evppi_gam_2nd_and_higher_revision$evppi[201] * c(1, discounted_population_size)

# Utilities
evppi_gam_utilities <- evppi(parameter =  
                               c( "state_utility_post_thr",              
                                  "state_utility_post_1st_rev",          
                                  "state_utility_post_2nd_rev"), 
                             input = input_parameters, he = hips_bcea, method = 'gam')
evpi_table["Utilities", c("Per person", "Population")] <- evppi_gam_utilities$evppi[201] * c(1, discounted_population_size)


# Rank proportions EVPPI/EVPI to identify most influential parameters
info.rank(parameter = colnames(input_parameters), input = input_parameters, 
          hips_bcea, xlim = c(0,1), wtp = 20000, howManyPars = 10)



