# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom April 2023

# Set seed for random number generation
set.seed(2345295)

# Target population is 65 year old males

# Load necessary libraries
# Install the latest version of BCEA from github
# devtools::install_github("n8thangreen/BCEA")
# install.packages("devtools")
# install.packages("roxygen2")
# install.packages("readxl")
# install.packages("ggplot2")

library(readxl)
library(BCEA)
library(ggplot2)

# Load necessary functions from the package
devtools::load_all()
# Update the documentation (if necessary)
# roxygen2::roxygenise()

# Define global parameters for the data directory
# Data could be separate location than GitHub project if data were sensitive
data_directory <- "data"



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

# Generate the net benefit using the above input parameters
model_output <- generate_net_benefit(input_parameters)

# Average costs
rowMeans(model_output$total_costs)
# Average QALYS
rowMeans(model_output$total_qalys)

# Calculate an ICER without BCEA
incremental_costs <- t(model_output$total_costs) - t(model_output$total_costs)[, 1]
incremental_qalys <- t(model_output$total_qalys) - t(model_output$total_qalys)[, 1]
ICER <- colMeans(incremental_costs) / colMeans(incremental_qalys)


# Now use the BCEA package to analyse the results___
# Note costs and QALYs need to be transposed in this example for BCEA to run
hips_bcea <- bcea(e = t(model_output$total_qalys), 
                  c = t(model_output$total_costs), ref = 1, 
                  interventions = implant_names) 
#get summary statistics
summary(hips_bcea, wtp = 20000)
#plot the cost-effectiveness plane
ceplane.plot(hips_bcea, wtp = 20000, xlim = c(-0.05, 0.05), ylim = c(-1800, 1800))
#plot a CEAC
hips_multi_ce <- multi.ce(hips_bcea)
ceac.plot(hips_multi_ce, graph = "ggplot",
          line = list(color = c("red", "green", "blue", "orange")),
        pos = c(0, 0.50))

# Rank proportions EVPPI/EVPI to identify most influential parameters
# Used lambda = 30000 as little uncertainty at lambda = 20000
hips_inp <- createInputs(inputs = input_parameters)
info.rank(inp = hips_inp,
          hips_bcea, xlim = c(0, 1), wtp = 30000, howManyPars = NA, graph = "base",
          mai = c(1.36, 2.2, 1, 1)) # Adjusted the margins to keep param names on screen


##################################################################################
## VoI analysis - takes a few minutes ############################################
##################################################################################

evpi_table <- matrix(nrow = 4, ncol = 2)
rownames(evpi_table) <- c("Total", "1st revision probabilities", "2nd and higher revision probabilities",
                          "Utilities")
colnames(evpi_table) <- c("Per person", "Population")

# Assume implant decision remains relevant for at least 10 years
# And that there are 160000 THR per year
# https://www.njrcentre.org.uk/njrcentre/Patients/Joint-replacement-statistics#:~:text=In%20England%20and%20Wales%20there,the%20practice%20is%20growing%20rapidly.
# This is regardless of age and gender but assume EVPI same for each category
technology_horizon <- 10
discounted_population_size <- sum((1/1.035)^(0:(technology_horizon - 1))) * 160000 

# Total EVPI
evpi_table["Total", c("Per person", "Population")] <-  hips_bcea$evi[201] * c(1, discounted_population_size)


# Log rate of first revision (same EVPPI as if calculating probability)
# Use GP instead of GAM as 4 parameters
evppi_gp_1st_revision <- evppi(he = hips_bcea,
                               param_idx = c("log_rate_1st_revision_cemented",
                                             "log_rate_1st_revision_uncemented",    
                                             "log_rate_1st_revision_hybrid",        
                                             "log_rate_1st_revision_reverse_hybrid"), 
                               input = input_parameters,  method = 'gp')
evpi_table["1st revision probabilities", c("Per person", "Population")] <- evppi_gam_1st_revision$evppi[201] * c(1, discounted_population_size)


# Log rate 2nd and higher revision probabilities
evppi_gam_2nd_and_higher_revision <- evppi(he = hips_bcea,
                                           param_idx = c("log_rate_2nd_revision",
                                                         "log_rate_higher_revision"), 
                                           input = input_parameters, method = 'gam')
evpi_table["2nd and higher revision probabilities", c("Per person", "Population")] <- evppi_gam_2nd_and_higher_revision$evppi[201] * c(1, discounted_population_size)

# Utilities
evppi_gam_utilities <- evppi(he = hips_bcea, 
                             param_idx = c( "state_utility_post_thr",              
                                            "state_utility_post_1st_rev",          
                                            "state_utility_post_2nd_rev"), 
                             input = input_parameters,  method = 'gam')
evpi_table["Utilities", c("Per person", "Population")] <- evppi_gam_utilities$evppi[201] * c(1, discounted_population_size)




