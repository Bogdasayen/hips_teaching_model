% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_transition_matrices.R
\name{generate_transition_matrices}
\alias{generate_transition_matrices}
\title{Generate transition matrices}
\usage{
generate_transition_matrices(input_parameters)
}
\arguments{
\item{input_parameters}{Matrix with row for each sample and column 
for each parameter, with values samples for the model input parameters}
}
\value{
Array of n_cycles by n_implants by n_samples by n_states by n_states, 
with values equal to transition probabilities. These correspond to a n_state
by n_state transition matrix for each cycle, implant and sample.
}
\description{
Function to generate transition matrices. Used internall by 
generate_net_benefit() function
}
\examples{
# First sample the input parameters
input_parameters <- generate_input_parameters(n_samples)

transition_matrices <- generate_transition_matrices(input_parameters)

# The transition matrix for first cycle, first implant and first sample
transition_matrices[1, 1, 1, , ]
#' # The transition matrix for first cycle, third implant and second sample
transition_matrices[1, 3, 2, , ]
# Sampled transitions from the first state for first cycle and third implant
transition_matrices[1, 3, , 1, , ]
}
