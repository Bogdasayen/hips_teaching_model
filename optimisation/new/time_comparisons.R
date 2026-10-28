# HIPS Teaching model based on published data in 
# Fawsitt 2019 https://pubmed.ncbi.nlm.nih.gov/30832968/
# Howard Thom February 2022
# Modified by Qian Xin October 2026

# Benchmarks the optimisation techniques described in the tutorial and writes
# the speed tables to results/benchmark.csv.
#
# Run from the package root:  source("optimisation/new/time_comparisons.R")
# Needs: microbenchmark, Rcpp (and a C++ compiler), future, future.apply

# Target population is 65 year old males

library(microbenchmark)
library(readxl)
devtools::load_all()                       # baseline functions in R/

# The optimised variants live outside the package; change this one line if
# they move
optimisation_dir <- "optimisation/new"
source(file.path(optimisation_dir, "generate_transition_matrices_vectorised.R"))
source(file.path(optimisation_dir, "generate_transition_matrices_cpp.R"))
source(file.path(optimisation_dir, "generate_net_benefit_vectorised.R"))
source(file.path(optimisation_dir, "generate_net_benefit_cpp.R"))
source(file.path(optimisation_dir, "generate_net_benefit_parallel.R"))
Rcpp::sourceCpp(file.path(optimisation_dir, "transition_matrices.cpp"))
Rcpp::sourceCpp(file.path(optimisation_dir, "markov_loop.cpp"))

model_settings <- list(
  n_samples = 1000, n_cycles = 30, cycle_length = 1, initial_age = 65,
  discount_rate = 0.035,
  state_names    = c("post_thr", "post_1st_rev", "post_2nd_rev", "dead"),
  implant_names  = c("cemented", "uncemented", "hybrid", "reverse_hybrid"),
  data_directory = "data")

n_times <- 100 # number of times each function is run
time_it <- function(expr, times = n_times) {
  e <- substitute(expr); env <- parent.frame()
  med <- suppressWarnings(median(microbenchmark(eval(e, env), times = times)$time)) / 1e9
  # Functions faster than a millisecond are below the clock's resolution on
  # some systems; time 100 calls together and divide
  if (med < 1e-3) {
    med <- median(microbenchmark(for (k in 1:100) eval(e, env), times = times)$time) / 1e9 / 100
  }
  med
}


# The functions of the model, by stage and technique. Stages with no entry
# for a technique use the baseline, because the technique does not apply:
# state costs and QALYs are already single vectorised assignments, and input
# parameters is dominated by reading the workbook.
implementations <- list(
  input_parameters    = list(baseline = generate_input_parameters),
  transition_matrices = list(baseline   = generate_transition_matrices,
                             vectorised = generate_transition_matrices_vectorised,
                             cpp        = generate_transition_matrices_cpp),
  state_costs         = list(baseline = generate_state_costs),
  state_qalys         = list(baseline = generate_state_qalys),
  net_benefit         = list(baseline   = generate_net_benefit,
                             vectorised = generate_net_benefit_vectorised,
                             cpp        = generate_net_benefit_cpp))

# Pick the implementation of a stage for a technique, falling back to baseline
pick <- function(stage, technique) {
  impl <- implementations[[stage]]
  if (!is.null(impl[[technique]])) impl[[technique]] else impl$baseline
}

# Run the whole model with one technique applied wherever it exists
run_model <- function(technique, model_settings) {
  input_parameters    <- pick("input_parameters",    technique)(model_settings)
  transition_matrices <- pick("transition_matrices", technique)(input_parameters, model_settings)
  state_costs         <- pick("state_costs",         technique)(input_parameters, model_settings)
  state_qalys         <- pick("state_qalys",         technique)(input_parameters, model_settings)
  pick("net_benefit", technique)(transition_matrices, state_costs, state_qalys,
                                 input_parameters, model_settings)
}

#  Check every implementation gives the same result as the baseline 
set.seed(123456); reference <- run_model("baseline", model_settings)
for (technique in c("vectorised", "cpp")) {
  set.seed(123456); out <- run_model(technique, model_settings)
  stopifnot(isTRUE(all.equal(out, reference, tolerance = 1e-10)))
}

# Shared inputs for the per-function timings
set.seed(123456)
input_parameters    <- generate_input_parameters(model_settings)
transition_matrices <- generate_transition_matrices(input_parameters, model_settings)
state_costs         <- generate_state_costs(input_parameters, model_settings)
state_qalys         <- generate_state_qalys(input_parameters, model_settings)

# each function, each technique 
call_stage <- function(stage, fun) {
  switch(stage,
         input_parameters    = fun(model_settings),
         transition_matrices = fun(input_parameters, model_settings),
         state_costs         = fun(input_parameters, model_settings),
         state_qalys         = fun(input_parameters, model_settings),
         net_benefit         = fun(transition_matrices, state_costs, state_qalys,
                                   input_parameters, model_settings))
}
part1 <- do.call(rbind, lapply(names(implementations), function(stage) {
  do.call(rbind, lapply(names(implementations[[stage]]), function(technique) {
    fun <- implementations[[stage]][[technique]]
    data.frame(part = "1. Individual functions", stage = stage, technique = technique,
               median_seconds = time_it(call_stage(stage, fun)))
  }))
}))
part1$reduction_percent <- NA
for (stage in unique(part1$stage)) {
  rows <- part1$stage == stage
  base <- part1$median_seconds[rows & part1$technique == "baseline"]
  part1$reduction_percent[rows] <- if (base > 0) 100 * (1 - part1$median_seconds[rows] / base) else 0
}

# the whole model, one technique throughout
part2 <- do.call(rbind, lapply(c("baseline", "vectorised", "cpp"), function(technique) {
  data.frame(part = "2. Whole model, one run", stage = "all", technique = technique,
             median_seconds = time_it(run_model(technique, model_settings)))
}))
part2$reduction_percent <- 100 * (1 - part2$median_seconds / part2$median_seconds[1])

#  many complete runs, sequential vs parallel
part3_settings <- expand.grid(
  model  = c("baseline", "vectorised"),
  n_runs = c(8, 100),
  stringsAsFactors = FALSE)

part3 <- do.call(rbind, lapply(seq_len(nrow(part3_settings)), function(k) {
  s  <- part3_settings[k, ]
  nb <- if (s$model == "baseline") generate_net_benefit else generate_net_benefit_vectorised
  tm <- if (s$model == "baseline") generate_transition_matrices else generate_transition_matrices_vectorised
  rows <- do.call(rbind, lapply(c(FALSE, TRUE), function(parallel) {
    data.frame(part = sprintf("3. %d complete runs, %s model", s$n_runs, s$model),
               stage = "all",
               technique = if (parallel) "parallel (multisession)" else "sequential",
               median_seconds = time_it(
                 run_model_parallel(s$n_runs, model_settings,
                                    net_benefit_fun = nb, transition_fun = tm,
                                    parallel = parallel),
                 times = 3))
  }))
  rows$reduction_percent <- 100 * (1 - rows$median_seconds / rows$median_seconds[1])
  rows
}))
#  Write out
results <- rbind(part1, part2, part3)
results$median_seconds    <- round(results$median_seconds, 4)
results$reduction_percent <- round(results$reduction_percent)
results$machine <- paste(Sys.info()[["sysname"]], R.version.string,
                         parallel::detectCores(), "cores")
results$date <- as.character(Sys.Date())

dir.create("results", showWarnings = FALSE)
write.csv(results, "results/benchmark.csv", row.names = FALSE)
print(results[, 1:5], row.names = FALSE)

