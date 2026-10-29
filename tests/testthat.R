# Runs the test suite: devtools::test() during development, or R CMD check
library(testthat)
library(HIPS)
test_check("HIPS")
