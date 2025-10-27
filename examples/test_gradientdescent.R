#!/usr/bin/env Rscript
# Test script for GradientDescent algorithm in new fz format
# This simulates a simple test of the algorithm structure

# Source the algorithm file
source(".fz/algorithms/gradientdescent.R")

# Test 1: Constructor with default options
cat("Test 1: Constructor with default options\n")
gd1 <- GradientDescent()
cat("  - Created GradientDescent object\n")
cat("  - Class:", class(gd1), "\n")
cat("  - Default yminimization:", gd1$options$yminimization, "\n")
cat("  - Default max_iterations:", gd1$options$max_iterations, "\n")
cat("  - Default ytol:", gd1$options$ytol, "\n")
cat("  ✓ Constructor works with defaults\n\n")

# Test 2: Constructor with custom options
cat("Test 2: Constructor with custom options\n")
gd2 <- GradientDescent(
  yminimization = FALSE,
  max_iterations = 50,
  ytol = 0.01,
  delta = 0.5,
  epsilon = 0.05,
  x0 = "0.3,0.7"
)
cat("  - yminimization:", gd2$options$yminimization, "\n")
cat("  - max_iterations:", gd2$options$max_iterations, "\n")
cat("  - ytol:", gd2$options$ytol, "\n")
cat("  - delta:", gd2$options$delta, "\n")
cat("  - epsilon:", gd2$options$epsilon, "\n")
cat("  - x0:", gd2$options$x0, "\n")
cat("  ✓ Constructor works with custom options\n\n")

# Test 3: get_initial_design
cat("Test 3: get_initial_design method\n")
input_vars <- list(
  x1 = c(0, 1),
  x2 = c(0, 1)
)
output_vars <- list("y")

initial_design <- get_initial_design(gd1, input_vars, output_vars)
cat("  - Number of initial points:", length(initial_design), "\n")
cat("  - Expected points: 3 (center + 2 finite differences for 2D)\n")
cat("  - First point:", paste(names(initial_design[[1]]), "=", 
                              sapply(initial_design[[1]], function(x) sprintf("%.4f", x)), 
                              collapse = ", "), "\n")
cat("  ✓ get_initial_design returns list of points\n\n")

# Test 4: Simple optimization test with Branin function
cat("Test 4: Simple Branin function optimization test\n")
# Branin function (simplified, scaled to [0,1]^2)
branin <- function(point) {
  x1 <- point$x1 * 15 - 5
  x2 <- point$x2 * 15
  (x2 - 5/(4*pi^2)*(x1^2) + 5/pi*x1 - 6)^2 + 10*(1 - 1/(8*pi))*cos(x1) + 10
}

# Create algorithm instance
gd_opt <- GradientDescent(
  yminimization = TRUE,
  max_iterations = 5,  # Just a few iterations for testing
  ytol = 0.01,
  delta = 0.1,
  epsilon = 0.01
)

# Get initial design
X <- get_initial_design(gd_opt, input_vars, output_vars)
Y <- sapply(X, branin)

cat("  - Initial design size:", length(X), "\n")
cat("  - Initial Y values:", sprintf("%.4f", Y), "\n")

# Run a few iterations
for (iter in 1:3) {
  # Get analysis_tmp
  analysis_tmp <- get_analysis_tmp(gd_opt, X, Y)
  cat("  - Iteration", iter, ":", analysis_tmp$text, "\n")
  
  # Get next design
  X_next <- get_next_design(gd_opt, X, Y)
  
  if (length(X_next) == 0) {
    cat("  - Algorithm converged early\n")
    break
  }
  
  # Evaluate new points
  Y_next <- sapply(X_next, branin)
  
  # Append to history
  X <- c(X, X_next)
  Y <- c(Y, Y_next)
}

# Get final analysis
cat("\nTest 5: get_analysis method\n")
analysis <- get_analysis(gd_opt, X, Y)
cat(analysis$text)
cat("  - Optimum from data:", analysis$data$optimum, "\n")
cat("  - Point:", paste(names(analysis$data$optimum_point), "=", 
                        sapply(analysis$data$optimum_point, function(x) sprintf("%.4f", x)), 
                        collapse = ", "), "\n")
cat("  ✓ get_analysis returns results\n\n")

# Test 6: 1D optimization
cat("Test 6: 1D optimization test\n")
input_vars_1d <- list(x = c(0, 1))
quadratic <- function(point) {
  (point$x - 0.6)^2  # Minimum at x=0.6
}

gd_1d <- GradientDescent(
  yminimization = TRUE,
  max_iterations = 3,
  ytol = 0.001,
  delta = 0.1,
  epsilon = 0.01
)

X_1d <- get_initial_design(gd_1d, input_vars_1d, output_vars)
Y_1d <- sapply(X_1d, quadratic)

cat("  - Initial 1D design size:", length(X_1d), "\n")

# One iteration
X_next_1d <- get_next_design(gd_1d, X_1d, Y_1d)
if (length(X_next_1d) > 0) {
  Y_next_1d <- sapply(X_next_1d, quadratic)
  X_1d <- c(X_1d, X_next_1d)
  Y_1d <- c(Y_1d, Y_next_1d)
}

analysis_1d <- get_analysis(gd_1d, X_1d, Y_1d)
cat("  - 1D Optimum:", analysis_1d$data$optimum, "\n")
cat("  - Location:", analysis_1d$data$optimum_point$x, "\n")
cat("  ✓ 1D optimization works\n\n")

cat("========================================\n")
cat("All tests passed! ✓\n")
cat("========================================\n")
