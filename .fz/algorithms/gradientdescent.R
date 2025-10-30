
#title: First-order local optimization algorithm
#author: Yann Richet
#help: http://en.wikipedia.org/wiki/Gradient_descent
#tags: optimization
#type: optimization
#options: yminimization=true;ytol=0.1;max_iterations=100;delta=1;epsilon=0.01;target=Inf;x0=
#options.help: yminimization='Minimize output value ?';ytol='Convergence precision on output value';max_iterations='Maximum number of iterations';delta='Gradient step factor (initial value)';epsilon='Relative finite difference step';target='Output target limit (early convergence)';x0='Starting input values (comma separated)'

# Constructor for GradientDescent S3 class
GradientDescent <- function(...) {
  # Get options from ... arguments
  opts <- list(...)
  
  # Create object with initial state
  # Use an environment for mutable state (idiomatic S3 pattern)
  state <- new.env(parent = emptyenv())
  state$i <- 0
  state$input <- NULL
  state$x0 <- NULL
  state$delta <- NULL
  
  obj <- list(
    options = list(
      yminimization = isTRUE(as.logical(
        ifelse(is.null(opts$yminimization), TRUE, opts$yminimization)
      )),
      ytol = as.numeric(
        ifelse(is.null(opts$ytol), 0.1, opts$ytol)
      ),
      max_iterations = as.integer(
        ifelse(is.null(opts$max_iterations), 100, opts$max_iterations)
      ),
      delta = as.numeric(
        ifelse(is.null(opts$delta), 1, opts$delta)
      ),
      epsilon = as.numeric(
        ifelse(is.null(opts$epsilon), 0.01, opts$epsilon)
      ),
      target = as.numeric(
        ifelse(is.null(opts$target), Inf, opts$target)
      ),
      x0 = if (is.null(opts$x0) || opts$x0 == '') NULL else as.numeric(strsplit(as.character(opts$x0), ",")[[1]])
    ),
    state = state  # Environment for mutable state
  )
  
  # Adjust target based on minimization/maximization
  if (!obj$options$yminimization) {
    if (isTRUE(obj$options$target == -Inf)) {
      obj$options$target <- Inf
    }
  }
  if (obj$options$yminimization) {
    if (isTRUE(obj$options$target == Inf)) {
      obj$options$target <- -Inf
    }
  }
  
  # Initialize mutable delta in state
  state$delta <- obj$options$delta
  
  # Set S3 class
  class(obj) <- "GradientDescent"
  
  return(obj)
}

# Generic function definitions (if not already defined)
if (!exists("get_initial_design")) {
  get_initial_design <- function(obj, ...) UseMethod("get_initial_design")
}

if (!exists("get_next_design")) {
  get_next_design <- function(obj, ...) UseMethod("get_next_design")
}

if (!exists("get_analysis")) {
  get_analysis <- function(obj, ...) UseMethod("get_analysis")
}

if (!exists("get_analysis_tmp")) {
  get_analysis_tmp <- function(obj, ...) UseMethod("get_analysis_tmp")
}

# Method: get_initial_design
get_initial_design.GradientDescent <- function(obj, input_variables, output_variables) {
  # Store input variables in mutable state
  # input_variables is a named list: list(var1 = c(min, max), var2 = c(min, max))
  obj$state$input <- input_variables
  obj$state$i <- 0
  
  d <- length(input_variables)
  
  # Generate or use provided initial point x0
  if (!is.null(obj$options$x0)) {
    x0 <- rep(obj$options$x0, d)
    if (length(x0) > d) {
      x0 <- x0[1:d]
    }
    x0 <- matrix(x0, ncol = d)
    colnames(x0) <- names(input_variables)
    x0 <- to01(x0, obj$state$input)
    colnames(x0) <- names(input_variables)
  } else {
    x0 <- matrix(runif(d), ncol = d)
    colnames(x0) <- names(input_variables)
  }
  
  if (ncol(x0) > 1) {
    x0 <- x0[1, ]
  }
  obj$state$x0 <- x0
  
  # Ask for finite differences around x0
  x <- askfinitedifferences(x0, obj$options$epsilon)
  colnames(x) <- names(input_variables)
  
  # Convert from [0,1] to actual bounds and return as list of points
  x_real <- from01(x, obj$state$input)
  
  # Convert matrix to list of points, preserving names
  points <- list()
  for (i in 1:nrow(x_real)) {
    point <- as.list(x_real[i, ])
    names(point) <- names(obj$state$input)
    points[[i]] <- point
  }
  
  return(points)
}

# Method: get_next_design
get_next_design.GradientDescent <- function(obj, X, Y) {
  # Check max iterations
  if (obj$state$i >= obj$options$max_iterations) {
    return(list())  # Empty list signals finished
  }
  
  # Filter out NULL/NA values
  Y_valid <- Y[!sapply(Y, is.null) & !is.na(Y)]
  Y_valid <- unlist(Y_valid)
  
  # Check target convergence
  if (obj$options$yminimization) {
    if (length(Y_valid) > 0 && min(Y_valid) < obj$options$target) {
      return(list())
    }
  } else {
    if (length(Y_valid) > 0 && max(Y_valid) > obj$options$target) {
      return(list())
    }
  }
  
  # Convert X from list of points to matrix
  X_mat <- points_to_matrix(X, names(obj$state$input))
  
  # Convert to [0,1] space
  X_01 <- to01(X_mat, obj$state$input)
  
  d <- ncol(X_01)
  n <- nrow(X_01)
  
  # Get previous d+1 points
  prevXn <- X_01[(n - d):n, , drop = FALSE]
  prevYn <- Y_valid[(n - d):n]
  
  # Check convergence based on ytol
  if (obj$state$i > 0) {
    if (abs(Y_valid[n - d] - Y_valid[n - d - 1 - d]) < obj$options$ytol) {
      return(list())
    }
    
    # Adjust delta if not improving
    if (obj$options$yminimization) {
      if (Y_valid[n - d] >= Y_valid[n - d - 1 - d]) {
        obj$state$delta <- obj$state$delta / 2
        prevXn <- X_01[(n - d - d - 1):(n - d - 1), , drop = FALSE]
        prevYn <- Y_valid[(n - d - d - 1):(n - d - 1)]
      }
    }
    if (!obj$options$yminimization) {
      if (Y_valid[n - d] <= Y_valid[n - d - 1 - d]) {
        obj$state$delta <- obj$state$delta / 2
        prevXn <- X_01[(n - d - d - 1):(n - d - 1), , drop = FALSE]
        prevYn <- Y_valid[(n - d - d - 1):(n - d - 1)]
      }
    }
  }
  
  if (d == 1) {
    prevXn <- matrix(prevXn, ncol = 1)
  }
  
  # Calculate gradient (normalized)
  grad_norm <- gradient(prevXn, prevYn) / (max(Y_valid) - min(Y_valid))
  
  # Adjust delta if gradient is too large
  if (max(abs(grad_norm)) * obj$state$delta > 1) {
    obj$state$delta <- obj$state$delta / max(abs(grad_norm))
  }
  
  # Calculate next point
  if (obj$options$yminimization) {
    xnext <- prevXn[1, ] - grad_norm * obj$state$delta
  } else {
    xnext <- prevXn[1, ] + grad_norm * obj$state$delta
  }
  xnext <- t(xnext)
  
  # Keep xnext in [0,1] bounds with reflection
  for (t in 1:d) {
    while ((xnext[t] > 1.0) | (xnext[t] < 0)) {
      if (xnext[t] > 1.0) {
        xnext[t] <- 2.0 - xnext[t]
      }
      if (xnext[t] < 0.0) {
        xnext[t] <- 0.0 - xnext[t]
      }
    }
  }
  
  obj$state$i <- obj$state$i + 1
  
  # Ask for finite differences around xnext
  x <- askfinitedifferences(xnext, obj$options$epsilon)
  colnames(x) <- names(obj$state$input)
  
  # Convert from [0,1] to actual bounds
  x_real <- from01(x, obj$state$input)
  
  # Convert matrix to list of points, preserving names
  points <- list()
  for (i in 1:nrow(x_real)) {
    point <- as.list(x_real[i, ])
    names(point) <- names(obj$state$input)
    points[[i]] <- point
  }
  
  return(points)
}

# Method: get_analysis
get_analysis.GradientDescent <- function(obj, X, Y) {
  analysis_dict <- list(text = "", data = list())
  
  # Filter out NULL/NA values
  Y_valid <- Y[!sapply(Y, is.null) & !is.na(Y)]
  Y_valid <- unlist(Y_valid)
  
  if (length(Y_valid) < 1) {
    analysis_dict$text <- "No valid results to analyze"
    analysis_dict$data <- list(valid_samples = 0)
    return(analysis_dict)
  }
  
  # Convert X from list of points to matrix
  X_mat <- points_to_matrix(X, names(obj$state$input))
  
  # Find optimum
  if (isTRUE(obj$options$yminimization)) {
    m <- min(Y_valid)
    m.ix <- which.min(Y_valid)
  } else {
    m <- max(Y_valid)
    m.ix <- which.max(Y_valid)
  }
  m.ix <- m.ix[1]
  x_opt <- X_mat[m.ix, ]
  
  d <- ncol(X_mat)
  
  # Store data
  analysis_dict$data <- list(
    optimum = m,
    optimum_point = as.list(x_opt),
    n_evaluations = length(Y_valid),
    iterations = obj$state$i
  )
  
  # Create text summary
  opt_type <- if (obj$options$yminimization) "minimum" else "maximum"
  analysis_dict$text <- sprintf(
    "Gradient Descent Results:
  Iterations: %d
  Total evaluations: %d
  %s: %.6f
  Found at: %s
",
    obj$state$i,
    length(Y_valid),
    opt_type,
    m,
    paste(paste(names(x_opt), "=", sprintf("%.6f", x_opt), sep = ""), collapse = "; ")
  )
  
  # Try to create HTML with plot
  tryCatch({
    # Calculate color based on distance from optimum
    if (obj$options$yminimization) {
      red <- (as.matrix(Y_valid) - min(Y_valid)) / (max(Y_valid) - min(Y_valid))
    } else {
      red <- (max(Y_valid) - as.matrix(Y_valid)) / (max(Y_valid) - min(Y_valid))
    }
    
    # Create plot
    png_file <- tempfile(fileext = ".png")
    png(png_file, width = 600, height = 600, bg = "transparent")
    
    if (d > 1) {
      pairs(cbind(X_mat, Y_valid), 
            col = rgb(r = red, g = 0, b = 1 - red),
            labels = c(names(X_mat), "Y"))
    } else {
      plot(x = X_mat[, 1], y = Y_valid,
           xlab = names(X_mat),
           ylab = "Y",
           col = rgb(r = red, g = 0, b = 1 - red),
           pch = 20)
    }
    
    dev.off()
    
    # Convert to base64
    if (requireNamespace("base64enc", quietly = TRUE)) {
      img_base64 <- base64enc::base64encode(png_file)
      
      html_output <- sprintf(
        '<div>
  <h3>Gradient Descent Results</h3>
  <p><strong>%s:</strong> %.6f</p>
  <p><strong>Found at:</strong> %s</p>
  <p><strong>Iterations:</strong> %d</p>
  <p><strong>Total evaluations:</strong> %d</p>
  <img src="data:image/png;base64,%s" alt="Optimization plot" style="max-width:600px;"/>
</div>',
        opt_type,
        m,
        paste(paste(names(x_opt), "=", sprintf("%.6f", x_opt), sep = ""), collapse = "; "),
        obj$state$i,
        length(Y_valid),
        img_base64
      )
      analysis_dict$html <- html_output
    }
    
    # Clean up temp file
    unlink(png_file)
  }, error = function(e) {
    # If plotting fails, just skip it
  })
  
  return(analysis_dict)
}

# Method: get_analysis_tmp
get_analysis_tmp.GradientDescent <- function(obj, X, Y) {
  # Filter out NULL/NA values
  Y_valid <- Y[!sapply(Y, is.null) & !is.na(Y)]
  Y_valid <- unlist(Y_valid)
  
  if (length(Y_valid) < 1) {
    return(list(
      text = sprintf("  Progress: Iteration %d, no valid samples yet", obj$state$i),
      data = list(iteration = obj$state$i, valid_samples = 0)
    ))
  }
  
  # Find current best
  if (isTRUE(obj$options$yminimization)) {
    current_best <- min(Y_valid)
  } else {
    current_best <- max(Y_valid)
  }
  
  opt_type <- if (obj$options$yminimization) "min" else "max"
  
  return(list(
    text = sprintf(
      "  Progress: Iteration %d, %d evaluations, current %s=%.6f",
      obj$state$i,
      length(Y_valid),
      opt_type,
      current_best
    ),
    data = list(
      iteration = obj$state$i,
      n_evaluations = length(Y_valid),
      current_best = current_best
    )
  ))
}

# Helper function: points_to_matrix (not a method, internal use only)
points_to_matrix <- function(X, input_names) {
  X_list <- lapply(X, function(point) {
    unlist(point[input_names])
  })
  X_mat <- do.call(rbind, X_list)
  if (is.null(dim(X_mat))) {
    # Single row case
    X_mat <- matrix(X_mat, nrow = 1)
  }
  colnames(X_mat) <- input_names
  return(X_mat)
}

# Helper function: askfinitedifferences (not a method, internal use only)
askfinitedifferences <- function(x, epsilon) {
  xd <- matrix(x, nrow = 1)
  for (i in 1:length(x)) {
    xdi <- as.array(x)
    if (xdi[i] + epsilon > 1.0) {
      xdi[i] <- xdi[i] - epsilon
    } else {
      xdi[i] <- xdi[i] + epsilon
    }
    xd <- rbind(xd, matrix(xdi, nrow = 1))
  }
  return(xd)
}

# Helper function: gradient (not a method, internal use only)
gradient <- function(xd, yd) {
  d <- ncol(xd)
  grad <- rep(0, d)
  for (i in 1:d) {
    grad[i] <- (yd[i + 1] - yd[1]) / (xd[i + 1, i] - xd[1, i])
  }
  return(grad)
}

# Helper function: from01 (not a method, internal use only)
from01 <- function(X, inp) {
  nX <- colnames(X)
  if (is.null(nX)) nX <- names(X)
  for (i in 1:ncol(X)) {
    namei <- nX[i]
    bounds <- inp[[namei]]
    # Validate bounds structure
    if (!is.numeric(bounds) || length(bounds) != 2) {
      stop(paste("Invalid bounds for variable", namei, ": expected c(min, max)"))
    }
    min_val <- bounds[1]
    max_val <- bounds[2]
    X[, i] <- X[, i] * (max_val - min_val) + min_val
  }
  return(X)
}

# Helper function: to01 (not a method, internal use only)
to01 <- function(X, inp) {
  nX <- colnames(X)
  if (is.null(nX)) nX <- names(X)
  for (i in 1:ncol(X)) {
    namei <- nX[i]
    bounds <- inp[[namei]]
    # Validate bounds structure
    if (!is.numeric(bounds) || length(bounds) != 2) {
      stop(paste("Invalid bounds for variable", namei, ": expected c(min, max)"))
    }
    min_val <- bounds[1]
    max_val <- bounds[2]
    X[, i] <- (X[, i] - min_val) / (max_val - min_val)
  }
  return(X)
}
