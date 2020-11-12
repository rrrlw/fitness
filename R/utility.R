#####VARIABLE CLASS#####
close_to_integer <- function(x, eps = 1e-6) {
  return(is.integer(x) |
         abs(x - round(x)) < eps)
}

#####PARAMETER VALIDITY#####
# basic validity checks for rough mt fuji landscape
# error = TRUE will throw an error if invalid; error = FALSE will return TRUE/FALSE based on validity
basic_checks_rmf <- function(n_dim, n_vals, fit_table, rand_func,
                             error = TRUE) {
  # n_dim should be integer >= 1
  if (!is.null(n_dim)) {
    valid_dim(n_dim, "n_dim", error = error)
  }
  
  # n_vals should be integer >= 2
  if (!is.null(n_vals)) {
    valid_vals(n_vals, "n_vals", error = error)
  }
  
  # fit_table should be n_vals-by-n_dim numeric matrix with no missing values
  if (!is.null(fit_table)) {
    valid_fit_table_rmf(fit_table,
                        n_dim, n_vals,
                        "fit_table", error = error)
  }
  
  # rand_func should be a function with no parameters (`formals` to check the latter)
  if (!is.null(rand_func)) {
    valid_rand_func_rmf(rand_func, "rand_func", error = error)
  }
}

# valid number of dimensions in fitness landscape
# conditions:
#   - integer
#   - positive
valid_dim <- function(n_dim, param_name, error = TRUE) {
  ans <- TRUE
  
  # integer check
  if (!close_to_integer(n_dim)) {
    if (error) {
      stop(paste("Parameter", param_name, "must be an integer.",
                 "Passed value =", n_dim))
    } else {
      ans <- FALSE
    }
  # positive check
  } else if (n_dim < 1) {
    if (error) {
      stop(paste("Parameter", param_name, "must be greater than or equal to 1.",
                 "Passed value =", n_dim))
    } else {
      ans <- FALSE
    }
  }
  
  return(ans)
}

# valid alphabet length in fitness landscape
# conditions:
#   - integer
#   - greater than or equal to 2
valid_vals <- function(n_vals, param_name, error = TRUE) {
  ans <- TRUE
  
  # integer check
  if (!close_to_integer(n_vals)) {
    if (error) {
      stop(paste("Parameter", param_name, "must be an integer.",
                 "Passed value =", n_vals))
    } else {
      ans <- FALSE
    }
  # greater than or equal to 2 check
  } else if (n_vals < 2) {
    if (error) {
      stop(paste("Parameter", param_name, "must be greater than or equal to 2.",
                 "Passed value =", n_vals))
    } else {
      ans <- FALSE
    }
  }
  
  return(ans)
}

# valid fitness table for Rough Mt. Fuji fitness landscape
# conditions:
#   - matrix
#   - contains numeric values only
#   - n_vals rows
#   - n_dim cols
#   - no missing values
valid_fit_table_rmf <- function(fit_table,
                                n_dim, n_vals,
                                param_name, error = TRUE) {
  ans <- TRUE
  
  # matrix check
  if (!is.matrix(fit_table)) {
    if (error) {
      stop(paste("Parameter", param_name, "must be a matrix. Passed object has",
                 "class =", class(fit_table)))
    } else {
      ans <- FALSE
    }
  # num rows check
  } else if (nrow(fit_table) != n_vals) {
    if (error) {
      stop(paste0("Parameter ", param_name, " should have same number of rows (",
                  nrow(fit_table), ") as there are values in the fitness",
                 " landscape alphabet (", n_vals, ")."))
    } else {
      ans <- FALSE
    }
  # num cols check
  } else if (ncol(fit_table) != n_dim) {
    if (error) {
      stop(paste0("Parameter ", param_name, " should have same number of cols (",
                  ncol(fit_table), ") as there are dimensions in the fitness",
                  " landscape (", n_dim, ")."))
    } else {
      ans <- FALSE
    }
  # numeric check
  } else if (!is.numeric(fit_table)) {
    if (error) {
      stop(paste("Parameter", param_name, "must contain numeric values. Passed",
                 "object contains values with class =", class(fit_table[1, 1])))
    } else {
      ans <- FALSE
    }
  # missing values check
  } else if (sum(complete.cases(fit_table)) < nrow(fit_table)) {
    if (error) {
      stop(paste("Parameter", param_name, "should have no missing values. Passed",
                 "matrix has missing values in rows:",
                 which(!complete.cases(fit_table))))
    } else {
      ans <- FALSE
    }
  }
  
  return(ans)
}

# valid PRNG function for rough Mt. Fuji fitness landscape
# conditions:
#   - class function
#   - takes no parameters
# should also return numeric values, but can't check that w/o potentially messing up seed & reproducibility?
valid_rand_func_rmf <- function(rand_func, param_name, error = TRUE) {
  ans <- TRUE
  
  # check class function
  if (!("function" %in% class(rand_func))) {
    if (error) {
      stop(paste("Parameter", rand_func, "should be a function. Passed parameter",
                 "has class =", class(rand_func)))
    } else {
      ans <- FALSE
    }
  # check if function takes any parameters
  } else if (length(formals(rand_func)) > 0) {
    if (error) {
      stop(paste("Parameter", rand_func, "should be a function that takes zero",
                 "parameters. Passed function takes",
                 length(formals(rand_func)),
                 "parameters."))
    } else {
      ans <- FALSE
    }
  }
  
  return(ans)
}

#####FITLAND UTILITY#####
vec_neighbors <- function(vec, rmf) {
  
}

element_neighbors <- function(element, rmf) {
  
}