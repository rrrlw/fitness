#####PARAMETERS#####
check_n_gene <- function(n_gene) {
  if (!close_to_int(n_gene) | n_gene < 1) {
    stop(paste("n_gene parameter must be a positive integer; passed value =",
               n_gene))
  }
}

check_n_allele <- function(n_allele) {
  if (!close_to_int(n_allele) | n_allele < 1) {
    stop(paste("n_allele parameter must be a positive integer; passed value =",
               n_allele))
  }
}

check_wt_fit <- function(wt_fit) {
  if (!is.numeric(wt_fit)) {
    stop(paste("wt_fit parameter must be numeric; passed value =",
               wt_fit))
  }
}

check_max_fit <- function(max_fit, wt_fit) {
  check_wt_fit(max_fit)
  
  if (max_fit < wt_fit) {
    stop(paste("maximum fitness must be greater than wild-type fitness in",
               "stickbreaking model; passed max fitness =", max_fit, "and",
               "passed wild-type fitness =", wt_fit))
  }
}

check_mut_fit <- function(mut_fit, n_gene, n_allele) {
  if (!is.function(mut_fit) & !is.matrix(mut_fit)) {
    stop(paste("mut_fit parameter must be a matrix; passed object has class =",
               class(mut_fit)))
  }
  if (is.function(mut_fit)) {
    if (length(formals(mut_fit)) > 0) {
      stop(paste("if mut_fit parameter is a function, it must take 0 parameters;",
                 "passed function takes", length(formals(mut_fit)), "parameters"))
    }
  } else if (is.matrix(mut_fit)) {
    # make sure correct number of dimensions
    if (ncol(mut_fit) != n_gene | nrow(mut_fit) != n_allele) {
      desired <- paste(n_allele, n_gene, sep = "x")
      actual <- paste(nrow(mut_fit), ncol(mut_fit), sep = "x")
      stop(paste("mut_fit matrix must have dimensions", desired, "while passed",
                 "matrix has dimensions", actual))
    }
    
    # make sure numeric
    if (!is.numeric(mut_fit)) {
      stop(paste("mut_fit matrix must be composed of numeric values only"))
    }
    
    # make sure no missing vals
    if (!all(complete.cases(mut_fit))) {
      stop(paste("mut_fit matrix parameter cannot have missing values; passed",
                 "matrix has at least 1 missing value"))
    }
  }
}

# valid PRNG function for rough Mt. Fuji fitness landscape
# conditions:
#   - class function
#   - takes no parameters
# should also return numeric values, but can't check that w/o potentially messing up seed & reproducibility?
check_rand_func_rmf <- function(rand_func) {
  # check class function
  if (!("function" %in% class(rand_func))) {
    stop(paste("Parameter", rand_func, "should be a function. Passed parameter",
               "has class =", class(rand_func)))
  }
  
  # check if function takes any parameters
  if (length(formals(rand_func)) > 0) {
    stop(paste("Parameter", rand_func, "should be a function that takes zero",
               "parameters. Passed function takes",
               length(formals(rand_func)),
               "parameters."))
  }
}

# valid fitness table for landscape
# conditions:
#   - matrix
#   - contains numeric values only
#   - n_vals rows
#   - n_dim cols
#   - no missing values
check_fit_table <- function(fit_table,
                            n_cols, n_rows) {
  # matrix check
  if (!is.matrix(fit_table)) {
    stop(paste("Parameter", param_name, "must be a matrix. Passed object has",
               "class =", class(fit_table)))
  # num rows check
  } else if (nrow(fit_table) != n_vals) {
    stop(paste0("Parameter ", param_name, " should have same number of rows (",
                nrow(fit_table), ") as there are values in the fitness",
                " landscape alphabet (", n_vals, ")."))
  # num cols check
  } else if (ncol(fit_table) != n_dim) {
    stop(paste0("Parameter ", param_name, " should have same number of cols (",
                ncol(fit_table), ") as there are dimensions in the fitness",
                " landscape (", n_dim, ")."))
  # numeric check
  } else if (!is.numeric(fit_table)) {
    stop(paste("Parameter", param_name, "must contain numeric values. Passed",
               "object contains values with class =", class(fit_table[1, 1])))
  # missing values check
  } else if (sum(complete.cases(fit_table)) < nrow(fit_table)) {
    stop(paste("Parameter", param_name, "should have no missing values. Passed",
               "matrix has missing values in rows:",
               which(!complete.cases(fit_table))))
  }
}
