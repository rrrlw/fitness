# constructor for RMF_FitLand class object
# eventually allow fit_table and rand_term to be functions?
# transition matrix will be handled separately (backend of evolve() or equivalent function - S3 generic?)
new_RMF_FitLand <- function(n_dim, n_vals,
                            fit_table,
                            rand_func = runif(1)) {
  # conduct basic checks to check if RMF parameters are valid
  basic_checks_rmf(n_dim, n_vals, fit_table, rand_func,
                   error = TRUE)
  
  # set variables directly
  ans_n_dim <- n_dim
  ans_n_vals <- n_vals
  
  ## calculate deterministic component of fitness table
  # setup matrix (long format)
  fit_mat <- matrix(NA, ncol = n_dim + 1, nrow = n_vals ^ n_dim)
  
  # setup first n_dim columns w/ appropriate lattice coordinates
  curr_period <- n_vals ^ (n_dim - 1)
  curr_rep <- 1
  for (i in seq_len(ncol(fit_mat) - 1)) {
    fit_mat[, i] <- rep(seq_len(n_vals),
                        each = curr_period,
                        times = curr_rep)
    
    curr_period <- curr_period / n_vals
    curr_rep <- curr_rep * n_vals
  }
  
  # fill out last column w/ appropriate fitness values
  for (curr_row in seq_len(nrow(fit_mat))) {
    
    # go through each column and add up fitness contribution
    fit_counter <- 0
    for (curr_col in seq_len(ncol(fit_mat) - 1)) {
      fit_counter <- fit_counter + fit_table[fit_mat[curr_row, curr_col], curr_col]
    }
    
    # add in random component of fitness table
    rand_contr <- rand_func()
    fit_mat[curr_row, ncol(fit_mat)] <- fit_counter + rand_contr
  }
  
  # convert from matrix to array format
  fits <- array(data = NA, dim = rep(n_vals, n_dim))
  for (curr_row in seq_len(nrow(fit_mat))) {
    curr_coord <- fit_mat[curr_row, seq_len(ncol(fit_mat) - 1)]
    fits[t(curr_coord)] <- fit_mat[curr_row, ncol(fit_mat)]
  }
  
  # return RMF_FitLand object (don't calculate transition matrix until necessary)
  structure(list(n_dim = n_dim,
                 n_vals = n_vals,
                 fits = fits),
            class = "RMF_FitLand")
}

# if error = TRUE, throw error if invalid, silent if valid
# if error = FALSE, return FALSE if invalid, return TRUE if valid
validate_RMF_FitLand <- function(rmf, error = TRUE) {
  # basic checks
  basic_checks_rmf(n_dim = rmf$n_dim, n_vals = rmf$n_vals,
                   fit_table = NULL,
                   rand_func = NULL,
                   error = error)
  
  # add complex checks
}

RMF_FitLand <- function(n_dim, n_vals, fitness_table, noise_func) {
  
}

# TRUE if it can pass validation tests; FALSE otherwise
is.RMF_FitLand <- function(rmf) {
  return(validate_RMF_FitLand(rmf, error = FALSE))
}
