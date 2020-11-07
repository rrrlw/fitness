# constructor for RMF_FitLand class object
# eventually allow fit_table and rand_term to be functions?
new_RMF_FitLand <- function(n_dim, n_vals,
                            fit_table, rand_func) {
  # set variables directly
  ans_n_dim <- n_dim
  ans_n_vals <- n_vals
  ans_t_m <- NULL   # don't calculate transition matrix until needed
  
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
    fit_mat[curr_row, ncol(fit_mat)] <- fit_counter +
      rand_contr
  }
  
  # convert from matrix to array format
  fits <- array(data = NA, dim = rep(n_vals, n_dim))
  for (curr_row in seq_len(nrow(fit_mat))) {
    curr_coord <- fit_mat[curr_row, seq_len(ncol(fit_mat) - 1)]
    fits[t(curr_coord)] <- fit_mat[curr_row, ncol(fit_mat)]
  }
  
  # return RMF_FitLand object (don't calculate transition matrix until necessary)
  structure(list(n_vals = n_vals,
                 n_dim = n_dim,
                 fits = fits,
                 t_m = NULL),
            class = "RMF_FitLand")
}

# calculate transition matrix
calc_rmf_tm <- function(rmf, method = "") {
  
}

validate_RMF_FitLand <- function(rmf) {
  
}

RMF_FitLand <- function(n_dim, n_vals, fitness_table, noise_func) {
  
}

is.RMF_FitLand <- function(rmf) {
  
}