# function that calculates lattice distance between two genotypes/locations
#   in fitness landscapes
lattice_dist <- function(coord1, coord2) {

}

double_equals <- function(x, y, eps = 1e-6) {
  return(abs(x - y) < eps)
}

close_to_int <- function(x, eps = 1e-6) {
  return(double_equals(x, round(x), eps = eps))
}

# n is an integer
# modified base conversion (away from 10) for 1-based array indexing
single_to_multi <- function(n, n_dim, n_vals) {
  ans <- integer(n_dim)
  counter <- length(ans)
  n <- n - 1
  while (n > 0 & counter > 0) {
    ans[counter] <- n %% n_vals
    n <- as.integer(n / n_vals)
    counter <- counter - 1
  }
  return(rev(ans + 1))
}

# n is an integer vector
# modified base conversion (back to 10) for 1-based array indexing
# n_dim must equal length(n)
multi_to_single <- function(n, n_vals) {
  ans <- 0
  counter <- length(n)
  curr_val <- 1
  n <- rev(n - 1)
  while (counter > 0) {
    ans <- ans + n[counter] * curr_val
    counter <- counter - 1
    curr_val <- curr_val * n_vals
  }
  return(ans + 1)
}

#####SETUP MATRIX FUNCTIONS#####
# generate setup matrix for the generate_* functions
# last column is empty w/ NAs
setup_matrix <- function(n_dim, n_val) {
  # setup matrix (long format)
  fit_mat <- matrix(NA, ncol = n_dim + 1, nrow = n_val ^ n_dim)
  
  # setup first n_dim columns w/ appropriate lattice coordinates
  curr_period <- n_val ^ (n_dim - 1)
  curr_rep <- 1
  for (i in seq_len(ncol(fit_mat) - 1)) {
    fit_mat[, i] <- rep(seq_len(n_val),
                        each = curr_period,
                        times = curr_rep)
    
    curr_period <- curr_period / n_val
    curr_rep <- curr_rep * n_val
  }
  
  return(fit_mat)
}

setup_matrix_to_array <- function(setup_mat) {
  # extract basic parameters
  n_dim <- ncol(setup_mat) - 1
  n_val <- length(unique(setup_mat[, 1]))
  
  # convert to array
  ans <- array(data = NA, dim = rep(n_val, n_dim))
  for (curr_row in seq_len(n_val ^ n_dim)) {
    curr_coord <- setup_mat[curr_row, seq_len(n_dim)]
    ans[t(curr_coord)] <- setup_mat[curr_row, ncol(setup_mat)]
  }
  
  # return array
  return(ans)
}
