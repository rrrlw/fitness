# n_dim = number of genes
# n_val = allele values per gene
# wt_fit = wild-type fitness value (assumed to be equal to rep(1, n_dim))
# mut_fit = if numeric matrix, generate landscape based on wt + appropriate muts
#           if function (no params), generate numeric vector, then do same
generate_add <- function(n_gene, n_allele, wt_fit, mut_fit) {
  ## validate parameters
  check_n_gene(n_gene)
  check_n_allele(n_allele)
  check_wt_fit(wt_fit)
  check_mut_fit(mut_fit)
  
  ## setup matrix with lattice coordinates
  temp_mat <- setup_matrix(n_dim = n_gene, n_val = n_allele)
  val_col <- ncol(temp_mat)
  
  ## fill in last column
  # get mutations as matrix if not in that form
  mut_matrix <- NULL
  if (is.matrix(mut_fit)) {
    mut_matrix <- mut_fit
  } else { # assume it's a function
    mut_matrix <- matrix(NA, nrow = n_allele, ncol = n_gene)
    mut_matrix[1, ] <- 0 # wild-type
    for (curr_gene in seq_len(ncol(mut_matrix))) {
      for (curr_allele in seq(from = 2, to = nrow(mut_matrix), by = 1)) {
        # not wild-type
        mut_matrix[curr_allele, curr_gene] <- mut_fit()
      }
    }
  }
  
  ## setup fitness table based on mut_matrix
  for (curr_genotype in seq_len(nrow(temp_mat))) {
    # start at wild-type fitness
    fit <- wt_fit
    
    # keep adjusting based on mutations
    for (curr_gene in seq(from = 1, to = n_gene, by = 1)) {
      fit <- fit + mut_matrix[temp_mat[curr_genotype, curr_gene], curr_gene]
    }
    
    # set value in setup matrix
    temp_mat[curr_genotype, val_col] <- fit
  }
  
  ## setup params for FitLand object
  fit_table <- setup_matrix_to_array(temp_mat)
  dims <- rep(n_allele, n_gene)
  type <- "add"
  params <- list(dims = dims,
                 type = type)
  
  ## convert to FitLand object
  new_FitLand(fit_table = fit_table, params = params)
}

# same parameters as generate_add
generate_mult <- function(n_gene, n_allele, wt_fit, mut_fit) {
  ## validate parameters
  check_n_gene(n_gene)
  check_n_allele(n_allele)
  check_wt_fit(wt_fit)
  check_mut_fit(mut_fit)
  
  ## setup matrix with lattice coordinates
  temp_mat <- setup_matrix(n_dim = n_gene, n_val = n_allele)
  val_col <- ncol(temp_mat)
  
  ## fill in last column
  # get mutations as matrix if not in that form
  mut_matrix <- NULL
  if (is.matrix(mut_fit)) {
    mut_matrix <- mut_fit
  } else { # assume it's a function
    mut_matrix <- matrix(NA, nrow = n_allele, ncol = n_gene)
    mut_matrix[1, ] <- 0 # wild-type
    for (curr_gene in seq_len(ncol(mut_matrix))) {
      for (curr_allele in seq(from = 2, to = nrow(mut_matrix), by = 1)) {
        # not wild-type
        mut_matrix[curr_allele, curr_gene] <- mut_fit()
      }
    }
  }
  
  ## setup fitness table based on mut_matrix
  for (curr_genotype in seq_len(nrow(temp_mat))) {
    # start at wild-type fitness
    fit <- wt_fit
    
    # keep adjusting based on mutations
    for (curr_gene in seq(from = 1, to = n_gene, by = 1)) {
      fit <- fit * (1 + mut_matrix[temp_mat[curr_genotype, curr_gene], curr_gene])
    }
    
    # set value in setup matrix
    temp_mat[curr_genotype, val_col] <- fit
  }
  
  ## setup params for FitLand object
  fit_table <- setup_matrix_to_array(temp_mat)
  dims <- rep(n_allele, n_gene)
  type <- "mult"
  params <- list(dims = dims,
                 type = type)
  
  ## convert to FitLand object
  new_FitLand(fit_table = fit_table, params = params)
}

# n_dim = length of string
# n_val = number of possible letters in each string position (allow names?)
# fitness = if array, direct copy; if function (0 or 1 loc param), then generate
# noise = if array, direct add; if function (no params), then generate
generate_rmf <- function(n_dim, n_val, fitness, noise) {
  
}

generate_nk <- function(...) {
  
}

generate_sb <- function(...) {
  
}

generate_correlated <- function(...) {
  
}