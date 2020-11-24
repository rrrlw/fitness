#####ADDITIVE LANDSCAPE#####
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
  new_FitLand(fit_arr = fit_table, params = params)
}

#####MULTIPLICATIVE LANDSCAPE#####
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
  new_FitLand(fit_arr = fit_table, params = params)
}

#####ROUGH MT. FUJI MODEL#####
# n_gene = length of string
# n_allele = number of possible letters in each string position (allow names?)
# fitness = if table of fitnesse, use as reference; if function (no params), then generate
# noise = if array, direct add; if function (no params), then generate
generate_rmf <- function(n_gene, n_allele, fitness, noise) {
  ## validate parameters
  check_n_gene(n_gene)
  check_n_allele(n_allele)
  if ("function" %in% class(fitness)) {
    check_rand_func_rmf(fitness)
  } else {
    check_fit_table(fitness, n_cols = n_gene, n_rows = n_allele)
  }
  if ("function" %in% class(noise)) {
    check_rand_func_rmf(noise)
  } else {
    dim(noise) <- c(prod(dim(noise)), 1)
    check_fit_table(noise, n_cols = 1, n_rows = n_allele ^ n_gene)
  }
  
  ## setup matrix with lattice coordinates
  fit_mat <- setup_matrix(n_dim = n_gene, n_val = n_allele)
  val_col <- ncol(fit_mat)
  
  ## generate fitness table if `fitness` is a function
  fit_table <- NULL
  if ("function" %in% class(fitness)) {
    ans <- vapply(X = seq_len(n_gene * n_allele),
                  FUN.VALUE = numeric(1),
                  FUN = function(i) fitness(),
                  USE.NAMES = FALSE)
    dim(ans) <- c(n_allele, n_gene)
    
    fit_table <- ans
  } else {
    fit_table <- fitness
  }
  
  ## fill out last column w/ appropriate fitness values
  for (curr_row in seq_len(nrow(fit_mat))) {
    # go through each column and add up fitness contribution
    fit_counter <- 0
    for (curr_col in seq_len(ncol(fit_mat) - 1)) {
      fit_counter <- fit_counter + fit_table[fit_mat[curr_row, curr_col], curr_col]
    }
    
    ## add in random component of fitness table
    rand_contr <- ifelse("function" %in% class(noise),
                         noise(),
                         noise[curr_row, 1])
    fit_mat[curr_row, ncol(fit_mat)] <- fit_counter + rand_contr
  }
  
  ## setup params for FitLand object
  fit_table <- setup_matrix_to_array(fit_mat)
  dims <- rep(n_allele, n_gene)
  type <- "rmf"
  params <- list(dims = dims,
                 type = type)
  
  ## convert to FitLand object
  new_FitLand(fit_arr = fit_table, params = params)
}

#####NK MODEL#####
generate_nk <- function(...) {
  
}

#####STICKBREAKING MODEL#####
generate_sb <- function(...) {
  
}

#####CORRELATED LANDSCAPES#####
generate_correlated <- function(...) {
  
}