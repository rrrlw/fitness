# n_dim = number of genes
# n_val = allele values per gene
# wt_fit = wild-type fitness value (assumed to be equal to rep(1, n_dim))
# mut_fit = if numeric matrix, generate landscape based on wt + appropriate muts
#           if function (no params), generate numeric vector, then do same
generate_add <- function(n_gene, n_allele, wt_fit, mut_fit) {
  ## validate parameters
  if (!close_to_int(n_gene) | n_gene < 1) {
    stop(paste("n_gene parameter must be a positive integer; passed value =",
               n_gene))
  }
  if (!close_to_int(n_allele) | n_allele < 1) {
    stop(paste("n_allele parameter must be a positive integer; passed value =",
               n_allele))
  }
  if (!is.numeric(wt_fit)) {
    stop(paste("wt_fit parameter must be numeric; passed value =",
               wt_fit))
  }
  if (!is.function(mut_fit) & !is.matrix(mut_fit)) {
    stop(paste("mut_fit parameter must be a matrix; passed object has class =",
               class(mut_fit)))
  }
  if (is.function(mut_fit) & length(formals(mut_fit)) > 0) {
    stop(paste("if mut_fit parameter is a function, it must take 0 parameters;",
               "passed function takes", length(formals(mut_fit)), "parameters"))
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
  
  # setup fitness table based on mut_matrix
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
  
  # setup params for FitLand object
  fit_table <- setup_matrix_to_array(temp_mat)
  dims <- rep(n_allele, n_gene)
  type <- "add"
  params = list(dims = dims,
                type = type)
  
  # convert to FitLand object
  new_FitLand(fit_table = fit_table, params = params)
}

generate_mult <- function(...) {
  
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