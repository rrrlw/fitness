#####ADDITIVE LANDSCAPE#####
# n_dim = number of genes
# n_val = allele values per gene
# wt_fit = wild-type fitness value (assumed to be equal to rep(1, n_dim))
# mut_fit = if numeric matrix, generate landscape based on wt + appropriate muts
#           if function (no params), generate numeric vector, then do same
#' Additive Model for Fitness Landscapes
#' 
#' Generates a fitness landscape using the additive model, as characterized in
#' Nagel et al. (2012) <doi:10.1534.genetics.111.132134>.
#' 
#' @param n_gene number of genes
#' @param n_allele number of alleles per gene
#' @param wt_fit wild-type genotype fitness
#' @param mut_fit if of class \code{matrix}, contains numeric fitness values
#'   for each allele and gene combination (each column is a gene, each row is
#'   an allele); if a \code{function}, must have zero parameters and return a
#'   single \code{numeric} value upon each call (it will be called repeatedly
#'   when generating the landscape)
#' @return fitness landscape stored in \code{FitLand} class object
#' @export
#' @examples
#' num_genes <- 4
#' num_alleles <- 3
#' landscape <- generate_add(num_genes,
#'                           num_alleles,
#'                           mut_fit = function() {stats::rnorm(1)})
#'                           
#' print(landscape)
generate_add <- function(n_gene, n_allele, wt_fit = 0, mut_fit) {
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
  params <- list(n_gene = n_gene,
                 n_allele = n_allele,
                 wt = wt_fit,
                 mut = mut_matrix,
                 type = "add")
  
  ## convert to FitLand object
  new_FitLand(fit_arr = fit_table, params = params)
}

#####MULTIPLICATIVE LANDSCAPE#####
# same parameters as generate_add
#' Multiplicative Model for Fitness Landscapes
#' 
#' Generates a fitness landscape using the multiplicative model, as
#' characterized in Nagel et al. (2012) <doi:10.1534.genetics.111.132134>.
#' 
#' @inheritParams generate_add
#' @return fitness landscape stored in \code{FitLand} class object
#' @export
#' @examples
#' num_genes <- 4
#' num_alelles <- 3
#' landscape <- generate_mult(num_genes,
#'                            num_alleles,
#'                            mut_fit = function() {stats::runif(1, max = 2)})
#'                           
#' print(landscape)
generate_mult <- function(n_gene, n_allele, wt_fit = 1, mut_fit) {
  ## validate parameters
  check_n_gene(n_gene)
  check_n_allele(n_allele)
  check_wt_fit(wt_fit)
  check_mut_fit(mut_fit, n_gene, n_allele)
  
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
  params <- list(n_gene = n_gene,
                 n_allele = n_allele,
                 wt = wt_fit,
                 mut = mut_matrix,
                 type = "mult")
  
  ## convert to FitLand object
  new_FitLand(fit_arr = fit_table, params = params)
}

#####ROUGH MT. FUJI MODEL#####
# n_gene = length of string
# n_allele = number of possible letters in each string position (allow names?)
# fitness = if table of fitnesse, use as reference; if function (no params), then generate
# noise = if array, direct add; if function (no params), then generate
#' Rough Mount Fuji Model for Fitness Landscapes
#' 
#' Generates a fitness landscape using the rough Mt. Fuji model, as
#' characterized in Aita et al. (2000)
#' <doi:10.1002/(SICI)1097-0282(200007)54:1<64::AID-BIP70>3.0.CO;2-R>.
#' 
#' @inheritParams generate_add
#' @param fitness if of class \code{matrix}, contains the deterministic fitness
#'   of each allele-gene combination (each column is a gene, each row is an
#'   allele); if of class \code{function}, must require zero parameters,
#'   return a single \code{numeric} value, and will be repeatedly called for
#'   each genotype to determine deterministic fitness
#' @param noise if of class \code{array}, contains the (pseudo)random
#'   contribution to each genotype; if of class \code{function}, must require
#'   zero parameters, return a single \code{numeric} value, and will be
#'   repeatedly called for each genotype to determine (pseudo)random fitness
#'   contribution
#' @return fitness landscape stored in \code{FitLand} class object
#' @export
#' @examples
#' num_genes <- 4
#' num_alelles <- 3
#' landscape <- generate_rmf(num_genes,
#'                           num_alleles,
#'                           fitness = function(){stats::rnorm(1, mean = 25)},
#'                           noise = function(){stats::rnorm(1, mean = 2)})
#'                           
#' print(landscape)
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
  fit_arr <- setup_matrix_to_array(fit_mat)
  params <- list(n_gene = n_gene,
                 n_allele = n_allele,
                 noise_func = noise,
                 type = "rmf")
  
  ## convert to FitLand object
  new_FitLand(fit_arr = fit_arr, params = params)
}

#####NK MODEL#####
# (need to fix Roxygen documentation below)
#' NK Model for Fitness Landscapes
#' 
#' Generates a fitness landscape using the NK model, as characterized in
#' Kauffman & Weinberger (1989) <doi:10.1016/s0022-5193(89)80019-0>.
#' 
#' @inheritParams generate_add
#' @param n number of genes
#' @param k number of other loci that influence fitness of each gene;
#'   determines landscape ruggedness
#' @param fitness if of class \code{array}, contains \code{numeric} fitness
#'   values for each (\code{k+1})-length genotype; if a \code{function}, takes
#'   a single parameter (\code{numeric} vector of length \code{k+1}) and
#'   returns a single \code{numeric} fitness value
#' @return fitness landscape stored in \code{FitLand} class object
#' @export
generate_nk <- function(n, k, n_allele = 2, fitness) {
  ## validate parameters
  # n - integer greater than 1
  if (!close_to_integer_limit(x = n, min_val = 1.5)) {
    stop(paste("parameter n must be an integer greater than 1, passed value =",
               n, "with class", class(n)),
         call. = FALSE)
  }
  
  # k - 0 <= k < n
  if (!close_to_integer_limit(x = k, min_val = -0.5, max_val = n - 0.5)) {
    stop(paste("parameter k must be an integer between 0 and n-1 (inclusive),",
               "passed value =", k, "with class", class(k)),
         call. = FALSE)
  }
  
  # n_allele
  if (!close_to_integer_limit(x = n_allele, min_val = 0.5)) {
    stop(paste("parameter n_allele must be a positive integer, passed value =",
               n_allele, "with class", class(n_allele)),
         call. = FALSE)
  }
  
  # fitness (array or function)
  if (is.function(fitness)) {
    if (length(formals(fitness)) != 1) {
      stop(paste("when fitness parameter is a function, it must have a single",
                 "parameter, passed function takes",
                 length(formals(fitness)), "parameters"),
           call. = FALSE)
    }
  } else if (is.array(fitness)) {
    if (!is.numeric(fitness)) {
      stop("when fitness parameter is an array, it must contain numeric",
           "values, passed array does not contain numeric values",
           call. = FALSE)
    } else if (!all.equal(dim(fitness), rep(n_allele, k + 1))) {
      stop(paste("when k =", k, "and n_allele =", n_allele, "fitness parameter",
                 "array must have dim =", rep(n_allele, k + 1), "but passed",
                 "array has dim =", dim(fitness)),
           call. = FALSE)
    }
  } else {
    stop(paste("fitness parameter must be a function or array, passed value",
               "has class =", class(fitness)),
         call. = FALSE)
  }
  
  ## setup fitness function
  fitness_func <- NULL
  if (!("function" %in% class(fitness))) {
    fitness_func <- function(index_vec) {
      fit_mat[t(index_vec)]
    }
  } else {
    fitness_func <- fitness
  }
  
  ## setup matrix with lattice coordinates
  fit_mat <- setup_matrix(n_dim = n, n_val = n_allele)
  val_col <- ncol(fit_mat)
  
  ## fill out last column w/ appropriate fitness values
  for (curr_row in seq_len(nrow(fit_mat))) {
    coords <- as.numeric(
      fit_mat[curr_row, seq_len(val_col - 1)]
    )
    
    counter_sum <- 0
    for (start in seq_len(length(coords))) {
      curr_coords <- nk_cycle_substr(index_vec = coords,
                                     start = start,
                                     len = k + 1)
      
      counter_sum <- counter_sum + fitness_func(curr_coords)
    }
    
    # store calculated fitness value
    fit_mat[curr_row, val_col] <- counter_sum
  }
  
  ## setup params for FitLand object
  fit_table <- setup_matrix_to_array(fit_mat)
  params <- list(n = n,
                 k = k,
                 n_allele = n_allele,
                 fit_calc = fitness_func,
                 type = "nk")
  if (is.array(fitness)) {
    params$fit_table <- fitness
  } else {
    params$fit_func <- fitness_func
  }
  
  ## convert to FitLand object
  new_FitLand(fit_arr = fit_table,
              params = params)
}

#####STICKBREAKING MODEL#####
# based on: Nagel et al. Genetics. 2012; 190: 655-667.
#' Stickbreaking Model for Fitness Landscapes
#' 
#' Generates a fitness landscape using the stickbreaking model, as
#' characterized in Nagel et al. (2012) <doi:10.1534.genetics.111.132134>.
#' 
#' @inheritParams generate_add
#' @param max_fit maximum fitness value attainable in this landscape
#' @return fitness landscape stored in \code{FitLand} class object
#' @export
#' @examples
#' n_gene <- 4
#' n_allele <- 2
#' wt_fit <- 1
#' max_fit <- 2
#' 
#' landscape <- generate_sb(n_gene = n_gene,
#'                          n_allele = n_allele,
#'                          wt_fit = wt_fit,
#'                          mut_fit = function() {stats::runif(1)},
#'                          max_fit = max_fit)
#'
#' print(landscape)
generate_sb <- function(n_gene, n_allele, wt_fit, mut_fit, max_fit) {
  ## validate parameters
  check_n_gene(n_gene)
  check_n_allele(n_allele)
  check_wt_fit(wt_fit)
  check_mut_fit(mut_fit, n_gene, n_allele)
  check_max_fit(max_fit, wt_fit)
  
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
  d <- max_fit - wt_fit
  for (curr_genotype in seq_len(nrow(temp_mat))) {
    # start at wild-type fitness
    fit <- wt_fit
    mult_counter <- 1
    
    # keep adjusting based on mutations
    for (curr_gene in seq(from = 1, to = n_gene, by = 1)) {
      mult_counter <- mult_counter *
        (1 - mut_matrix[temp_mat[curr_genotype, curr_gene], curr_gene])
    }
    
    # set value in setup matrix
    fit <- fit + d * (1 - mult_counter)
    temp_mat[curr_genotype, val_col] <- fit
  }
  
  ## setup params for FitLand object
  fit_table <- setup_matrix_to_array(temp_mat)
  params <- list(n_gene = n_gene,
                 n_allele = n_allele,
                 wt_fit = wt_fit,
                 max_fit = max_fit,
                 mut = mut_matrix,
                 type = "sb")
  
  ## convert to FitLand object
  new_FitLand(fit_arr = fit_table, params = params)
}

#####CORRELATED LANDSCAPES#####
generate_correlated <- function(...) {
  
}
