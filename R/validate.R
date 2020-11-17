# S3 generic (not exported)
check_fitland <- function(x, ...) {
  UseMethod("check_fitland")
}

# check validity of rough Mt. Fuji fitness landscape
check_fitland.RMF_FitLand <- function(x) {
  # check class
  
  # check dimensions (1 dim or more; all dims should be equal)
}

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

check_mut_fit <- function(mut_fit) {
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
}
