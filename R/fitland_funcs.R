#####DEPRECATED AFTER NEW FITLAND STRUCTURE#####
# S3 generic (should be exported)
# get dimension of a fitness landscape
get_dim <- function(x, ...) {
  UseMethod("get_dim")
}

# for rough Mt. Fuji
get_dim.RMF_FitLand <- function(x) {
  return(length(dim(x$fits)))
}

# S3 generic (should be exported)
# get alphabet length of a fitness landscape
get_nalpha <- function(x, ...) {
  UseMethod("get_nalpha")
}

get_nalpha.RMF_FitLand <- function(x) {
  return(dim(x$fits)[1])
}

# S3 generic (should be exported)
# get alphabet for a fitness landscape
get_alphabet <- function(x, ...) {
  UseMethod("get_alphabet")
}

get_alphabet.RMF_FitLand <- function(x) {
  return(dimnames(x$fits)[[1]])
}
