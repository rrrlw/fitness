#####CONSTRUCTOR & VALIDATOR#####
# FitLand constructor
new_FitLand <- function(fit_table, params) {
  # assume all parameters are valid
  ans <- list()
  ans$fit_table <- fit_table
  ans$params <- params
  
  # return FitLand object
  structure(ans, class = "FitLand")
}

# FitLand validator
validate_FitLand <- function(x) {
  # which elements is params guaranteed to have?
  # - dims (integer vector of dimensions)
  # - type (character vector naming type of landscape)
  # - other params specific to type
  
  # match fit_table dimensions to params values
  
}

#####HELPERS#####
# add AFTER generators have been built (just a wrapper for them)