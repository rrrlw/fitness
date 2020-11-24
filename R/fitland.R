#####CONSTRUCTOR & VALIDATOR#####
# FitLand constructor
new_FitLand <- function(fit_arr, params) {
  # assume all parameters are valid
  ans <- list()
  ans$fit_arr <- fit_arr
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

#####HELPER(S)#####
FitLand <- function(fit_array, params = list()) {
  validate_FitLand(
    new_FitLand(fit_array, params)
  )
}
