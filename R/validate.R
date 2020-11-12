# S3 generic (not exported)
check_fitland <- function(x, ...) {
  UseMethod("check_fitland")
}

# check validity of rough Mt. Fuji fitness landscape
check_fitland.RMF_FitLand <- function(x) {
  # check class
  
  # check dimensions (1 dim or more; all dims should be equal)
}
