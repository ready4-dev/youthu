#' Predictors lookup table
#' 
#' A lookup table of the short name and long name of each predictor used in the models included with the youthu package.
#' 
#' A tibble
#' 
#' \describe{
#'   \item{short_name_chr}{Short name (a character vector)}
#'   \item{long_name_chr}{Long name (a character vector)}
#'   \item{min_val_dbl}{Min value (a double vector)}
#'   \item{max_val_dbl}{Max value (a double vector)}
#'   \item{class_chr}{Class (a character vector)}
#'   \item{increment_dbl}{Increment (a double vector)}
#'   \item{class_fn_chr}{Class function (a character vector)}
#'   \item{mdl_scaling_dbl}{Mdl scaling (a double vector)}
#'   \item{covariate_lgl}{Covariate (a logical vector)}
#' }
"predictors_lup"
