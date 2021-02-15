#' Calculate rmse tfmn
#' @description calculate_rmse_tfmn() is a Calculate function that calculates a numeric value. Specifically, this function implements an algorithm to calculate rmse tfmn. The function returns Rmse tfmn (a double vector).
#' @param y_dbl Y (a double vector)
#' @param yhat_dbl Yhat (a double vector)
#' @return Rmse tfmn (a double vector)
#' @rdname calculate_rmse_tfmn
#' @export 

#' @keywords internal
calculate_rmse_tfmn <- function (y_dbl, yhat_dbl) 
{
    rmse_tfmn_dbl <- sqrt(mean((1 - exp(-exp(yhat_dbl)) - y_dbl)^2))
    return(rmse_tfmn_dbl)
}
