calculate_rmse_tfmn <- function (y_dbl, yhat_dbl)
{
    rmse_tfmn_dbl <- sqrt(mean((1 - exp(-exp(yhat_dbl)) - y_dbl)^2))
    return(rmse_tfmn_dbl)
}
