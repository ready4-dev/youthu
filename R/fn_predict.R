#' Predict from mdl coefs
#' @description predict_from_mdl_coefs() is a Predict function that makes predictions from data using a specified statistical model. Specifically, this function implements an algorithm to predict from mdl coefs. The function returns Pred (a double vector).
#' @param smry_of_mdl_tb Smry of mdl (a tibble)
#' @param new_data_tb New data (a tibble)
#' @return Pred (a double vector)
#' @rdname predict_from_mdl_coefs
#' @export 
#' @importFrom dplyr filter pull
#' @importFrom purrr map
#' @importFrom stringr str_replace
#' @keywords internal
predict_from_mdl_coefs <- function (smry_of_mdl_tb, new_data_tb) 
{
    coef_tb <- smry_of_mdl_tb %>% dplyr::filter(!Parameter %in% 
        c("R2", "RMSE", "Sigma"))
    vecs_1_ls <- coef_tb$Parameter[-1] %>% purrr::map(~new_data_tb %>% 
        dplyr::pull(.x %>% stringr::str_replace(" ", "_")) * 
        coef_tb %>% dplyr::filter(Parameter == .x) %>% dplyr::pull(Estimate))
    pred_dbl <- exp(Reduce(`+`, vecs_1_ls) + coef_tb %>% dplyr::filter(Parameter == 
        "Intercept") %>% dplyr::pull(Estimate))
    return(pred_dbl)
}
