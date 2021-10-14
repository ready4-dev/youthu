predict_from_mdl_coefs <- function (smry_of_mdl_tb, new_data_tb) 
{
    coef_tb <- smry_of_mdl_tb %>% dplyr::filter(!Parameter %in% 
        c("R2", "RMSE", "Sigma"))
    vecs_1_ls <- coef_tb$Parameter[-1] %>% purrr::map(~new_data_tb %>% 
        dplyr::pull(.x %>% stringr::str_replace(" ", "_")) * 
        coef_tb %>% dplyr::filter(Parameter == .x) %>% dplyr::pull(Estimate))
    predd_dbl <- exp(Reduce(`+`, vecs_1_ls) + coef_tb %>% dplyr::filter(Parameter == 
        "Intercept") %>% dplyr::pull(Estimate))
    return(predd_dbl)
}
