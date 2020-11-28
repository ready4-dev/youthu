impute_adult_aqol6d_items_tb <- function (aqol6d_items_tb, domain_items_ls) 
{
    aqol6d_items_tb <- 1:length(domain_items_ls) %>% purrr::reduce(.init = aqol6d_items_tb, 
        ~{
            idx_int <- .y
            tb <- .x %>% dplyr::mutate(`:=`(!!rlang::sym(paste0("D", 
                names(domain_items_ls)[idx_int], "miss")), rowMeans(dplyr::select(., 
                domain_items_ls[[idx_int]]), na.rm = TRUE)))
            tb %>% dplyr::mutate_at(dplyr::vars(domain_items_ls[[idx_int]]), 
                ~.x %>% purrr::map2_dbl(!!rlang::sym(paste0("D", 
                  names(domain_items_ls)[idx_int], "miss")), 
                  ~{
                    ifelse(is.na(.x), .y, .x)
                  }))
            tb
        })
    return(aqol6d_items_tb)
}
impute_unscrd_adol_aqol6d_ds <- function (unscrd_aqol_ds_tb) 
{
    unscrd_aqol_ds_tb <- unscrd_aqol_ds_tb %>% dplyr::mutate(missing = rowSums(is.na(dplyr::select(., 
        paste0("Q", c(1:10))))))
    aqol_cases_to_imp_tb <- unscrd_aqol_ds_tb %>% dplyr::filter(missing < 
        10) %>% dplyr::select(-missing)
    aqol_cases_not_to_imp_tb <- unscrd_aqol_ds_tb %>% dplyr::filter(missing >= 
        10) %>% dplyr::select(-missing)
    imputed_aqol_tb <- mice::mice(aqol_cases_to_imp_tb, m = 1, 
        maxit = 50, meth = "pmm", seed = 1234)
    aqol_cases_to_imp_tb <- mice::complete(imputed_aqol_tb, "long") %>% 
        dplyr::select(-.imp, -.id)
    imputed_unscrd_aqol_ds_tb_tb <- data.frame(rbind(aqol_cases_to_imp_tb, 
        aqol_cases_not_to_imp_tb))
    return(imputed_unscrd_aqol_ds_tb_tb)
}
