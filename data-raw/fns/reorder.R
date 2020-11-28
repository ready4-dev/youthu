reorder_cndt_predrs_chr <- function (candidate_predrs_chr, data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", 
    method_1L_chr = "pearson") 
{
    data_mat <- as.matrix(data_tb %>% dplyr::select(c(dplyr::all_of(dep_var_nm_1L_chr), 
        dplyr::all_of(candidate_predrs_chr))))
    corr_ls <- Hmisc::rcorr(data_mat, type = method_1L_chr)
    reordered_cndt_predrs <- corr_ls$r %>% tibble::as_tibble(rownames = "var_nms_chr") %>% 
        dplyr::mutate(dplyr::across(where(is.numeric), abs)) %>% 
        dplyr::arrange(dplyr::desc(!!rlang::sym(dep_var_nm_1L_chr))) %>% 
        dplyr::filter(var_nms_chr != dep_var_nm_1L_chr) %>% dplyr::pull(var_nms_chr) %>% 
        as.vector()
    return(reordered_cndt_predrs)
}
reorder_tbs_for_target_cors <- function (tbs_ls, cor_dbl, cor_var_chr, id_var_to_rm_1L_chr = NA_character_) 
{
    n_fup_dbl <- nrow(tbs_ls[[2]])
    cor_mat <- matrix(cor_dbl, ncol = 2, nrow = 2)
    diag(cor_mat) <- 1
    mvdat_mat <- MASS::mvrnorm(n = n_fup_dbl, mu = c(0, 0), Sigma = cor_mat, 
        empirical = TRUE)
    rank_x_int <- rank(mvdat_mat[, 1], ties.method = "first")
    rank_y_int <- rank(mvdat_mat[, 2], ties.method = "first")
    matched_tb <- tbs_ls[[1]] %>% dplyr::slice(tbs_ls[[2]]$id) %>% 
        dplyr::arrange(!!rlang::sym(cor_var_chr[1])) %>% dplyr::slice(rank_x_int)
    unmatched_int <- setdiff(1:nrow(tbs_ls[[1]]), tbs_ls[[2]]$id)
    if (!identical(unmatched_int, character(0))) {
        tbs_ls[[1]] <- dplyr::bind_rows(matched_tb, tbs_ls[[1]] %>% 
            dplyr::slice(unmatched_int))
    }
    else {
        tbs_ls[[1]] <- matched_tb
    }
    tbs_ls[[2]] <- tbs_ls[[2]] %>% dplyr::arrange(!!rlang::sym(cor_var_chr[2])) %>% 
        dplyr::slice(rank_y_int)
    if (!is.na(id_var_to_rm_1L_chr)) 
        tbs_ls <- tbs_ls %>% purrr::map(~.x %>% dplyr::select(-!!rlang::sym(id_var_to_rm_1L_chr)))
    return(tbs_ls)
}
