#' Reorder cndt predrs
#' @description reorder_cndt_predrs_chr() is a Reorder function that reorders an object to conform to a pre-specified schema. Specifically, this function implements an algorithm to reorder cndt predrs character vector. The function is called for its side effects and does not return a value.
#' @param candidate_predrs_chr Candidate predrs (a character vector)
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param method_1L_chr Method (a character vector of length one), Default: 'pearson'
#' @return NA ()
#' @rdname reorder_cndt_predrs_chr
#' @export 
#' @importFrom dplyr select all_of mutate across arrange desc filter pull
#' @importFrom Hmisc rcorr
#' @importFrom tibble as_tibble
#' @importFrom rlang sym
#' @keywords internal
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
#' Reorder tibbles for target cors
#' @description reorder_tbs_for_target_cors() is a Reorder function that reorders an object to conform to a pre-specified schema. Specifically, this function implements an algorithm to reorder tibbles for target cors. The function returns Tibbles (a list).
#' @param tbs_ls Tibbles (a list)
#' @param cor_dbl Cor (a double vector)
#' @param cor_var_chr Cor var (a character vector)
#' @param id_var_to_rm_1L_chr Id var to rm (a character vector of length one), Default: 'NA'
#' @return Tibbles (a list)
#' @rdname reorder_tbs_for_target_cors
#' @export 
#' @importFrom MASS mvrnorm
#' @importFrom dplyr slice arrange bind_rows select
#' @importFrom rlang sym
#' @importFrom purrr map
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
