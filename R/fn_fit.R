#' Fit ts model with brm
#' @description fit_ts_model_with_brm() is a Fit function that fits a model of a specified type to a dataset Specifically, this function implements an algorithm to fit ts model with brm. The function returns Mdl (a list).
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one)
#' @param predr_vars_nms_chr Predr vars names (a character vector)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one)
#' @param backend_1L_chr Backend (a character vector of length one), Default: getOption("brms.backend", "rstan")
#' @param link_1L_chr Link (a character vector of length one), Default: 'identity'
#' @param iters_1L_int Iters (an integer vector of length one), Default: 4000
#' @param seed_1L_int Seed (an integer vector of length one), Default: 1000
#' @return Mdl (a list)
#' @rdname fit_ts_model_with_brm
#' @export 
#' @importFrom brms brm
#' @importFrom purrr map_chr
#' @keywords internal
fit_ts_model_with_brm <- function (data_tb, dep_var_nm_1L_chr, predr_vars_nms_chr, id_var_nm_1L_chr, 
    backend_1L_chr = getOption("brms.backend", "rstan"), link_1L_chr = "identity", 
    iters_1L_int = 4000L, seed_1L_int = 1000L) 
{
    mdl_ls <- brms::brm(formula = as.formula(paste0(dep_var_nm_1L_chr, 
        " ~ ", purrr::map_chr(predr_vars_nms_chr, ~paste0(.x, 
            "_baseline + ", .x, "_change + ")) %>% paste0(collapse = ""), 
        "(1|", id_var_nm_1L_chr, ")")), backend = backend_1L_chr, 
        data = data_tb, family = gaussian(link = link_1L_chr), 
        iter = iters_1L_int, seed = seed_1L_int)
    return(mdl_ls)
}
