#' Update column with difference
#' @description update_col_with_diff() is an Update function that edits an object, while preserving core object attributes. Specifically, this function implements an algorithm to update column with difference. The function is called for its side effects and does not return a value.
#' @param ds_tb Dataset (a tibble)
#' @param var_nm_1L_chr Variable name (a character vector of length one)
#' @param fn Function (a function)
#' @param abs_mean_diff_1L_dbl Absolute mean difference (a double vector of length one)
#' @param diff_sd_1L_dbl Difference standard deviation (a double vector of length one)
#' @param multiplier_1L_dbl Multiplier (a double vector of length one)
#' @param min_1L_dbl Minimum (a double vector of length one)
#' @param max_1L_dbl Maximum (a double vector of length one)
#' @param integer_1L_lgl Integer (a logical vector of length one)
#' @return New (a dataset)
#' @rdname update_col_with_diff
#' @export 
#' @importFrom rlang exec sym
#' @importFrom dplyr mutate
#' @importFrom purrr map2_dbl
#' @keywords internal
update_col_with_diff <- function (ds_tb, var_nm_1L_chr, fn, abs_mean_diff_1L_dbl, diff_sd_1L_dbl, 
    multiplier_1L_dbl, min_1L_dbl, max_1L_dbl, integer_1L_lgl) 
{
    args_ls <- append(list(nrow(ds_tb)), list(abs_mean_diff_1L_dbl, 
        diff_sd_1L_dbl)) %>% unname()
    diff_dbl <- (rlang::exec(.fn = fn, !!!args_ls)) * multiplier_1L_dbl
    new_ds <- ds_tb %>% dplyr::mutate(`:=`(!!rlang::sym(var_nm_1L_chr), 
        !!rlang::sym(var_nm_1L_chr) %>% purrr::map2_dbl(diff_dbl, 
            ~{
                new_total_1L_dbl <- .x + .y
                new_total_1L_dbl <- min(max(new_total_1L_dbl, 
                  min_1L_dbl), max_1L_dbl)
                new_total_1L_xx <- ifelse(integer_1L_lgl, as.integer(new_total_1L_dbl), 
                  new_total_1L_dbl)
                new_total_1L_xx
            })))
    return(new_ds)
}
#' Update multiplier columns with differences
#' @description update_multpl_cols_with_diffs() is an Update function that edits an object, while preserving core object attributes. Specifically, this function implements an algorithm to update multiplier columns with differences. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param var_nms_chr Variable names (a character vector)
#' @param fns_ls Functions (a list)
#' @param abs_mean_diff_dbl Absolute mean difference (a double vector)
#' @param diff_sd_dbl Difference standard deviation (a double vector)
#' @param multiplier_dbl Multiplier (a double vector)
#' @param min_dbl Minimum (a double vector)
#' @param max_dbl Maximum (a double vector)
#' @param integer_lgl Integer (a logical vector)
#' @return Updated dataset (a tibble)
#' @rdname update_multpl_cols_with_diffs
#' @export 
#' @importFrom purrr reduce map pluck
#' @importFrom rlang exec
#' @keywords internal
update_multpl_cols_with_diffs <- function (ds_tb, var_nms_chr, fns_ls, abs_mean_diff_dbl, diff_sd_dbl, 
    multiplier_dbl, min_dbl, max_dbl, integer_lgl) 
{
    args_ls_ls <- list(var_nms_chr = var_nms_chr, fns_ls = fns_ls, 
        abs_mean_diff_dbl = abs_mean_diff_dbl, diff_sd_dbl = diff_sd_dbl, 
        multiplier_dbl = multiplier_dbl, min_dbl = min_dbl, max_dbl = max_dbl, 
        integer_lgl = integer_lgl)
    updated_ds_tb <- 1:length(args_ls_ls[[1]]) %>% purrr::reduce(.init = ds_tb, 
        ~{
            idx_1L_int <- .y
            args_ls <- args_ls_ls %>% purrr::map(~.x %>% purrr::pluck(idx_1L_int)) %>% 
                unname()
            rlang::exec(.fn = update_col_with_diff, ds_tb = .x, 
                !!!args_ls)
        })
    return(updated_ds_tb)
}
