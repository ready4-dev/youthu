update_col_with_diff <- function(ds_tb,
                                 var_nm_1L_chr,
                                 fn,
                                 abs_mean_diff_1L_dbl, # NB - Population mean, not sample mean
                                 diff_sd_1L_dbl,
                                 multiplier_1L_dbl,
                                 min_1L_dbl,
                                 max_1L_dbl,
                                 integer_1L_lgl){
  args_ls <- append(list(nrow(ds_tb)), list(abs_mean_diff_1L_dbl,
                                            diff_sd_1L_dbl)) %>% unname()
  diff_dbl <- (rlang::exec(.fn = fn, !!!args_ls))*multiplier_1L_dbl
  new_ds <- ds_tb %>%
    dplyr::mutate(!!rlang::sym(var_nm_1L_chr) := !!rlang::sym(var_nm_1L_chr) %>% purrr::map2_dbl(diff_dbl, ~{
      new_total_1L_dbl <- .x + .y
      new_total_1L_dbl <- min(max(new_total_1L_dbl,min_1L_dbl),max_1L_dbl)
      new_total_1L_xx <- ifelse(integer_1L_lgl, as.integer(new_total_1L_dbl), new_total_1L_dbl)
      new_total_1L_xx
    })
    )
  return(new_ds)
}
update_multpl_cols_with_diffs <- function(ds_tb,
                                          var_nms_chr,
                                          fns_ls,
                                          abs_mean_diff_dbl, # NB - Population mean, not sample mean
                                          diff_sd_dbl,
                                          multiplier_dbl,
                                          min_dbl,
                                          max_dbl,
                                          integer_lgl){

  args_ls_ls <- list(var_nms_chr = var_nms_chr,
                     fns_ls = fns_ls,
                     abs_mean_diff_dbl = abs_mean_diff_dbl,
                     diff_sd_dbl = diff_sd_dbl,
                     multiplier_dbl = multiplier_dbl,
                     min_dbl = min_dbl,
                     max_dbl = max_dbl,
                     integer_lgl = integer_lgl)
  updated_ds_tb  <- 1:length(args_ls_ls[[1]]) %>%
    purrr::reduce(.init = ds_tb,
                  ~{
                    idx_1L_int <- .y
                    args_ls <- args_ls_ls %>% purrr::map(~.x %>% purrr::pluck(idx_1L_int)) %>% unname()
                    rlang::exec(.fn = update_col_with_diff,
                                ds_tb = .x,
                                !!!args_ls)
                  })
  return(updated_ds_tb)
}
