force_min_max_and_int_cnstrs <- function (tbl_tb, var_names_chr, min_max_ls, discrete_lgl) 
{
    tbl_tb <- purrr::reduce(1:length(var_names_chr), .init = tbl_tb, 
        ~{
            idx_dbl <- .y
            .x %>% dplyr::mutate(`:=`(!!rlang::sym(var_names_chr[.y]), 
                !!rlang::sym(var_names_chr[.y]) %>% purrr::map_dbl(~{
                  min(ifelse(discrete_lgl[idx_dbl], round(.x), 
                    .x), min_max_ls[[idx_dbl]][2]) %>% max(min_max_ls[[idx_dbl]][1])
                })))
        })
    return(tbl_tb)
}
force_vec_to_sum_to_int <- function (vec_int, target_1L_int, item_ranges_dbl_ls) 
{
    extras_int <- target_1L_int - sum(vec_int)
    if (extras_int != 0) {
        increment_int <- ifelse(extras_int > 0, 1, -1)
        vec_int <- purrr::reduce(1:abs(extras_int), .init = vec_int, 
            ~{
                new_vect_int <- .x
                idx_lgl <- purrr::map2_lgl(new_vect_int, item_ranges_dbl_ls, 
                  ~{
                    (.y[1] <= (.x + increment_int)) & ((.x + 
                      increment_int) <= .y[2])
                  })
                possible_idx_dbl <- (1:length(new_vect_int))[idx_lgl]
                idx_dbl <- ifelse(length(possible_idx_dbl) == 
                  1, possible_idx_dbl, sample(possible_idx_dbl, 
                  size = 1))
                new_vect_int[idx_dbl] <- new_vect_int[idx_dbl] + 
                  increment_int
                new_vect_int
            })
    }
    return(vec_int)
}
