#' Force min max and integer vector cnstrs
#' @description force_min_max_and_int_cnstrs() is a Force function that checks if a specified local or global environmental condition is met and if not, updates the specified environment to comply with the condition. Specifically, this function implements an algorithm to force min max and integer vector cnstrs. The function returns Table (a tibble).
#' @param tbl_tb Table (a tibble)
#' @param var_names_chr Var names (a character vector)
#' @param min_max_ls Min max (a list)
#' @param discrete_lgl Discrete (a logical vector)
#' @return Table (a tibble)
#' @rdname force_min_max_and_int_cnstrs
#' @export 
#' @importFrom purrr reduce map_dbl
#' @importFrom dplyr mutate
#' @importFrom rlang sym
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
#' Force vec to sum to
#' @description force_vec_to_sum_to_int() is a Force function that checks if a specified local or global environmental condition is met and if not, updates the specified environment to comply with the condition. Specifically, this function implements an algorithm to force vec to sum to integer vector. The function returns Vec (an integer vector).
#' @param vec_int Vec (an integer vector)
#' @param target_1L_int Target (an integer vector of length one)
#' @param item_ranges_dbl_ls Item ranges (a list of double vectors)
#' @return Vec (an integer vector)
#' @rdname force_vec_to_sum_to_int
#' @export 
#' @importFrom purrr reduce map2_lgl
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
