#' Rename from nmd vec
#' @description rename_from_nmd_vec() is a Rename function that renames elements of an object based on a pre-speccified schema. Specifically, this function implements an algorithm to rename from nmd vec. The function returns Renamed data (a tibble).
#' @param data_tb Data (a tibble)
#' @param nmd_vec_chr Nmd vec (a character vector)
#' @param vec_nms_as_new_1L_lgl Vec names as new (a logical vector of length one), Default: T
#' @return Renamed data (a tibble)
#' @rdname rename_from_nmd_vec
#' @export 
#' @importFrom purrr reduce
#' @importFrom dplyr rename
#' @importFrom rlang sym
rename_from_nmd_vec <- function (data_tb, nmd_vec_chr, vec_nms_as_new_1L_lgl = T) 
{
    if (vec_nms_as_new_1L_lgl) {
        renamed_data_tb <- purrr::reduce(1:length(nmd_vec_chr), 
            .init = data_tb, ~dplyr::rename(.x, `:=`(!!rlang::sym(names(nmd_vec_chr)[.y]), 
                nmd_vec_chr[.y] %>% as.vector())))
    }
    else {
        renamed_data_tb <- purrr::reduce(1:length(nmd_vec_chr), 
            .init = data_tb, ~dplyr::rename(.x, `:=`(!!rlang::sym(nmd_vec_chr[.y]), 
                names(nmd_vec_chr)[.y] %>% as.vector())))
    }
    return(renamed_data_tb)
}
