#' Transform dataset for comparison
#' @description transform_ds_for_cmprsn() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dataset for comparison. Function argument ds_tb specifies the object to be updated. Argument cmprsn_var_nm_1L_chr provides the object to be updated. The function returns Dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param cmprsn_var_nm_1L_chr Comparison variable name (a character vector of length one)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'UID_chr'
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'Timepoint_chr'
#' @param cmprsn_groups_chr Comparison groups (a character vector), Default: c("Intervention", "Control")
#' @return Dataset (a tibble)
#' @rdname transform_ds_for_cmprsn
#' @export 
#' @importFrom dplyr group_by mutate n ungroup filter pull case_when select
#' @importFrom rlang sym
#' @keywords internal
transform_ds_for_cmprsn <- function (ds_tb, cmprsn_var_nm_1L_chr, id_var_nm_1L_chr = "UID_chr", 
    round_var_nm_1L_chr = "Timepoint_chr", cmprsn_groups_chr = c("Intervention", 
        "Control")) 
{
    ds_tb <- ds_tb %>% dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) %>% 
        dplyr::mutate(n_measures_int = dplyr::n()) %>% dplyr::ungroup() %>% 
        dplyr::filter(n_measures_int == 2)
    UIDs_with_FUP_data_chr <- ds_tb %>% dplyr::pull(!!rlang::sym(id_var_nm_1L_chr)) %>% 
        unique()
    group_1_IDs <- sample(UIDs_with_FUP_data_chr, floor(length(UIDs_with_FUP_data_chr)/2))
    group_2_IDs <- setdiff(UIDs_with_FUP_data_chr, group_1_IDs)[1:length(group_1_IDs)]
    ds_tb <- ds_tb %>% dplyr::mutate(`:=`(!!rlang::sym(cmprsn_var_nm_1L_chr), 
        dplyr::case_when(!!rlang::sym(id_var_nm_1L_chr) %in% 
            group_1_IDs ~ cmprsn_groups_chr[1], !!rlang::sym(id_var_nm_1L_chr) %in% 
            group_2_IDs ~ cmprsn_groups_chr[2], T ~ NA_character_))) %>% 
        dplyr::select(-n_measures_int) %>% dplyr::filter(!!rlang::sym(cmprsn_var_nm_1L_chr) %in% 
        cmprsn_groups_chr)
    return(ds_tb)
}
