#' Transform dataset for cmprsn
#' @description transform_ds_for_cmprsn() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dataset for cmprsn. Function argument ds_tb specifies the object to be updated. Argument cmprsn_var_nm_1L_chr provides the object to be updated. The function returns Dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param cmprsn_var_nm_1L_chr Cmprsn var name (a character vector of length one)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'UID_chr'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'Timepoint_chr'
#' @param cmprsn_groups_chr Cmprsn groups (a character vector), Default: c("Intervention", "Control")
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
#' Transform dataset to predn dataset
#' @description transform_ds_to_predn_ds() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dataset to predn dataset. Function argument data_tb specifies the object to be updated. Argument predr_vars_nms_chr provides the object to be updated. The function returns Data (a tibble).
#' @param data_tb Data (a tibble)
#' @param predr_vars_nms_chr Predr vars names (a character vector)
#' @param tfmn_1L_chr Tfmn (a character vector of length one)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round bl value (a character vector of length one), Default: 'Baseline'
#' @param predictors_lup Predictors (a lookup table), Default: NULL
#' @return Data (a tibble)
#' @rdname transform_ds_to_predn_ds
#' @export 
#' @importFrom utils data
#' @importFrom dplyr mutate
#' @importFrom rlang sym exec
#' @importFrom purrr reduce map_dbl
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom TTU transform_tb_to_mdl_inp
#' @keywords internal
transform_ds_to_predn_ds <- function (data_tb, predr_vars_nms_chr, tfmn_1L_chr, dep_var_nm_1L_chr = "aqol6d_total_w", 
    id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round", 
    round_bl_val_1L_chr = "Baseline", predictors_lup = NULL) 
{
    if (is.null(predictors_lup)) 
        utils::data("predictors_lup", envir = environment())
    data_tb <- data_tb %>% dplyr::mutate(`:=`(!!rlang::sym(dep_var_nm_1L_chr), 
        NA_real_))
    data_tb <- purrr::reduce(predr_vars_nms_chr, .init = data_tb, 
        ~{
            predr_cls_fn <- eval(parse(text = ready4fun::get_from_lup_obj(predictors_lup, 
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = .y, 
                target_var_nm_1L_chr = "class_fn_chr", evaluate_lgl = F)))
            dplyr::mutate(.x, `:=`(!!rlang::sym(.y), !!rlang::sym(.y) %>% 
                rlang::exec(.fn = predr_cls_fn)))
        })
    data_tb <- data_tb %>% TTU::transform_tb_to_mdl_inp(dep_var_nm_1L_chr = dep_var_nm_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, id_var_nm_1L_chr = id_var_nm_1L_chr, 
        round_var_nm_1L_chr = round_var_nm_1L_chr, round_bl_val_1L_chr = round_bl_val_1L_chr, 
        drop_all_msng_1L_lgl = F, scaling_fctr_dbl = purrr::map_dbl(predr_vars_nms_chr, 
            ~ready4fun::get_from_lup_obj(predictors_lup, target_var_nm_1L_chr = "mdl_scaling_dbl", 
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = .x, 
                evaluate_lgl = F)), ungroup_1L_lgl = T, add_cll_tfmn_1L_lgl = ifelse(tfmn_1L_chr == 
            "CLL", T, F))
    return(data_tb)
}
