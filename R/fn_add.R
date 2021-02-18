#' Add aqol6d predn to dataset
#' @description add_aqol6d_predn_to_ds() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add aqol6d predn to dataset. Function argument data_tb specifies the object to be updated. The function returns Updated (a tibble).
#' @param data_tb Data (a tibble)
#' @param model_mdl PARAM_DESCRIPTION
#' @param tfmn_1L_chr Tfmn (a character vector of length one)
#' @param predr_vars_nms_chr Predr vars names (a character vector), Default: NULL
#' @param utl_var_nm_1L_chr Utl var name (a character vector of length one), Default: NULL
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round bl value (a character vector of length one), Default: 'Baseline'
#' @param utl_cls_fn Utl class (a function), Default: firstbounce_aqol6d_adol
#' @param predictors_lup Predictors (a lookup table), Default: NULL
#' @return Updated (a tibble)
#' @rdname add_aqol6d_predn_to_ds
#' @export 
#' @importFrom utils data
#' @importFrom purrr flatten_chr map_chr
#' @importFrom stringr str_replace
#' @importFrom TTU add_utility_predn_to_ds
#' @importFrom dplyr rename select left_join
#' @importFrom rlang sym
#' @importFrom tidyselect all_of
add_aqol6d_predn_to_ds <- function (data_tb, model_mdl, tfmn_1L_chr, predr_vars_nms_chr = NULL, 
    utl_var_nm_1L_chr = NULL, id_var_nm_1L_chr = "fkClientID", 
    round_var_nm_1L_chr = "round", round_bl_val_1L_chr = "Baseline", 
    utl_cls_fn = firstbounce_aqol6d_adol, predictors_lup = NULL) 
{
    if (is.null(predictors_lup)) 
        utils::data("predictors_lup", envir = environment())
    if (!is.null(names(predr_vars_nms_chr))) {
        data_tb <- rename_from_nmd_vec(data_tb, nmd_vec_chr = predr_vars_nms_chr, 
            vec_nms_as_new_1L_lgl = T)
    }
    terms_ls <- model_mdl$terms
    mdl_dep_var_1L_chr <- terms_ls[[2]] %>% as.character()
    mdl_predr_terms_chr <- terms_ls[[3]] %>% as.character()
    mdl_predr_terms_chr <- mdl_predr_terms_chr %>% strsplit(split = " +") %>% 
        purrr::flatten_chr()
    mdl_predr_terms_chr <- mdl_predr_terms_chr[mdl_predr_terms_chr != 
        "+"]
    mdl_predr_terms_chr <- mdl_predr_terms_chr %>% purrr::map_chr(~stringr::str_replace(.x, 
        "_baseline", "") %>% stringr::str_replace("_change", 
        "")) %>% unique()
    original_ds_vars_chr <- names(data_tb)[!names(data_tb) %in% 
        c(mdl_predr_terms_chr, ifelse(!is.null(utl_var_nm_1L_chr), 
            utl_var_nm_1L_chr, mdl_dep_var_1L_chr))]
    updated_tb <- data_tb %>% transform_ds_to_predn_ds(predr_vars_nms_chr = mdl_predr_terms_chr, 
        tfmn_1L_chr = tfmn_1L_chr, dep_var_nm_1L_chr = mdl_dep_var_1L_chr, 
        id_var_nm_1L_chr = id_var_nm_1L_chr, round_var_nm_1L_chr = round_var_nm_1L_chr, 
        round_bl_val_1L_chr = round_bl_val_1L_chr, predictors_lup = predictors_lup) %>% 
        TTU::add_utility_predn_to_ds(model_mdl = model_mdl, tfmn_1L_chr = tfmn_1L_chr, 
            dep_var_nm_1L_chr = mdl_dep_var_1L_chr, predr_vars_nms_chr = mdl_predr_terms_chr, 
            utl_cls_fn = firstbounce_aqol6d_adol, rmv_tfmd_dep_var_1L_lgl = T)
    if (!is.null(utl_var_nm_1L_chr)) {
        updated_tb <- updated_tb %>% dplyr::rename(`:=`(!!rlang::sym(utl_var_nm_1L_chr), 
            tidyselect::all_of(mdl_dep_var_1L_chr)))
    }
    if (!is.null(names(predr_vars_nms_chr))) {
        updated_tb <- rename_from_nmd_vec(updated_tb, nmd_vec_chr = predr_vars_nms_chr, 
            vec_nms_as_new_1L_lgl = F)
    }
    if ("aqol6d_total_w_CLL_cloglog" %in% names(updated_tb)) 
        updated_tb <- updated_tb %>% dplyr::select(-aqol6d_total_w_CLL_cloglog)
    names_to_incl_chr <- c(names(updated_tb), setdiff(names(data_tb), 
        names(updated_tb)))
    updated_tb <- dplyr::left_join(data_tb %>% dplyr::select(tidyselect::all_of(original_ds_vars_chr)), 
        updated_tb)
    updated_tb <- updated_tb %>% dplyr::select(tidyselect::all_of(names_to_incl_chr[names_to_incl_chr %in% 
        names(updated_tb)]))
    return(updated_tb)
}
#' Add change in dataset var
#' @description add_change_in_ds_var() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add change in dataset var. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round bl value (a character vector of length one), Default: 'Baseline'
#' @param change_var_nm_1L_chr Change var name (a character vector of length one)
#' @param var_nm_1L_chr Var name (a character vector of length one)
#' @param arrange_by_id_lgl Arrange by id (a logical vector), Default: T
#' @return Updated dataset (a tibble)
#' @rdname add_change_in_ds_var
#' @export 
#' @importFrom dplyr group_by mutate case_when lag ungroup arrange
#' @importFrom rlang sym
#' @keywords internal
add_change_in_ds_var <- function (ds_tb, id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round", 
    round_bl_val_1L_chr = "Baseline", change_var_nm_1L_chr, var_nm_1L_chr, 
    arrange_by_id_lgl = T) 
{
    updated_ds_tb <- ds_tb %>% dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) %>% 
        dplyr::mutate(`:=`(!!rlang::sym(change_var_nm_1L_chr), 
            dplyr::case_when(!!rlang::sym(round_var_nm_1L_chr) == 
                round_bl_val_1L_chr ~ 0, T ~ (as.numeric(!!rlang::sym(var_nm_1L_chr)) - 
                dplyr::lag(as.numeric(!!rlang::sym(var_nm_1L_chr))))))) %>% 
        dplyr::ungroup()
    if (arrange_by_id_lgl) 
        updated_ds_tb <- updated_ds_tb %>% dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr))
    return(updated_ds_tb)
}
#' Add costs by tmpt
#' @description add_costs_by_tmpt() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add costs by tmpt. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param round_var_nm_1L_chr Round var name (a character vector of length one)
#' @param round_lvls_chr Round lvls (a character vector), Default: c("Baseline", "Follow-up")
#' @param costs_mean_dbl Costs mean (a double vector)
#' @param costs_sd_dbl Costs sd (a double vector)
#' @param extra_cost_args_ls Extra cost arguments (a list), Default: list(costs_var_nm_1L_chr = "costs_dbl")
#' @param fn Function (a function), Default: add_costs_from_gamma_dist
#' @return Updated dataset (a tibble)
#' @rdname add_costs_by_tmpt
#' @export 
#' @importFrom purrr pmap_dfr
#' @importFrom dplyr filter
#' @importFrom rlang sym exec
#' @keywords internal
add_costs_by_tmpt <- function (ds_tb, round_var_nm_1L_chr, round_lvls_chr = c("Baseline", 
    "Follow-up"), costs_mean_dbl, costs_sd_dbl, extra_cost_args_ls = list(costs_var_nm_1L_chr = "costs_dbl"), 
    fn = add_costs_from_gamma_dist) 
{
    updated_ds_tb <- purrr::pmap_dfr(list(round_lvls_chr, costs_mean_dbl, 
        costs_sd_dbl), ~{
        args_ls <- list(ds_tb %>% dplyr::filter(!!rlang::sym(round_var_nm_1L_chr) == 
            ..1), ..2, ..3)
        if (!is.null(extra_cost_args_ls)) 
            args_ls <- append(args_ls, extra_cost_args_ls)
        rlang::exec(.fn = fn, !!!args_ls)
    })
    return(updated_ds_tb)
}
#' Add costs from gamma dist
#' @description add_costs_from_gamma_dist() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add costs from gamma dist. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param costs_mean_dbl Costs mean (a double vector)
#' @param costs_sd_dbl Costs sd (a double vector)
#' @param costs_var_nm_1L_chr Costs var name (a character vector of length one), Default: 'costs_dbl'
#' @return Updated dataset (a tibble)
#' @rdname add_costs_from_gamma_dist
#' @export 
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @keywords internal
add_costs_from_gamma_dist <- function (ds_tb, costs_mean_dbl, costs_sd_dbl, costs_var_nm_1L_chr = "costs_dbl") 
{
    updated_ds_tb <- dplyr::mutate(ds_tb, `:=`(!!rlang::sym(costs_var_nm_1L_chr), 
        make_costs_vec_from_gamma_dist(n_int = nrow(ds_tb), costs_mean_dbl = costs_mean_dbl, 
            costs_sd_dbl = costs_sd_dbl)))
    return(updated_ds_tb)
}
#' Add dates from dist
#' @description add_dates_from_dist() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add dates from dist. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param bl_start_date_dtm PARAM_DESCRIPTION
#' @param bl_end_date_dtm PARAM_DESCRIPTION
#' @param duration_args_ls Duration arguments (a list)
#' @param duration_fn Duration (a function), Default: rnorm
#' @param date_var_nm_1L_chr Date var name (a character vector of length one), Default: 'date_psx'
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round bl value (a character vector of length one), Default: 'Baseline'
#' @param origin_1L_chr Origin (a character vector of length one), Default: '1970-01-01'
#' @return Updated dataset (a tibble)
#' @rdname add_dates_from_dist
#' @export 
#' @importFrom rlang exec sym
#' @importFrom dplyr mutate case_when n group_by lag ungroup select everything
#' @importFrom lubridate days
#' @keywords internal
add_dates_from_dist <- function (ds_tb, bl_start_date_dtm, bl_end_date_dtm, duration_args_ls, 
    duration_fn = rnorm, date_var_nm_1L_chr = "date_psx", id_var_nm_1L_chr = "fkClientID", 
    round_var_nm_1L_chr = "round", round_bl_val_1L_chr = "Baseline", 
    origin_1L_chr = "1970-01-01") 
{
    args_ls <- append(list(n = nrow(ds_tb)), duration_args_ls)
    days_of_fup_int <- rlang::exec(.fn = duration_fn, !!!args_ls) %>% 
        round(0) %>% as.integer()
    updated_ds_tb <- ds_tb %>% dplyr::mutate(duration_prd = dplyr::case_when(!!rlang::sym(round_var_nm_1L_chr) != 
        round_bl_val_1L_chr ~ lubridate::days(days_of_fup_int), 
        T ~ lubridate::days(0))) %>% dplyr::mutate(`:=`(!!rlang::sym(date_var_nm_1L_chr), 
        dplyr::case_when(!!rlang::sym(round_var_nm_1L_chr) == 
            round_bl_val_1L_chr ~ as.Date(sample(as.numeric(bl_start_date_dtm):as.numeric(bl_end_date_dtm), 
            dplyr::n(), replace = T), origin = origin_1L_chr)))) %>% 
        dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) %>% dplyr::mutate(`:=`(!!rlang::sym(date_var_nm_1L_chr), 
        dplyr::case_when(!!rlang::sym(round_var_nm_1L_chr) == 
            round_bl_val_1L_chr ~ !!rlang::sym(date_var_nm_1L_chr), 
            T ~ dplyr::lag(!!rlang::sym(date_var_nm_1L_chr)) + 
                duration_prd))) %>% dplyr::ungroup() %>% dplyr::select(!!rlang::sym(id_var_nm_1L_chr), 
        !!rlang::sym(round_var_nm_1L_chr), !!rlang::sym(date_var_nm_1L_chr), 
        duration_prd, dplyr::everything())
    return(updated_ds_tb)
}
#' Add diffs by group and tmpt
#' @description add_diffs_by_group_and_tmpt() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add diffs by group and tmpt. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble), Default: trial_ds_tb
#' @param cmprsn_var_nm_1L_chr Cmprsn var name (a character vector of length one), Default: 'study_arm_chr'
#' @param cmprsn_group_match_val_chr Cmprsn group match value (a character vector), Default: c("Intervention")
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param timepoint_match_val_1L_chr Timepoint match value (a character vector of length one), Default: 'Follow-up'
#' @param match_idx_var_nm_1L_chr Match index var name (a character vector of length one), Default: 'match_idx_int'
#' @param var_nms_chr Var names (a character vector)
#' @param fns_ls Functions (a list)
#' @param abs_mean_diff_dbl Abs mean diff (a double vector)
#' @param diff_sd_dbl Diff sd (a double vector)
#' @param multiplier_dbl Multiplier (a double vector)
#' @param min_dbl Min (a double vector)
#' @param max_dbl Max (a double vector)
#' @param integer_lgl Integer (a logical vector)
#' @return Updated dataset (a tibble)
#' @rdname add_diffs_by_group_and_tmpt
#' @export 
#' @importFrom dplyr filter bind_rows arrange
#' @importFrom rlang sym
#' @keywords internal
add_diffs_by_group_and_tmpt <- function (ds_tb = trial_ds_tb, cmprsn_var_nm_1L_chr = "study_arm_chr", 
    cmprsn_group_match_val_chr = c("Intervention"), round_var_nm_1L_chr = "round", 
    timepoint_match_val_1L_chr = "Follow-up", match_idx_var_nm_1L_chr = "match_idx_int", 
    var_nms_chr, fns_ls, abs_mean_diff_dbl, diff_sd_dbl, multiplier_dbl, 
    min_dbl, max_dbl, integer_lgl) 
{
    unchanged_vals_tb <- ds_tb %>% dplyr::filter(!(!!rlang::sym(cmprsn_var_nm_1L_chr) == 
        cmprsn_group_match_val_chr & !!rlang::sym(round_var_nm_1L_chr) == 
        timepoint_match_val_1L_chr))
    updated_vals_tb <- ds_tb %>% dplyr::filter(!!rlang::sym(cmprsn_var_nm_1L_chr) == 
        cmprsn_group_match_val_chr & !!rlang::sym(round_var_nm_1L_chr) == 
        timepoint_match_val_1L_chr) %>% update_multpl_cols_with_diffs(var_nms_chr = var_nms_chr, 
        fns_ls = fns_ls, abs_mean_diff_dbl = abs_mean_diff_dbl, 
        diff_sd_dbl = diff_sd_dbl, multiplier_dbl = multiplier_dbl, 
        min_dbl = min_dbl, max_dbl = max_dbl, integer_lgl = integer_lgl)
    updated_ds_tb <- dplyr::bind_rows(unchanged_vals_tb, updated_vals_tb) %>% 
        dplyr::arrange(!!rlang::sym(match_idx_var_nm_1L_chr))
    return(updated_ds_tb)
}
#' Add qalys
#' @description add_qalys() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add qalys. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param cmprsn_var_nm_1L_chr Cmprsn var name (a character vector of length one), Default: 'study_arm_chr'
#' @param duration_var_nm_1L_chr Duration var name (a character vector of length one), Default: 'duration_prd'
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param match_idx_var_nm_1L_chr Match index var name (a character vector of length one), Default: 'match_idx_int'
#' @param qalys_var_nm_1L_chr Qalys var name (a character vector of length one), Default: 'qalys_dbl'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param utl_change_var_nm_1L_chr Utl change var name (a character vector of length one), Default: 'utl_change_dbl'
#' @param utl_var_nm_1L_chr Utl var name (a character vector of length one), Default: 'utility_dbl'
#' @param reshape_1L_lgl Reshape (a logical vector of length one), Default: T
#' @return Updated dataset (a tibble)
#' @rdname add_qalys
#' @export 
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @importFrom lubridate years
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect all_of
#' @keywords internal
add_qalys <- function (ds_tb, cmprsn_var_nm_1L_chr = "study_arm_chr", duration_var_nm_1L_chr = "duration_prd", 
    id_var_nm_1L_chr = "fkClientID", match_idx_var_nm_1L_chr = "match_idx_int", 
    qalys_var_nm_1L_chr = "qalys_dbl", round_var_nm_1L_chr = "round", 
    utl_change_var_nm_1L_chr = "utl_change_dbl", utl_var_nm_1L_chr = "utility_dbl", 
    reshape_1L_lgl = T) 
{
    updated_ds_tb <- ds_tb %>% dplyr::mutate(`:=`(!!rlang::sym(qalys_var_nm_1L_chr), 
        (!!rlang::sym(utl_var_nm_1L_chr) - (!!rlang::sym(utl_change_var_nm_1L_chr) * 
            0.5)) * (duration_prd/lubridate::years(1))))
    if (reshape_1L_lgl) {
        vars_to_spread_chr <- names(updated_ds_tb)[!names(updated_ds_tb) %in% 
            c(cmprsn_var_nm_1L_chr, id_var_nm_1L_chr, match_idx_var_nm_1L_chr, 
                round_var_nm_1L_chr)]
        updated_ds_tb <- updated_ds_tb %>% tidyr::pivot_wider(names_from = !!rlang::sym(round_var_nm_1L_chr), 
            values_from = tidyselect::all_of(vars_to_spread_chr))
    }
    return(updated_ds_tb)
}
#' Add qalys to dataset
#' @description add_qalys_to_ds() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add qalys to dataset. Function argument ds_tb specifies the object to be updated. The function returns Dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param ds_smry_ls Dataset smry (a list)
#' @return Dataset (a tibble)
#' @rdname add_qalys_to_ds
#' @export 
#' @importFrom purrr map reduce
add_qalys_to_ds <- function (ds_tb, ds_smry_ls) 
{
    args_ls_ls <- purrr::map(c(ds_smry_ls$predr_var_nms, ds_smry_ls$utl_var_nm_1L_chr), 
        ~list(change_var_nm_1L_chr = paste0(.x, "_change_dbl"), 
            var_nm_1L_chr = .x))
    ds_tb <- purrr::reduce(1:length(args_ls_ls), .init = ds_tb, 
        ~add_change_in_ds_var(.x, var_nm_1L_chr = args_ls_ls[[.y]]$var_nm_1L_chr, 
            change_var_nm_1L_chr = args_ls_ls[[.y]]$change_var_nm_1L_chr)) %>% 
        add_qalys(utl_change_var_nm_1L_chr = paste0(ds_smry_ls$utl_var_nm_1L_chr, 
            "_change_dbl"), utl_var_nm_1L_chr = ds_smry_ls$utl_var_nm_1L_chr, 
            duration_var_nm_1L_chr = "duration_prd", qalys_var_nm_1L_chr = "qalys_dbl")
    return(ds_tb)
}
