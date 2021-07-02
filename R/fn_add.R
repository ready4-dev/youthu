#' Add change in dataset variable
#' @description add_change_in_ds_var() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add change in dataset variable. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round baseline value (a character vector of length one), Default: 'Baseline'
#' @param change_var_nm_1L_chr Change variable name (a character vector of length one)
#' @param var_nm_1L_chr Variable name (a character vector of length one)
#' @param arrange_by_id_lgl Arrange by identity (a logical vector), Default: T
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
#' Add costs by time point
#' @description add_costs_by_tmpt() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add costs by time point. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one)
#' @param round_lvls_chr Round levels (a character vector), Default: c("Baseline", "Follow-up")
#' @param costs_mean_dbl Costs mean (a double vector)
#' @param costs_sd_dbl Costs standard deviation (a double vector)
#' @param extra_cost_args_ls Extra cost arguments (a list), Default: list(costs_var_nm_1L_chr = "costs_dbl")
#' @param fn Function (a function), Default: add_costs_from_gamma_dstr
#' @return Updated dataset (a tibble)
#' @rdname add_costs_by_tmpt
#' @export 
#' @importFrom purrr pmap_dfr
#' @importFrom dplyr filter
#' @importFrom rlang sym exec
#' @keywords internal
add_costs_by_tmpt <- function (ds_tb, round_var_nm_1L_chr, round_lvls_chr = c("Baseline", 
    "Follow-up"), costs_mean_dbl, costs_sd_dbl, extra_cost_args_ls = list(costs_var_nm_1L_chr = "costs_dbl"), 
    fn = add_costs_from_gamma_dstr) 
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
#' Add costs from gamma distribution
#' @description add_costs_from_gamma_dstr() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add costs from gamma distribution. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param costs_mean_dbl Costs mean (a double vector)
#' @param costs_sd_dbl Costs standard deviation (a double vector)
#' @param costs_var_nm_1L_chr Costs variable name (a character vector of length one), Default: 'costs_dbl'
#' @return Updated dataset (a tibble)
#' @rdname add_costs_from_gamma_dstr
#' @export 
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @keywords internal
add_costs_from_gamma_dstr <- function (ds_tb, costs_mean_dbl, costs_sd_dbl, costs_var_nm_1L_chr = "costs_dbl") 
{
    updated_ds_tb <- dplyr::mutate(ds_tb, `:=`(!!rlang::sym(costs_var_nm_1L_chr), 
        make_costs_vec_from_gamma_dstr(n_int = nrow(ds_tb), costs_mean_dbl = costs_mean_dbl, 
            costs_sd_dbl = costs_sd_dbl)))
    return(updated_ds_tb)
}
#' Add dates from distribution
#' @description add_dates_from_dstr() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add dates from distribution. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param bl_start_date_dtm Baseline start date (a date vector)
#' @param bl_end_date_dtm Baseline end date (a date vector)
#' @param duration_args_ls Duration arguments (a list)
#' @param duration_fn Duration (a function), Default: stats::rnorm
#' @param date_var_nm_1L_chr Date variable name (a character vector of length one), Default: 'date_psx'
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'round'
#' @param round_bl_val_1L_chr Round baseline value (a character vector of length one), Default: 'Baseline'
#' @param origin_1L_chr Origin (a character vector of length one), Default: '1970-01-01'
#' @return Updated dataset (a tibble)
#' @rdname add_dates_from_dstr
#' @export 
#' @importFrom stats rnorm
#' @importFrom rlang exec sym
#' @importFrom dplyr mutate case_when n group_by lag ungroup select everything
#' @importFrom lubridate days
#' @keywords internal
add_dates_from_dstr <- function (ds_tb, bl_start_date_dtm, bl_end_date_dtm, duration_args_ls, 
    duration_fn = stats::rnorm, date_var_nm_1L_chr = "date_psx", 
    id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round", 
    round_bl_val_1L_chr = "Baseline", origin_1L_chr = "1970-01-01") 
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
#' Add differences by group and time point
#' @description add_diffs_by_group_and_tmpt() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add differences by group and time point. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble), Default: trial_ds_tb
#' @param cmprsn_var_nm_1L_chr Comparison variable name (a character vector of length one), Default: 'study_arm_chr'
#' @param cmprsn_group_match_val_chr Comparison group match value (a character vector), Default: c("Intervention")
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'round'
#' @param timepoint_match_val_1L_chr Timepoint match value (a character vector of length one), Default: 'Follow-up'
#' @param match_idx_var_nm_1L_chr Match index variable name (a character vector of length one), Default: 'match_idx_int'
#' @param var_nms_chr Variable names (a character vector)
#' @param fns_ls Functions (a list)
#' @param abs_mean_diff_dbl Absolute mean difference (a double vector)
#' @param diff_sd_dbl Difference standard deviation (a double vector)
#' @param multiplier_dbl Multiplier (a double vector)
#' @param min_dbl Minimum (a double vector)
#' @param max_dbl Maximum (a double vector)
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
#' Add Quality Adjusted Life Years
#' @description add_qalys() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add quality adjusted life years. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param cmprsn_var_nm_1L_chr Comparison variable name (a character vector of length one), Default: 'study_arm_chr'
#' @param duration_var_nm_1L_chr Duration variable name (a character vector of length one), Default: 'duration_prd'
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param match_idx_var_nm_1L_chr Match index variable name (a character vector of length one), Default: 'match_idx_int'
#' @param qalys_var_nm_1L_chr Quality Adjusted Life Years variable name (a character vector of length one), Default: 'qalys_dbl'
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'round'
#' @param utl_change_var_nm_1L_chr Utility change variable name (a character vector of length one), Default: 'utl_change_dbl'
#' @param utl_var_nm_1L_chr Utility variable name (a character vector of length one), Default: 'utility_dbl'
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
            0.5)) * (!!rlang::sym(duration_var_nm_1L_chr)/lubridate::years(1))))
    if (reshape_1L_lgl) {
        vars_to_spread_chr <- names(updated_ds_tb)[!names(updated_ds_tb) %in% 
            c(cmprsn_var_nm_1L_chr, id_var_nm_1L_chr, match_idx_var_nm_1L_chr, 
                round_var_nm_1L_chr)]
        updated_ds_tb <- updated_ds_tb %>% tidyr::pivot_wider(names_from = !!rlang::sym(round_var_nm_1L_chr), 
            values_from = tidyselect::all_of(vars_to_spread_chr))
    }
    return(updated_ds_tb)
}
#' Add Quality Adjusted Life Years to dataset
#' @description add_qalys_to_ds() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add quality adjusted life years to dataset. Function argument ds_tb specifies the object to be updated. The function returns Dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param predn_ds_ls Prediction dataset (a list)
#' @return Dataset (a tibble)
#' @rdname add_qalys_to_ds
#' @export 
#' @importFrom purrr map reduce
add_qalys_to_ds <- function (ds_tb, predn_ds_ls) 
{
    if (is.null(predn_ds_ls$ds_ls$predr_vars_nms_chr)) 
        predn_ds_ls$ds_ls$predr_vars_nms_chr <- predn_ds_ls$mdl_ls$predictors_lup$short_name_chr
    ds_smry_ls <- predn_ds_ls$ds_ls
    args_ls_ls <- purrr::map(c(ds_smry_ls$predr_vars_nms_chr, 
        ds_smry_ls$utl_var_nm_1L_chr), ~list(change_var_nm_1L_chr = paste0(.x, 
        "_change_dbl"), var_nm_1L_chr = .x))
    ds_tb <- purrr::reduce(1:length(args_ls_ls), .init = ds_tb, 
        ~add_change_in_ds_var(.x, var_nm_1L_chr = args_ls_ls[[.y]]$var_nm_1L_chr, 
            change_var_nm_1L_chr = args_ls_ls[[.y]]$change_var_nm_1L_chr)) %>% 
        add_qalys(utl_change_var_nm_1L_chr = paste0(ds_smry_ls$utl_var_nm_1L_chr, 
            "_change_dbl"), utl_var_nm_1L_chr = ds_smry_ls$utl_var_nm_1L_chr, 
            duration_var_nm_1L_chr = "duration_prd", qalys_var_nm_1L_chr = "qalys_dbl")
    return(ds_tb)
}
#' Add utility prediction
#' @description add_utl_predn() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add utility prediction. Function argument data_tb specifies the object to be updated. The function returns Updated (a tibble).
#' @param data_tb Data (a tibble)
#' @param predn_ds_ls Prediction dataset (a list)
#' @param deterministic_1L_lgl Deterministic (a logical vector of length one), Default: T
#' @param force_min_max_1L_lgl Force minimum maximum (a logical vector of length one), Default: T
#' @param key_1L_chr Key (a character vector of length one), Default: NULL
#' @param make_from_tbl_1L_lgl Make from table (a logical vector of length one), Default: T
#' @param model_mdl Model (a model), Default: NULL
#' @param new_data_is_1L_chr New data is (a character vector of length one), Default: 'Simulated'
#' @param server_1L_chr Server (a character vector of length one), Default: 'dataverse.harvard.edu'
#' @param utl_cls_fn Utility class (a function), Default: NULL
#' @return Updated (a tibble)
#' @rdname add_utl_predn
#' @export 
#' @importFrom TTU add_utl_predn_to_new_ds
#' @importFrom ready4fun get_from_lup_obj
add_utl_predn <- function (data_tb, predn_ds_ls, deterministic_1L_lgl = T, force_min_max_1L_lgl = T, 
    key_1L_chr = NULL, make_from_tbl_1L_lgl = T, model_mdl = NULL, 
    new_data_is_1L_chr = "Simulated", server_1L_chr = "dataverse.harvard.edu", 
    utl_cls_fn = NULL) 
{
    id_var_nm_1L_chr = predn_ds_ls$ds_ls$id_var_nm_1L_chr
    predr_vars_nms_chr = predn_ds_ls$ds_ls$predr_vars_nms_chr
    round_var_nm_1L_chr = predn_ds_ls$ds_ls$round_var_nm_1L_chr
    round_bl_val_1L_chr = predn_ds_ls$ds_ls$round_bl_val_1L_chr
    utl_var_nm_1L_chr = predn_ds_ls$ds_ls$utl_var_nm_1L_chr
    mdl_meta_data_ls = predn_ds_ls$mdl_ls$mdl_meta_data_ls
    mdls_lup = predn_ds_ls$mdl_ls$mdls_lup
    mdl_nm_1L_chr = predn_ds_ls$mdl_ls$mdl_nm_1L_chr
    if (is.null(model_mdl)) 
        model_mdl <- get_model(mdls_lup, mdl_nm_1L_chr = mdl_nm_1L_chr, 
            make_from_tbl_1L_lgl = make_from_tbl_1L_lgl, mdl_meta_data_ls = mdl_meta_data_ls, 
            server_1L_chr = server_1L_chr, key_1L_chr = key_1L_chr)
    updated_tb <- TTU::add_utl_predn_to_new_ds(data_tb, mdl_nm_1L_chr = mdl_nm_1L_chr, 
        id_var_nm_1L_chr = id_var_nm_1L_chr, analysis_1L_chr = ready4fun::get_from_lup_obj(mdls_lup, 
            match_value_xx = mdl_nm_1L_chr, match_var_nm_1L_chr = "mdl_nms_chr", 
            target_var_nm_1L_chr = "source_chr", evaluate_lgl = F), 
        ingredients_ls = get_mdl_metadata(mdls_lup, mdl_nm_1L_chr = mdl_nm_1L_chr), 
        model_mdl = model_mdl, new_data_is_1L_chr = new_data_is_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, round_var_nm_1L_chr = round_var_nm_1L_chr, 
        round_bl_val_1L_chr = round_bl_val_1L_chr, utl_var_nm_1L_chr = utl_var_nm_1L_chr)
    return(updated_tb)
}
