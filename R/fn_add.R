#' Add adol6d scores
#' @description add_adol6d_scores() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add adol6d scores. Function argument unscored_aqol_tb specifies the object to be updated. The function returns Transformed Assessment of Quality of Life (a tibble).
#' @param unscored_aqol_tb Unscored Assessment of Quality of Life (a tibble)
#' @param prefix_1L_chr Prefix (a character vector of length one), Default: 'aqol6d_q'
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param wtd_aqol_var_nm_1L_chr Wtd Assessment of Quality of Life var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @return Transformed Assessment of Quality of Life (a tibble)
#' @rdname add_adol6d_scores
#' @export 
#' @importFrom dplyr select starts_with rename_all inner_join rename
#' @importFrom stringr str_replace
#' @importFrom rlang sym
#' @keywords internal
add_adol6d_scores <- function (unscored_aqol_tb, prefix_1L_chr = "aqol6d_q", id_var_nm_1L_chr = "fkClientID", 
    wtd_aqol_var_nm_1L_chr = "aqol6d_total_w") 
{
    complete_ds_tb <- unscored_aqol_tb
    unscored_aqol_tb <- unscored_aqol_tb %>% dplyr::select(id_var_nm_1L_chr, 
        dplyr::starts_with(unname(prefix_1L_chr)))
    old_nms_chr <- names(unscored_aqol_tb)
    names(unscored_aqol_tb) <- c("ID", paste0("Q", 1:20))
    unscored_aqol_tb <- impute_unscrd_adol_aqol6d_ds(unscored_aqol_tb)
    disvals_tb <- unscored_aqol_tb %>% add_itm_disv_to_aqol6d_itms_tb(disvalues_lup_tb = make_adol_aqol6d_disv_lup(), 
        pfx_1L_chr = "Q") %>% dplyr::select(ID, dplyr::starts_with("dv_")) %>% 
        dplyr::rename_all(~stringr::str_replace(.x, "dv_", "dv"))
    scored_aqol_tb <- add_aqol6d_adol_dim_scrg_eqs(disvals_tb)
    tfd_aqol_tb <- dplyr::inner_join(complete_ds_tb, scored_aqol_tb %>% 
        dplyr::rename(`:=`(!!rlang::sym(id_var_nm_1L_chr), ID), 
            `:=`(!!rlang::sym(wtd_aqol_var_nm_1L_chr), uaqol)))
    return(tfd_aqol_tb)
}
#' Add Assessment of Quality of Life Six Dimension adolescent dimension scoring equations
#' @description add_aqol6d_adol_dim_scrg_eqs() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add assessment of quality of life six dimension adolescent dimension scoring equations. Function argument unscored_aqol_tb specifies the object to be updated. The function returns Unscored Assessment of Quality of Life (a tibble).
#' @param unscored_aqol_tb Unscored Assessment of Quality of Life (a tibble)
#' @return Unscored Assessment of Quality of Life (a tibble)
#' @rdname add_aqol6d_adol_dim_scrg_eqs
#' @export 
#' @importFrom utils data
#' @importFrom dplyr mutate
#' @importFrom rlang parse_expr
#' @importFrom Hmisc label
#' @keywords internal
add_aqol6d_adol_dim_scrg_eqs <- function (unscored_aqol_tb) 
{
    utils::data("adol_dim_scalg_eqs_lup", envir = environment())
    for (var in adol_dim_scalg_eqs_lup$Dim_scal) {
        expression = adol_dim_scalg_eqs_lup[adol_dim_scalg_eqs_lup$Dim_scal == 
            var, ]$Equ
        unscored_aqol_tb <- unscored_aqol_tb %>% dplyr::mutate(`:=`(!!var, 
            !!rlang::parse_expr(expression)))
        Hmisc::label(unscored_aqol_tb[, var]) = adol_dim_scalg_eqs_lup[adol_dim_scalg_eqs_lup$Dim_scal == 
            var, ]$Label
    }
    return(unscored_aqol_tb)
}
#' Add Assessment of Quality of Life Six Dimension items to Assessment of Quality of Life Six Dimension tibbles
#' @description add_aqol6d_items_to_aqol6d_tbs_ls() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add assessment of quality of life six dimension items to assessment of quality of life six dimension tibbles list. Function argument aqol6d_tbs_ls specifies the object to be updated. The function returns Updated Assessment of Quality of Life Six Dimension tibbles (a list).
#' @param aqol6d_tbs_ls Assessment of Quality of Life Six Dimension tibbles (a list)
#' @param aqol_items_props_tbs_ls Assessment of Quality of Life items props tibbles (a list)
#' @param prefix_chr Prefix (a character vector)
#' @param aqol_tots_var_nms_chr Assessment of Quality of Life totals var names (a character vector)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param scaling_cnst_dbl Scaling cnst (a double vector), Default: 5
#' @return Updated Assessment of Quality of Life Six Dimension tibbles (a list)
#' @rdname add_aqol6d_items_to_aqol6d_tbs_ls
#' @export 
#' @importFrom purrr map2 map reduce map_int
#' @importFrom dplyr select mutate arrange left_join starts_with everything
#' @importFrom stats na.omit
#' @importFrom simstudy defData genData
#' @importFrom rlang sym
#' @importFrom tibble rowid_to_column
#' @keywords internal
add_aqol6d_items_to_aqol6d_tbs_ls <- function (aqol6d_tbs_ls, aqol_items_props_tbs_ls, prefix_chr, 
    aqol_tots_var_nms_chr, id_var_nm_1L_chr = "fkClientID", scaling_cnst_dbl = 5) 
{
    updated_aqol6d_tbs_ls <- purrr::map2(aqol6d_tbs_ls, aqol_items_props_tbs_ls, 
        ~{
            nbr_obs_1L_int <- nrow(.x) * scaling_cnst_dbl
            transposed_items_props_tb <- .y %>% dplyr::select(-Question) %>% 
                t()
            item_ranges_dbl_ls <- 1:ncol(transposed_items_props_tb) %>% 
                purrr::map(~c(1, length(transposed_items_props_tb[, 
                  .x] %>% stats::na.omit())))
            cat_probs_def_tbl <- purrr::reduce(1:ncol(transposed_items_props_tb), 
                .init = NULL, ~simstudy::defData(.x, varname = paste0("aqol6d_q", 
                  .y), formula = transposed_items_props_tb[, 
                  .y] %>% stats::na.omit() %>% as.vector() %>% 
                  format(digits = 10) %>% paste0(collapse = ";"), 
                  dist = "categorical"))
            items_tb <- simstudy::genData(nbr_obs_1L_int, cat_probs_def_tbl) %>% 
                dplyr::select(-id) %>% dplyr::mutate(`:=`(!!rlang::sym(unname(aqol_tots_var_nms_chr["cumulative"])), 
                rowSums(., na.rm = T))) %>% dplyr::arrange(!!rlang::sym(unname(aqol_tots_var_nms_chr["cumulative"]))) %>% 
                tibble::rowid_to_column("id")
            items_tb <- items_tb %>% dplyr::mutate(aqol6dU = calculate_adol_aqol6dU(items_tb, 
                prefix_1L_chr = prefix_chr["aqol_item"], id_var_nm_1L_chr = "id"))
            .x <- .x %>% dplyr::mutate(id = purrr::map_int(aqol6d_total_w, 
                ~which.min(abs(items_tb$aqol6dU - .x)))) %>% 
                dplyr::left_join(items_tb)
            updated_tb <- .x %>% dplyr::mutate(`:=`(!!rlang::sym(unname(aqol_tots_var_nms_chr["weighted"])), 
                aqol6dU)) %>% dplyr::select(-aqol6dU, -id) %>% 
                dplyr::select(!!rlang::sym(id_var_nm_1L_chr), 
                  dplyr::starts_with(prefix_chr[["aqol_item"]]), 
                  !!rlang::sym(unname(aqol_tots_var_nms_chr["cumulative"])), 
                  !!rlang::sym(unname(aqol_tots_var_nms_chr["weighted"])), 
                  dplyr::everything())
            updated_tb
        })
    return(updated_aqol6d_tbs_ls)
}
#' Add Assessment of Quality of Life Six Dimension predn to dataset
#' @description add_aqol6d_predn_to_ds() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add assessment of quality of life six dimension predn to dataset. Function argument data_tb specifies the object to be updated. The function returns Updated (a tibble).
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
#' @keywords internal
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
    if (!is.null(names(predr_vars_nms_chr))) {
        updated_tb <- rename_from_nmd_vec(updated_tb, nmd_vec_chr = predr_vars_nms_chr, 
            vec_nms_as_new_1L_lgl = F)
    }
    if (!is.null(utl_var_nm_1L_chr)) {
        updated_tb <- updated_tb %>% dplyr::rename(`:=`(!!rlang::sym(utl_var_nm_1L_chr), 
            tidyselect::all_of(mdl_dep_var_1L_chr)))
    }
    if ("aqol6d_total_w_CLL_cloglog" %in% names(updated_tb)) 
        updated_tb <- updated_tb %>% dplyr::select(-aqol6d_total_w_CLL_cloglog)
    updated_tb <- dplyr::left_join(data_tb %>% dplyr::select(original_ds_vars_chr), 
        updated_tb) %>% dplyr::select(c(names(updated_tb), setdiff(names(data_tb), 
        names(updated_tb))))
    return(updated_tb)
}
#' Add Assessment of Quality of Life Six Dimension Health Utility to Assessment of Quality of Life Six Dimension items
#' @description add_aqol6dU_to_aqol6d_items_tb() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add assessment of quality of life six dimension health utility to assessment of quality of life six dimension items tibble. Function argument aqol6d_items_tb specifies the object to be updated. The function returns Assessment of Quality of Life Six Dimension items (a tibble).
#' @param aqol6d_items_tb Assessment of Quality of Life Six Dimension items (a tibble)
#' @param coeffs_lup_tb Coeffs lookup table (a tibble), Default: aqol6d_from_8d_coeffs_lup_tb
#' @return Assessment of Quality of Life Six Dimension items (a tibble)
#' @rdname add_aqol6dU_to_aqol6d_items_tb
#' @export 
#' @importFrom dplyr pull mutate
#' @importFrom purrr map_dbl
add_aqol6dU_to_aqol6d_items_tb <- function (aqol6d_items_tb, coeffs_lup_tb = aqol6d_from_8d_coeffs_lup_tb) 
{
    coeff_dbl <- coeffs_lup_tb[match(c(paste0("vD", 1:6), "Constant"), 
        coeffs_lup_tb$var_name_chr), ] %>% dplyr::pull(coeff_dbl)
    aqol6d_items_tb <- aqol6d_items_tb %>% dplyr::mutate(aqol6dU = coeff_dbl[1] * 
        vD1 + coeff_dbl[2] * vD2 + coeff_dbl[3] * vD3 + coeff_dbl[4] * 
        vD4 + coeff_dbl[5] * vD5 + coeff_dbl[6] * vD6 + coeff_dbl[7]) %>% 
        dplyr::mutate(aqol6dU = aqol6dU %>% purrr::map_dbl(~ifelse(.x > 
            1, 1, .x)))
    return(aqol6d_items_tb)
}
#' Add Assessment of Quality of Life Six Dimension Health Utility to Assessment of Quality of Life Six Dimension tibbles
#' @description add_aqol6dU_to_aqol6d_tbs_ls() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add assessment of quality of life six dimension health utility to assessment of quality of life six dimension tibbles list. Function argument aqol6d_tbs_ls specifies the object to be updated. The function returns Assessment of Quality of Life Six Dimension tibbles (a list).
#' @param aqol6d_tbs_ls Assessment of Quality of Life Six Dimension tibbles (a list)
#' @param prefix_1L_chr Prefix (a character vector of length one), Default: 'aqol6d_q'
#' @param id_var_nm_1L_chr Id var name (a character vector of length one)
#' @return Assessment of Quality of Life Six Dimension tibbles (a list)
#' @rdname add_aqol6dU_to_aqol6d_tbs_ls
#' @export 
#' @importFrom purrr map
#' @importFrom dplyr mutate
#' @keywords internal
add_aqol6dU_to_aqol6d_tbs_ls <- function (aqol6d_tbs_ls, prefix_1L_chr = "aqol6d_q", id_var_nm_1L_chr) 
{
    aqol6d_tbs_ls <- aqol6d_tbs_ls %>% purrr::map(~.x %>% dplyr::mutate(aqol6dU = calculate_adol_aqol6dU(.x, 
        prefix_1L_chr = prefix_1L_chr, id_var_nm_1L_chr = id_var_nm_1L_chr)))
    return(aqol6d_tbs_ls)
}
#' Add cors and uts to Assessment of Quality of Life Six Dimension tibbles
#' @description add_cors_and_uts_to_aqol6d_tbs_ls() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add cors and uts to assessment of quality of life six dimension tibbles list. Function argument aqol6d_tbs_ls specifies the object to be updated. The function returns Assessment of Quality of Life Six Dimension tibbles (a list).
#' @param aqol6d_tbs_ls Assessment of Quality of Life Six Dimension tibbles (a list)
#' @param aqol_scores_pars_ls Assessment of Quality of Life scores parameters (a list)
#' @param aqol_items_props_tbs_ls Assessment of Quality of Life items props tibbles (a list)
#' @param temporal_cors_ls Temporal cors (a list)
#' @param prefix_chr Prefix (a character vector)
#' @param aqol_tots_var_nms_chr Assessment of Quality of Life totals var names (a character vector)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @return Assessment of Quality of Life Six Dimension tibbles (a list)
#' @rdname add_cors_and_uts_to_aqol6d_tbs_ls
#' @export 

#' @keywords internal
add_cors_and_uts_to_aqol6d_tbs_ls <- function (aqol6d_tbs_ls, aqol_scores_pars_ls, aqol_items_props_tbs_ls, 
    temporal_cors_ls, prefix_chr, aqol_tots_var_nms_chr, id_var_nm_1L_chr = "fkClientID") 
{
    aqol6d_tbs_ls <- reorder_tbs_for_target_cors(aqol6d_tbs_ls, 
        cor_dbl = temporal_cors_ls[[1]], cor_var_chr = rep(names(temporal_cors_ls)[1], 
            2), id_var_to_rm_1L_chr = "id") %>% add_uids_to_tbs_ls(prefix_1L_chr = prefix_chr[["uid"]], 
        id_var_nm_1L_chr = id_var_nm_1L_chr)
    aqol6d_tbs_ls <- aqol6d_tbs_ls %>% add_aqol6d_items_to_aqol6d_tbs_ls(aqol_items_props_tbs_ls = aqol_items_props_tbs_ls, 
        prefix_chr = prefix_chr, aqol_tots_var_nms_chr = aqol_tots_var_nms_chr, 
        id_var_nm_1L_chr = id_var_nm_1L_chr)
    return(aqol6d_tbs_ls)
}
#' Add costs by tmpt
#' @description add_costs_by_tmpt() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add costs by tmpt. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param round_var_nm_1L_chr Round var name (a character vector of length one)
#' @param round_lvls_chr Round lvls (a character vector), Default: c("Baseline", "Follow-up")
#' @param costs_mean_dbl Costs mean (a double vector)
#' @param costs_sd_dbl Costs sd (a double vector)
#' @param fn Function (a function), Default: add_costs_from_gamma_dist
#' @return Updated dataset (a tibble)
#' @rdname add_costs_by_tmpt
#' @export 
#' @importFrom purrr pmap_dfr
#' @importFrom dplyr filter
#' @importFrom rlang sym exec
#' @keywords internal
add_costs_by_tmpt <- function (ds_tb, round_var_nm_1L_chr, round_lvls_chr = c("Baseline", 
    "Follow-up"), costs_mean_dbl, costs_sd_dbl, fn = add_costs_from_gamma_dist) 
{
    updated_ds_tb <- purrr::pmap_dfr(list(round_lvls_chr, costs_mean_dbl, 
        costs_sd_dbl), ~{
        args_ls <- list(ds_tb %>% dplyr::filter(!!rlang::sym(round_var_nm_1L_chr) == 
            ..1), ..2, ..3)
        rlang::exec(.fn = fn, !!!args_ls)
    })
    return(updated_ds_tb)
}
#' Add costs from gamma dist
#' @description add_costs_from_gamma_dist() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add costs from gamma dist. Function argument ds_tb specifies the object to be updated. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param costs_mean_dbl Costs mean (a double vector)
#' @param costs_sd_dbl Costs sd (a double vector)
#' @param costs_var_nm_1L_dbl Costs var name (a double vector of length one), Default: 'costs_dbl'
#' @return Updated dataset (a tibble)
#' @rdname add_costs_from_gamma_dist
#' @export 
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @keywords internal
add_costs_from_gamma_dist <- function (ds_tb, costs_mean_dbl, costs_sd_dbl, costs_var_nm_1L_dbl = "costs_dbl") 
{
    updated_ds_tb <- dplyr::mutate(ds_tb, `:=`(!!rlang::sym(costs_var_nm_1L_dbl), 
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
#' Add dimension disvalue to Assessment of Quality of Life Six Dimension items
#' @description add_dim_disv_to_aqol6d_items_tb() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add dimension disvalue to assessment of quality of life six dimension items tibble. Function argument aqol6d_items_tb specifies the object to be updated. The function returns Assessment of Quality of Life Six Dimension items (a tibble).
#' @param aqol6d_items_tb Assessment of Quality of Life Six Dimension items (a tibble)
#' @param domain_items_ls Domain items (a list)
#' @param domains_chr Domains (a character vector)
#' @param dim_sclg_con_lup_tb Dimension sclg con lookup table (a tibble), Default: aqol6d_dim_sclg_con_lup_tb
#' @param itm_wrst_wghts_lup_tb Itm wrst wghts lookup table (a tibble), Default: aqol6d_adult_itm_wrst_wghts_lup_tb
#' @return Assessment of Quality of Life Six Dimension items (a tibble)
#' @rdname add_dim_disv_to_aqol6d_items_tb
#' @export 
#' @importFrom purrr reduce
#' @importFrom dplyr select mutate
#' @importFrom rlang sym exec
add_dim_disv_to_aqol6d_items_tb <- function (aqol6d_items_tb, domain_items_ls, domains_chr, dim_sclg_con_lup_tb = aqol6d_dim_sclg_con_lup_tb, 
    itm_wrst_wghts_lup_tb = aqol6d_adult_itm_wrst_wghts_lup_tb) 
{
    aqol6d_disu_fn_ls <- make_aqol6d_fns_ls(domain_items_ls)
    kD_dbl <- make_dim_sclg_cons_dbl(domains_chr = domains_chr, 
        dim_sclg_con_lup_tb = dim_sclg_con_lup_tb)
    w_dbl_ls <- make_item_wrst_wghts_ls_ls(domain_items_ls = domain_items_ls, 
        itm_wrst_wghts_lup_tb = itm_wrst_wghts_lup_tb)
    aqol6d_items_tb <- purrr::reduce(1:length(domain_items_ls), 
        .init = aqol6d_items_tb, ~{
            args_ls <- list(dvQs_tb = .x %>% dplyr::select(domain_items_ls[[.y]] %>% 
                paste0("dv_", .)), kD_1L_dbl = kD_dbl[.y], w_dbl = w_dbl_ls[[.y]])
            .x %>% dplyr::mutate(`:=`(!!rlang::sym(paste0("dvD", 
                .y)), rlang::exec(aqol6d_disu_fn_ls[[.y]], !!!args_ls)))
        })
    return(aqol6d_items_tb)
}
#' Add dimension scores to Assessment of Quality of Life Six Dimension items
#' @description add_dim_scores_to_aqol6d_items_tb() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add dimension scores to assessment of quality of life six dimension items tibble. Function argument aqol6d_items_tb specifies the object to be updated. The function returns Assessment of Quality of Life Six Dimension items (a tibble).
#' @param aqol6d_items_tb Assessment of Quality of Life Six Dimension items (a tibble)
#' @param domain_items_ls Domain items (a list)
#' @return Assessment of Quality of Life Six Dimension items (a tibble)
#' @rdname add_dim_scores_to_aqol6d_items_tb
#' @export 
#' @importFrom dplyr mutate across rename_with
#' @importFrom stringr str_replace
add_dim_scores_to_aqol6d_items_tb <- function (aqol6d_items_tb, domain_items_ls) 
{
    aqol6d_items_tb <- aqol6d_items_tb %>% dplyr::mutate(dplyr::across(paste0("dvD", 
        1:length(domain_items_ls)), .fns = list(vD = ~1 - .x), 
        .names = "{fn}_{col}")) %>% dplyr::rename_with(~stringr::str_replace(., 
        "vD_dvD", "vD"))
    return(aqol6d_items_tb)
}
#' Add itm disvalue to Assessment of Quality of Life Six Dimension itms
#' @description add_itm_disv_to_aqol6d_itms_tb() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add itm disvalue to assessment of quality of life six dimension itms tibble. Function argument aqol6d_items_tb specifies the object to be updated. The function returns Assessment of Quality of Life Six Dimension items (a tibble).
#' @param aqol6d_items_tb Assessment of Quality of Life Six Dimension items (a tibble)
#' @param disvalues_lup_tb Disvalues lookup table (a tibble), Default: aqol6d_adult_disv_lup_tb
#' @param pfx_1L_chr Prefix (a character vector of length one)
#' @return Assessment of Quality of Life Six Dimension items (a tibble)
#' @rdname add_itm_disv_to_aqol6d_itms_tb
#' @export 
#' @importFrom purrr reduce
#' @importFrom dplyr mutate across
#' @importFrom tidyselect all_of
add_itm_disv_to_aqol6d_itms_tb <- function (aqol6d_items_tb, disvalues_lup_tb = aqol6d_adult_disv_lup_tb, 
    pfx_1L_chr) 
{
    aqol6d_items_tb <- purrr::reduce(1:20, .init = aqol6d_items_tb, 
        ~{
            q_1L_chr <- paste0(pfx_1L_chr, .y)
            disu_dbl <- disvalues_lup_tb[.y, -1] %>% as.numeric()
            .x %>% dplyr::mutate(dplyr::across(tidyselect::all_of(q_1L_chr), 
                .fns = list(dv = ~disu_dbl[.x]), .names = "{fn}_{col}"))
        })
    return(aqol6d_items_tb)
}
#' Add labels to Assessment of Quality of Life Six Dimension
#' @description add_labels_to_aqol6d_tb() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add labels to assessment of quality of life six dimension tibble. Function argument aqol6d_tb specifies the object to be updated. The function returns Assessment of Quality of Life Six Dimension (a tibble).
#' @param aqol6d_tb Assessment of Quality of Life Six Dimension (a tibble)
#' @param labels_chr Labels (a character vector), Default: 'NA'
#' @return Assessment of Quality of Life Six Dimension (a tibble)
#' @rdname add_labels_to_aqol6d_tb
#' @export 
#' @importFrom Hmisc label
add_labels_to_aqol6d_tb <- function (aqol6d_tb, labels_chr = NA_character_) 
{
    if (is.na(labels_chr)) 
        labels_chr <- c(fkClientID = "Unique client identifier", 
            round = "Data measurement round", d_age = "Age", 
            d_gender = "Gender", d_sexual_ori_s = "Sexual orientation", 
            d_studying_working = "Work and study", c_p_diag_s = " Primary diagnosis", 
            c_clinical_staging_s = "Clinical stage", c_sofas = "SOFAS", 
            s_centre = "Clinic", d_agegroup = "Age group", d_sex_birth_s = "Sex at birth", 
            d_country_bir_s = "Country of birth", d_ATSI = "Aboriginal and Torres Strait Islander", 
            d_english_home = "English spoken at home", d_english_native = "English is native language", 
            d_relation_s = "Relationship status", aqol6d_total_w = "AQoL health utility", 
            phq9_total = "PHQ9", bads_total = "BADS", gad7_total = "GAD7", 
            oasis_total = "OASIS", scared_total = "SCARED", k6_total = "K6", 
            aqol6d_total_c = "AQoL unweighted total", aqol6d_q1 = "Household tasks", 
            aqol6d_q2 = "Getting around", aqol6d_q3 = "Mobility", 
            aqol6d_q4 = "Self care", aqol6d_q5 = "Enjoy close rels", 
            aqol6d_q6 = "Family rels", aqol6d_q7 = "Community involvement", 
            aqol6d_q8 = "Despair", aqol6d_q9 = "Worry", aqol6d_q10 = "Sad", 
            aqol6d_q11 = "Agitated", aqol6d_q12 = "Energy level", 
            aqol6d_q13 = "Control", aqol6d_q14 = "Coping", aqol6d_q15 = "Frequency of pain", 
            aqol6d_q16 = "Degree of pain", aqol6d_q17 = "Pain interference", 
            aqol6d_q18 = "Vision", aqol6d_q19 = "Hearing", aqol6d_q20 = "Communication", 
            aqol6d_subtotal_c_IL = "Unweighted Independent Living", 
            aqol6d_subtotal_c_REL = "Unweighted Relationships", 
            aqol6d_subtotal_c_MH = "Unweighted Mental Health", 
            aqol6d_subtotal_c_COP = "Unweighted Coping", aqol6d_subtotal_c_P = "Unweighted Pain", 
            aqol6d_subtotal_c_SEN = "Unweighted Sense", aqol6d_subtotal_w_IL = "Independent Living", 
            aqol6d_subtotal_w_REL = "Relationships", aqol6d_subtotal_w_MH = "Mental Health", 
            aqol6d_subtotal_w_COP = "Coping", aqol6d_subtotal_w_P = "Pain", 
            aqol6d_subtotal_w_SEN = "Sense")
    Hmisc::label(aqol6d_tb) = as.list(labels_chr[match(names(aqol6d_tb), 
        names(labels_chr))])
    return(aqol6d_tb)
}
#' Add uids to tibbles
#' @description add_uids_to_tbs_ls() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add uids to tibbles list. Function argument tbs_ls specifies the object to be updated. The function returns Tibbles (a list).
#' @param tbs_ls Tibbles (a list)
#' @param prefix_1L_chr Prefix (a character vector of length one)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @return Tibbles (a list)
#' @rdname add_uids_to_tbs_ls
#' @export 
#' @importFrom purrr map map_chr
#' @importFrom dplyr mutate arrange
#' @importFrom rlang sym
#' @importFrom tidyselect all_of
#' @importFrom stringr str_replace
#' @importFrom stats setNames
add_uids_to_tbs_ls <- function (tbs_ls, prefix_1L_chr, id_var_nm_1L_chr = "fkClientID") 
{
    participant_ids <- paste0(prefix_1L_chr, 1:nrow(tbs_ls[[1]])) %>% 
        sample(nrow(tbs_ls[[1]]))
    tbs_ls <- purrr::map(tbs_ls, ~{
        .x %>% dplyr::mutate(`:=`(!!rlang::sym(id_var_nm_1L_chr), 
            tidyselect::all_of(participant_ids[1:nrow(.x)]))) %>% 
            dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr) %>% 
                purrr::map_chr(~stringr::str_replace(.x, prefix_1L_chr, 
                  "")) %>% as.numeric())
    }) %>% stats::setNames(names(tbs_ls))
    return(tbs_ls)
}
#' Add unwtd dimension totals
#' @description add_unwtd_dim_tots() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add unwtd dimension totals. Function argument items_tb specifies the object to be updated. The function returns Items and domains (a tibble).
#' @param items_tb Items (a tibble)
#' @param domain_items_ls Domain items (a list)
#' @param domain_pfx_1L_chr Domain prefix (a character vector of length one)
#' @return Items and domains (a tibble)
#' @rdname add_unwtd_dim_tots
#' @export 
#' @importFrom purrr reduce
#' @importFrom dplyr mutate select
#' @importFrom rlang sym
add_unwtd_dim_tots <- function (items_tb, domain_items_ls, domain_pfx_1L_chr) 
{
    items_and_domains_tb <- purrr::reduce(1:length(domain_items_ls), 
        .init = items_tb, ~.x %>% dplyr::mutate(`:=`(!!rlang::sym(paste0(domain_pfx_1L_chr, 
            names(domain_items_ls)[.y])), rowSums(dplyr::select(., 
            domain_items_ls[[.y]])))))
    return(items_and_domains_tb)
}
#' Add wtd dimension totals
#' @description add_wtd_dim_tots() is an Add function that updates an object by adding data to that object. Specifically, this function implements an algorithm to add wtd dimension totals. Function argument unwtd_dim_tb specifies the object to be updated. The function returns Wtd and unwtd dimension (a tibble).
#' @param unwtd_dim_tb Unwtd dimension (a tibble)
#' @param domain_items_ls Domain items (a list)
#' @param domain_unwtd_pfx_1L_chr Domain unwtd prefix (a character vector of length one)
#' @param domain_wtd_pfx_1L_chr Domain wtd prefix (a character vector of length one)
#' @return Wtd and unwtd dimension (a tibble)
#' @rdname add_wtd_dim_tots
#' @export 
#' @importFrom utils data
#' @importFrom purrr map_dbl map2_dbl discard reduce
#' @importFrom dplyr filter pull select_if mutate
#' @importFrom rlang sym
#' @keywords internal
add_wtd_dim_tots <- function (unwtd_dim_tb, domain_items_ls, domain_unwtd_pfx_1L_chr, 
    domain_wtd_pfx_1L_chr) 
{
    utils::data("aqol6d_adult_disv_lup_tb", envir = environment())
    utils::data("aqol6d_domain_qs_lup_tb", envir = environment())
    min_vals_dbl <- purrr::map_dbl(domain_items_ls, ~length(.x)) %>% 
        unname()
    max_vals_dbl <- purrr::map2_dbl(domain_items_ls, names(domain_items_ls), 
        ~{
            paste0("Q", aqol6d_domain_qs_lup_tb %>% dplyr::filter(Domain_chr == 
                .y) %>% dplyr::pull(Question_dbl)) %>% purrr::map_dbl(~{
                tb <- aqol6d_adult_disv_lup_tb %>% dplyr::filter(Question_chr == 
                  .x) %>% dplyr::select_if(is.numeric)
                as.numeric(as.data.frame(tb)[1, ]) %>% purrr::discard(is.na) %>% 
                  length()
            }) %>% sum()
        }) %>% unname()
    wtd_and_unwtd_dim_tb <- purrr::reduce(1:length(domain_items_ls), 
        .init = unwtd_dim_tb, ~.x %>% dplyr::mutate(`:=`(!!rlang::sym(paste0(domain_wtd_pfx_1L_chr, 
            names(domain_items_ls)[.y])), (1 - (!!rlang::sym(paste0(domain_unwtd_pfx_1L_chr, 
            names(domain_items_ls)[.y])) - min_vals_dbl[.y])/(max_vals_dbl[.y] - 
            min_vals_dbl[.y])))))
    return(wtd_and_unwtd_dim_tb)
}
