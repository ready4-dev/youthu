#' Make balanced fake dataset
#' @description make_balanced_fake_ds() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make balanced fake dataset. The function returns Dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param match_on_vars_chr Match on variables (a character vector)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'UID_chr'
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'Timepoint_chr'
#' @param timepoint_bl_val_1L_chr Timepoint baseline value (a character vector of length one), Default: 'Baseline'
#' @param cmprsn_var_nm_1L_chr Comparison variable name (a character vector of length one), Default: 'study_arm_chr'
#' @param cmprsn_groups_chr Comparison groups (a character vector), Default: c("Intervention", "Control")
#' @return Dataset (a tibble)
#' @rdname make_balanced_fake_ds
#' @export 

#' @keywords internal
make_balanced_fake_ds <- function (ds_tb, match_on_vars_chr, id_var_nm_1L_chr = "UID_chr", 
    round_var_nm_1L_chr = "Timepoint_chr", timepoint_bl_val_1L_chr = "Baseline", 
    cmprsn_var_nm_1L_chr = "study_arm_chr", cmprsn_groups_chr = c("Intervention", 
        "Control")) 
{
    ds_tb <- ds_tb %>% transform_ds_for_cmprsn(id_var_nm_1L_chr = id_var_nm_1L_chr, 
        round_var_nm_1L_chr = round_var_nm_1L_chr, cmprsn_var_nm_1L_chr = cmprsn_var_nm_1L_chr, 
        cmprsn_groups_chr = cmprsn_groups_chr) %>% make_matched_ds_spine(round_var_nm_1L_chr = round_var_nm_1L_chr, 
        timepoint_bl_val_1L_chr = timepoint_bl_val_1L_chr, cmprsn_var_nm_1L_chr = cmprsn_var_nm_1L_chr, 
        active_arm_val_1L_chr = cmprsn_groups_chr[1], id_var_nm_1L_chr = id_var_nm_1L_chr, 
        match_on_vars_chr = match_on_vars_chr)
    return(ds_tb)
}
#' Make costs vector from gamma distribution
#' @description make_costs_vec_from_gamma_dstr() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make costs vector from gamma distribution. The function returns Costs (a double vector).
#' @param n_int N (an integer vector)
#' @param costs_mean_dbl Costs mean (a double vector)
#' @param costs_sd_dbl Costs standard deviation (a double vector)
#' @return Costs (a double vector)
#' @rdname make_costs_vec_from_gamma_dstr
#' @export 
#' @importFrom stats rgamma
#' @keywords internal
make_costs_vec_from_gamma_dstr <- function (n_int, costs_mean_dbl, costs_sd_dbl) 
{
    scale_1L_dbl <- costs_sd_dbl^2/costs_mean_dbl
    shape_1L_dbl <- costs_mean_dbl/scale_1L_dbl
    costs_dbl <- stats::rgamma(n_int, shape = shape_1L_dbl, scale = scale_1L_dbl)
    return(costs_dbl)
}
#' Make cost efns summary
#' @description make_cst_efns_smry() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make cost efns summary. The function returns Summary (a double vector).
#' @param ds_tb Dataset (a tibble)
#' @param idxs_int Indices (an integer vector)
#' @param change_types_chr Change types (a character vector), Default: 'dbl'
#' @param benefits_pfx_1L_chr Benefits prefix (a character vector of length one), Default: 'qalys_dbl'
#' @param benefits_var_nm_1L_chr Benefits variable name (a character vector of length one), Default: 'qalys'
#' @param costs_pfx_1L_chr Costs prefix (a character vector of length one), Default: 'costs_dbl'
#' @param costs_var_nm_1L_chr Costs variable name (a character vector of length one), Default: 'costs'
#' @param change_sfx_1L_chr Change suffix (a character vector of length one), Default: 'change'
#' @param change_vars_chr Change variables (a character vector), Default: 'NA'
#' @param cmprsn_groups_chr Comparison groups (a character vector), Default: c("Intervention", "Control")
#' @param cmprsn_var_nm_1L_chr Comparison variable name (a character vector of length one), Default: 'study_arm_chr'
#' @param round_fup_val_1L_chr Round follow-up value (a character vector of length one), Default: 'Follow-up'
#' @return Summary (a double vector)
#' @rdname make_cst_efns_smry
#' @export 
#' @importFrom tibble tibble
#' @importFrom tidyselect all_of
#' @importFrom purrr map
#' @importFrom dplyr group_by select all_of rename_with summarise across ungroup mutate
#' @importFrom rlang sym
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom tidyr pivot_wider
#' @keywords internal
make_cst_efns_smry <- function (ds_tb, idxs_int, change_types_chr = "dbl", benefits_pfx_1L_chr = "qalys_dbl", 
    benefits_var_nm_1L_chr = "qalys", costs_pfx_1L_chr = "costs_dbl", 
    costs_var_nm_1L_chr = "costs", change_sfx_1L_chr = "change", 
    change_vars_chr = NA_character_, cmprsn_groups_chr = c("Intervention", 
        "Control"), cmprsn_var_nm_1L_chr = "study_arm_chr", round_fup_val_1L_chr = "Follow-up") 
{
    if (!is.na(change_vars_chr[1])) {
        change_vars_with_sfx_chr <- paste0(change_vars_chr, "_", 
            change_sfx_1L_chr, "_", change_types_chr)
        replacements_chr <- paste0(change_vars_chr, "_", change_sfx_1L_chr)
    }
    else {
        change_vars_with_sfx_chr <- replacements_chr <- character(0)
    }
    selected_cols_chr <- c(costs_pfx_1L_chr, benefits_pfx_1L_chr, 
        change_vars_with_sfx_chr) %>% paste0("_", round_fup_val_1L_chr)
    rename_lup <- tibble::tibble(old_name_chr = tidyselect::all_of(selected_cols_chr), 
        new_name_chr = c(costs_var_nm_1L_chr, benefits_var_nm_1L_chr, 
            replacements_chr))
    new_nms_ls <- rename_lup$new_name_chr %>% purrr::map(~paste0(.x, 
        "_", cmprsn_groups_chr))
    summary_tb <- ds_tb[idxs_int, ] %>% dplyr::group_by(!!rlang::sym(cmprsn_var_nm_1L_chr)) %>% 
        dplyr::select(dplyr::all_of(c(selected_cols_chr, cmprsn_var_nm_1L_chr))) %>% 
        dplyr::rename_with(.cols = tidyselect::all_of(selected_cols_chr), 
            ~ready4fun::get_from_lup_obj(rename_lup, match_var_nm_1L_chr = "old_name_chr", 
                match_value_xx = .x, target_var_nm_1L_chr = "new_name_chr", 
                evaluate_lgl = F)) %>% dplyr::summarise(dplyr::across(.fns = mean)) %>% 
        dplyr::ungroup() %>% tidyr::pivot_wider(names_from = tidyselect::all_of(cmprsn_var_nm_1L_chr), 
        values_from = rename_lup$new_name_chr) %>% dplyr::mutate(`:=`(!!rlang::sym(paste0("diff_", 
        costs_var_nm_1L_chr, "_dbl")), !!rlang::sym(new_nms_ls[[1]][1]) - 
        !!rlang::sym(new_nms_ls[[1]][2])), `:=`(!!rlang::sym(paste0("diff_", 
        benefits_var_nm_1L_chr, "_dbl")), (!!rlang::sym(new_nms_ls[[2]][1]) - 
        !!rlang::sym(new_nms_ls[[2]][2]))), ICER = (!!rlang::sym(paste0("diff_", 
        costs_var_nm_1L_chr, "_dbl"))/!!rlang::sym(paste0("diff_", 
        benefits_var_nm_1L_chr, "_dbl")))) %>% dplyr::rename_with(~paste0(.x, 
        "_dbl")) %>% t()
    summary_dbl <- summary_tb[, 1]
    return(summary_dbl)
}
#' Make fake dataset one
#' @description make_fake_ds_one() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make fake dataset one. The function returns Fake data (a tibble).

#' @return Fake data (a tibble)
#' @rdname make_fake_ds_one
#' @export 
#' @importFrom youthvars transform_raw_ds_for_analysis
#' @importFrom dplyr select arrange rename mutate
#' @importFrom tibble as_tibble
#' @keywords internal
make_fake_ds_one <- function () 
{
    data("replication_popl_tb", package = "youthvars", envir = environment())
    fake_data_tb <- replication_popl_tb %>% youthvars::transform_raw_ds_for_analysis() %>% 
        dplyr::select(fkClientID, round, PHQ9, SOFAS) %>% dplyr::arrange(fkClientID) %>% 
        na.omit() %>% dplyr::rename(UID = fkClientID, Timepoint = round, 
        PHQ_total = PHQ9, SOFAS_total = SOFAS) %>% dplyr::mutate(SOFAS_total = as.integer(round(SOFAS_total, 
        0))) %>% tibble::as_tibble()
    return(fake_data_tb)
}
#' Make fake dataset two
#' @description make_fake_ds_two() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make fake dataset two. The function returns Matched dataset (a tibble).

#' @return Matched dataset (a tibble)
#' @rdname make_fake_ds_two
#' @export 
#' @importFrom youthvars transform_raw_ds_for_analysis
#' @importFrom dplyr filter pull
#' @importFrom lubridate days
#' @importFrom truncnorm rtruncnorm
#' @importFrom tibble tibble
#' @importFrom stats rnorm
#' @keywords internal
make_fake_ds_two <- function () 
{
    data("replication_popl_tb", package = "youthvars", envir = environment())
    seed_ds_tb <- replication_popl_tb %>% youthvars::transform_raw_ds_for_analysis() %>% 
        dplyr::filter(fkClientID %in% (replication_popl_tb %>% 
            dplyr::filter(round == "Baseline" & PHQ9 < 20) %>% 
            dplyr::pull(fkClientID)))
    ds_smry_ls <- list(bl_start_date_dtm = Sys.Date() - lubridate::days(300), 
        bl_end_date_dtm = Sys.Date() - lubridate::days(120), 
        cmprsn_var_nm_1L_chr = "study_arm_chr", cmprsn_groups_chr = c("Intervention", 
            "Control"), costs_mean_dbl = c(400, 1500), costs_sd_dbl = c(100, 
            220), costs_var_nm_1L_chr = "costs_dbl", date_var_nm_1L_chr = "date_psx", 
        duration_args_ls = list(a = 160, b = 220, mean = 180, 
            sd = 7), duration_fn = truncnorm::rtruncnorm, id_var_nm_1L_chr = "fkClientID", 
        predr_var_nms = c("PHQ9", "SOFAS"), round_var_nm_1L_chr = "round", 
        round_lvls_chr = c("Baseline", "Follow-up"), utl_var_nm_1L_chr = "AQoL6D_HU")
    sngl_grp_ds_tb <- make_sngl_grp_ds(seed_ds_tb, ds_smry_ls = ds_smry_ls)
    matched_ds_tb <- make_matched_ds(sngl_grp_ds_tb, cmprsn_smry_tb = tibble::tibble(var_nms_chr = c(ds_smry_ls$predr_var_nms, 
        ds_smry_ls$costs_var_nm_1L_chr), fns_ls = list(stats::rnorm, 
        stats::rnorm, stats::rnorm), abs_mean_diff_dbl = c(2, 
        2, 300), diff_sd_dbl = c(2, 2, 200), multiplier_dbl = c(-1, 
        -1, 1), min_dbl = c(0, 0, 0), max_dbl = c(27, 100, Inf), 
        integer_lgl = c(T, T, F)), ds_smry_ls = ds_smry_ls)
    return(matched_ds_tb)
}
#' Make fake trial dataset
#' @description make_fake_trial_ds() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make fake trial dataset. The function returns Updated dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'round'
#' @param round_lvls_chr Round levels (a character vector), Default: c("Baseline", "Follow-up")
#' @param match_on_vars_chr Match on variables (a character vector)
#' @param cmprsn_var_nm_1L_chr Comparison variable name (a character vector of length one), Default: 'study_arm_chr'
#' @param cmprsn_groups_chr Comparison groups (a character vector), Default: c("Intervention", "Control")
#' @param fns_ls Functions (a list)
#' @param var_nms_chr Variable names (a character vector)
#' @param abs_mean_diff_dbl Absolute mean difference (a double vector)
#' @param diff_sd_dbl Difference standard deviation (a double vector)
#' @param multiplier_dbl Multiplier (a double vector)
#' @param min_dbl Minimum (a double vector)
#' @param max_dbl Maximum (a double vector)
#' @param integer_lgl Integer (a logical vector)
#' @param match_idx_var_nm_1L_chr Match index variable name (a character vector of length one), Default: 'match_idx_int'
#' @return Updated dataset (a tibble)
#' @rdname make_fake_trial_ds
#' @export 

#' @keywords internal
make_fake_trial_ds <- function (ds_tb, id_var_nm_1L_chr = "fkClientID", round_var_nm_1L_chr = "round", 
    round_lvls_chr = c("Baseline", "Follow-up"), match_on_vars_chr, 
    cmprsn_var_nm_1L_chr = "study_arm_chr", cmprsn_groups_chr = c("Intervention", 
        "Control"), fns_ls, var_nms_chr, abs_mean_diff_dbl, diff_sd_dbl, 
    multiplier_dbl, min_dbl, max_dbl, integer_lgl, match_idx_var_nm_1L_chr = "match_idx_int") 
{
    updated_ds_tb <- ds_tb %>% make_balanced_fake_ds(id_var_nm_1L_chr = id_var_nm_1L_chr, 
        round_var_nm_1L_chr = round_var_nm_1L_chr, timepoint_bl_val_1L_chr = round_lvls_chr[1], 
        match_on_vars_chr = match_on_vars_chr, cmprsn_var_nm_1L_chr = cmprsn_var_nm_1L_chr, 
        cmprsn_groups_chr = cmprsn_groups_chr) %>% add_diffs_by_group_and_tmpt(cmprsn_var_nm_1L_chr = cmprsn_var_nm_1L_chr, 
        cmprsn_group_match_val_chr = cmprsn_groups_chr[1], round_var_nm_1L_chr = round_var_nm_1L_chr, 
        timepoint_match_val_1L_chr = round_lvls_chr[2], match_idx_var_nm_1L_chr = match_idx_var_nm_1L_chr, 
        var_nms_chr = var_nms_chr, fns_ls = fns_ls, abs_mean_diff_dbl = abs_mean_diff_dbl, 
        diff_sd_dbl = diff_sd_dbl, multiplier_dbl = multiplier_dbl, 
        min_dbl = min_dbl, max_dbl = max_dbl, integer_lgl = integer_lgl)
    return(updated_ds_tb)
}
#' Make health economic summary
#' @description make_hlth_ec_smry() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make health economic summary. The function returns He summary (a list).
#' @param ds_tb Dataset (a tibble)
#' @param predn_ds_ls Prediction dataset (a list)
#' @param wtp_dbl Willingness to pay (a double vector), Default: 50000
#' @param bootstrap_iters_1L_int Bootstrap iterations (an integer vector of length one), Default: 1000
#' @param benefits_pfx_1L_chr Benefits prefix (a character vector of length one), Default: 'qalys_dbl'
#' @param benefits_var_nm_1L_chr Benefits variable name (a character vector of length one), Default: 'qalys'
#' @param costs_var_nm_1L_chr Costs variable name (a character vector of length one), Default: 'costs'
#' @param change_sfx_1L_chr Change suffix (a character vector of length one), Default: 'change'
#' @return He summary (a list)
#' @rdname make_hlth_ec_smry
#' @export 
#' @importFrom boot boot
#' @importFrom purrr map_int
#' @importFrom BCEA bcea
make_hlth_ec_smry <- function (ds_tb, predn_ds_ls, wtp_dbl = 50000, bootstrap_iters_1L_int = 1000, 
    benefits_pfx_1L_chr = "qalys_dbl", benefits_var_nm_1L_chr = "qalys", 
    costs_var_nm_1L_chr = "costs", change_sfx_1L_chr = "change") 
{
    if (is.null(predn_ds_ls$ds_ls$predr_vars_nms_chr)) 
        predn_ds_ls$ds_ls$predr_vars_nms_chr <- predn_ds_ls$mdl_ls$predictors_lup$short_name_chr
    costs_pfx_1L_chr = predn_ds_ls$ds_ls$costs_var_nm_1L_chr
    cmprsn_groups_chr = predn_ds_ls$ds$cmprsn_groups_chr
    cmprsn_var_nm_1L_chr = predn_ds_ls$ds$cmprsn_var_nm_1L_chr
    round_fup_val_1L_chr = predn_ds_ls$ds$round_fup_val_1L_chr
    change_vars_chr <- c(predn_ds_ls$ds_ls$predr_vars_nms_chr, 
        predn_ds_ls$ds_ls$utl_var_nm_1L_chr)
    change_types_chr <- rep("dbl", length(change_vars_chr))
    costs_pfx_1L_chr <- predn_ds_ls$ds_ls$costs_var_nm_1L_chr
    bootstraps_ls <- boot::boot(ds_tb, make_cst_efns_smry, R = bootstrap_iters_1L_int, 
        benefits_pfx_1L_chr = benefits_pfx_1L_chr, costs_pfx_1L_chr = costs_pfx_1L_chr, 
        change_vars_chr = change_vars_chr, change_sfx_1L_chr = change_sfx_1L_chr, 
        change_types_chr = change_types_chr, cmprsn_groups_chr = cmprsn_groups_chr, 
        cmprsn_var_nm_1L_chr = cmprsn_var_nm_1L_chr, round_fup_val_1L_chr = round_fup_val_1L_chr, 
        benefits_var_nm_1L_chr = benefits_var_nm_1L_chr, costs_var_nm_1L_chr = costs_var_nm_1L_chr)
    costs_mat <- bootstraps_ls$t[, 1:2]
    benefits_mat <- bootstraps_ls$t[, 3:4]
    reordered_cmprsn_chr <- cmprsn_groups_chr[cmprsn_groups_chr %>% 
        purrr::map_int(~which(endsWith(names(bootstraps_ls$t0)[1:2], 
            paste0(.x, "_dbl"))))]
    ce_res_ls <- BCEA::bcea(e = benefits_mat, c = costs_mat, 
        ref = which(reordered_cmprsn_chr == cmprsn_groups_chr[1]), 
        interventions = reordered_cmprsn_chr, Kmax = wtp_dbl * 
            2)
    named_mat <- bootstraps_ls$t
    colnames(named_mat) <- names(bootstraps_ls$t0)
    he_smry_ls <- list(ce_res_ls = ce_res_ls, bootstraps_ls = bootstraps_ls, 
        benefits_mat = benefits_mat, costs_mat = costs_mat, named_mat = named_mat)
    return(he_smry_ls)
}
#' Make matched dataset
#' @description make_matched_ds() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make matched dataset. The function returns Matched dataset (a tibble).
#' @param sngl_grp_ds_tb Single group dataset (a tibble)
#' @param cmprsn_smry_tb Comparison summary (a tibble)
#' @param ds_smry_ls Dataset summary (a list)
#' @return Matched dataset (a tibble)
#' @rdname make_matched_ds
#' @export 
#' @importFrom dplyr mutate across
make_matched_ds <- function (sngl_grp_ds_tb, cmprsn_smry_tb, ds_smry_ls) 
{
    matched_ds_tb <- sngl_grp_ds_tb %>% make_fake_trial_ds(id_var_nm_1L_chr = ds_smry_ls$id_var_nm_1L_chr, 
        round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr, 
        round_lvls_chr = ds_smry_ls$round_lvls_chr, match_on_vars_chr = cmprsn_smry_tb$var_nms_chr, 
        cmprsn_var_nm_1L_chr = ds_smry_ls$cmprsn_var_nm_1L_chr, 
        cmprsn_groups_chr = ds_smry_ls$cmprsn_groups_chr, fns_ls = cmprsn_smry_tb$fns_ls, 
        var_nms_chr = cmprsn_smry_tb$var_nms_chr, abs_mean_diff_dbl = cmprsn_smry_tb$abs_mean_diff_dbl, 
        diff_sd_dbl = cmprsn_smry_tb$diff_sd_dbl, multiplier_dbl = cmprsn_smry_tb$multiplier_dbl, 
        min_dbl = cmprsn_smry_tb$min_dbl, max_dbl = cmprsn_smry_tb$max_dbl, 
        integer_lgl = cmprsn_smry_tb$integer_lgl)
    if (any(cmprsn_smry_tb$integer_lgl)) 
        matched_ds_tb <- matched_ds_tb %>% dplyr::mutate(dplyr::across(cmprsn_smry_tb$var_nms_chr[cmprsn_smry_tb$integer_lgl], 
            as.integer))
    return(matched_ds_tb)
}
#' Make matched dataset spine
#' @description make_matched_ds_spine() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make matched dataset spine. The function returns Matched dataset (a tibble).
#' @param ds_tb Dataset (a tibble)
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one), Default: 'Timepoint_chr'
#' @param timepoint_bl_val_1L_chr Timepoint baseline value (a character vector of length one), Default: 'Baseline'
#' @param cmprsn_var_nm_1L_chr Comparison variable name (a character vector of length one), Default: 'study_arm_chr'
#' @param active_arm_val_1L_chr Active arm value (a character vector of length one), Default: 'Intervention'
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'fkClientID'
#' @param match_on_vars_chr Match on variables (a character vector)
#' @return Matched dataset (a tibble)
#' @rdname make_matched_ds_spine
#' @export 
#' @importFrom dplyr filter mutate case_when arrange n select right_join
#' @importFrom rlang sym
#' @importFrom youthvars make_formula
#' @importFrom MatchIt matchit match.data
#' @importFrom purrr map_dfr
#' @keywords internal
make_matched_ds_spine <- function (ds_tb, round_var_nm_1L_chr = "Timepoint_chr", timepoint_bl_val_1L_chr = "Baseline", 
    cmprsn_var_nm_1L_chr = "study_arm_chr", active_arm_val_1L_chr = "Intervention", 
    id_var_nm_1L_chr = "fkClientID", match_on_vars_chr) 
{
    match_ds <- ds_tb %>% dplyr::filter(!!rlang::sym(round_var_nm_1L_chr) == 
        timepoint_bl_val_1L_chr) %>% dplyr::mutate(Intervention_lgl = dplyr::case_when(!!rlang::sym(cmprsn_var_nm_1L_chr) == 
        active_arm_val_1L_chr ~ T, T ~ F))
    matched_ls <- youthvars::make_formula("Intervention_lgl", 
        predictors_chr = match_on_vars_chr) %>% MatchIt::matchit(data = match_ds, 
        method = "nearest", ratio = 1)
    matched_ds_tb <- MatchIt::match.data(matched_ls)
    match_key_ds_tb <- c(F, T) %>% purrr::map_dfr(~matched_ds_tb %>% 
        dplyr::filter(Intervention_lgl == .x) %>% dplyr::arrange(distance) %>% 
        dplyr::mutate(match_idx_int = 1:dplyr::n())) %>% dplyr::arrange(match_idx_int) %>% 
        dplyr::select(!!rlang::sym(id_var_nm_1L_chr), study_arm_chr, 
            match_idx_int)
    matched_ds_tb <- dplyr::right_join(ds_tb, match_key_ds_tb) %>% 
        dplyr::arrange(match_idx_int)
    return(matched_ds_tb)
}
#' Make prediction metadata
#' @description make_predn_metadata_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make prediction metadata list. The function returns Prediction metadata (a list).
#' @param data_tb Data (a tibble)
#' @param cmprsn_groups_chr Comparison groups (a character vector), Default: NULL
#' @param cmprsn_var_nm_1L_chr Comparison variable name (a character vector of length one), Default: NULL
#' @param costs_var_nm_1L_chr Costs variable name (a character vector of length one), Default: NULL
#' @param id_var_nm_1L_chr Identity variable name (a character vector of length one), Default: 'UID'
#' @param mdl_meta_data_ls Model meta data (a list), Default: NULL
#' @param mdls_lup Models (a lookup table), Default: NULL
#' @param mdl_nm_1L_chr Model name (a character vector of length one), Default: NULL
#' @param msrmnt_date_var_nm_1L_chr Measurement date variable name (a character vector of length one), Default: NULL
#' @param predr_vars_nms_chr Predictor variables names (a character vector), Default: NULL
#' @param round_var_nm_1L_chr Round variable name (a character vector of length one)
#' @param round_bl_val_1L_chr Round baseline value (a character vector of length one)
#' @param utl_var_nm_1L_chr Utility variable name (a character vector of length one), Default: 'AQoL6D_HU'
#' @param server_1L_chr Server (a character vector of length one), Default: 'dataverse.harvard.edu'
#' @param key_1L_chr Key (a character vector of length one), Default: NULL
#' @return Prediction metadata (a list)
#' @rdname make_predn_metadata_ls
#' @export 
#' @importFrom TTU rename_from_nmd_vec transform_mdl_vars_with_clss
#' @importFrom purrr walk pluck
#' @importFrom dplyr pull filter
#' @importFrom ready4fun get_from_lup_obj
#' @importFrom assertthat assert_that
#' @importFrom rlang sym
#' @importFrom lubridate is.Date
#' @keywords internal
make_predn_metadata_ls <- function (data_tb, cmprsn_groups_chr = NULL, cmprsn_var_nm_1L_chr = NULL, 
    costs_var_nm_1L_chr = NULL, id_var_nm_1L_chr = "UID", mdl_meta_data_ls = NULL, 
    mdls_lup = NULL, mdl_nm_1L_chr = NULL, msrmnt_date_var_nm_1L_chr = NULL, 
    predr_vars_nms_chr = NULL, round_var_nm_1L_chr, round_bl_val_1L_chr, 
    utl_var_nm_1L_chr = "AQoL6D_HU", server_1L_chr = "dataverse.harvard.edu", 
    key_1L_chr = NULL) 
{
    if (!is.null(predr_vars_nms_chr)) {
        data_tb <- TTU::rename_from_nmd_vec(data_tb, nmd_vec_chr = predr_vars_nms_chr, 
            vec_nms_as_new_1L_lgl = T)
    }
    predictors_lup <- get_predictors_lup(mdl_meta_data_ls = mdl_meta_data_ls, 
        mdls_lup = mdls_lup, mdl_nm_1L_chr = mdl_nm_1L_chr, outp_is_abbrvs_tb = F, 
        server_1L_chr = server_1L_chr, key_1L_chr = key_1L_chr)
    data_tb <- data_tb %>% TTU::transform_mdl_vars_with_clss(predictors_lup = predictors_lup)
    purrr::walk(predictors_lup$short_name_chr, ~{
        vector_xx <- data_tb %>% dplyr::pull(.x)
        min_dbl <- ready4fun::get_from_lup_obj(predictors_lup, 
            match_value_xx = .x, match_var_nm_1L_chr = "short_name_chr", 
            target_var_nm_1L_chr = "min_val_dbl", evaluate_lgl = F)
        max_dbl <- ready4fun::get_from_lup_obj(predictors_lup, 
            match_value_xx = .x, match_var_nm_1L_chr = "short_name_chr", 
            target_var_nm_1L_chr = "max_val_dbl", evaluate_lgl = F)
        assertthat::assert_that(max(vector_xx) <= max_dbl & min(vector_xx) >= 
            min_dbl, msg = paste0(predr_vars_nms_chr %>% purrr::pluck(.x), 
            " variable must be within range of ", min_dbl, " and ", 
            max_dbl, "."))
    })
    assertthat::assert_that(length(data_tb %>% dplyr::pull(!!rlang::sym(round_var_nm_1L_chr)) %>% 
        unique()) == 2, msg = paste0(round_var_nm_1L_chr, " variable must have two values - one for each data collection round."))
    assertthat::assert_that(is.factor(data_tb %>% dplyr::pull(!!rlang::sym(round_var_nm_1L_chr))), 
        msg = paste0(round_var_nm_1L_chr, " variable must be a factor."))
    assertthat::assert_that(round_bl_val_1L_chr %in% levels(data_tb %>% 
        dplyr::pull(!!rlang::sym(round_var_nm_1L_chr))), msg = paste0("Levels of ", 
        round_var_nm_1L_chr, " variable must include one that is named ", 
        round_bl_val_1L_chr, "."))
    bl_ds_tb <- data_tb %>% dplyr::filter(!!rlang::sym(round_var_nm_1L_chr) == 
        round_bl_val_1L_chr)
    fup_ds_tb <- data_tb %>% dplyr::filter(!!rlang::sym(round_var_nm_1L_chr) != 
        round_bl_val_1L_chr)
    bl_uids_xx <- bl_ds_tb %>% dplyr::pull(id_var_nm_1L_chr)
    fup_uids_xx <- fup_ds_tb %>% dplyr::pull(id_var_nm_1L_chr)
    assertthat::assert_that(nrow(bl_ds_tb) == length(unique(bl_uids_xx)) & 
        nrow(fup_ds_tb) == length(unique(fup_uids_xx)), msg = paste0("At each time-point (the ", 
        round_var_nm_1L_chr, " variable) there must be one unique record identifier (the ", 
        id_var_nm_1L_chr, " variable)."))
    if (!is.null(msrmnt_date_var_nm_1L_chr)) 
        assertthat::assert_that(lubridate::is.Date(data_tb %>% 
            dplyr::pull(msrmnt_date_var_nm_1L_chr)), msg = paste0(msrmnt_date_var_nm_1L_chr, 
            " variable must be of date class."))
    round_vals_chr <- data_tb %>% dplyr::pull(!!rlang::sym(round_var_nm_1L_chr)) %>% 
        levels()
    predn_metadata_ls <- list(ds_ls = list(cmprsn_groups_chr = cmprsn_groups_chr, 
        cmprsn_var_nm_1L_chr = cmprsn_var_nm_1L_chr, costs_var_nm_1L_chr = costs_var_nm_1L_chr, 
        id_var_nm_1L_chr = id_var_nm_1L_chr, msrmnt_date_var_nm_1L_chr = msrmnt_date_var_nm_1L_chr, 
        predr_vars_nms_chr = predr_vars_nms_chr, round_var_nm_1L_chr = round_var_nm_1L_chr, 
        round_bl_val_1L_chr = round_bl_val_1L_chr, round_fup_val_1L_chr = round_vals_chr[round_vals_chr != 
            round_bl_val_1L_chr], utl_var_nm_1L_chr = utl_var_nm_1L_chr), 
        mdl_ls = list(mdl_meta_data_ls = mdl_meta_data_ls, mdls_lup = mdls_lup, 
            mdl_nm_1L_chr = mdl_nm_1L_chr, predictors_lup = predictors_lup))
    return(predn_metadata_ls)
}
#' Make single group dataset
#' @description make_sngl_grp_ds() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make single group dataset. The function returns Single group dataset (a tibble).
#' @param seed_ds_tb Seed dataset (a tibble), Default: NULL
#' @param ds_smry_ls Dataset summary (a list)
#' @return Single group dataset (a tibble)
#' @rdname make_sngl_grp_ds
#' @export 
#' @importFrom dplyr select mutate arrange
#' @importFrom stats na.omit
#' @importFrom tibble as_tibble
make_sngl_grp_ds <- function (seed_ds_tb = NULL, ds_smry_ls) 
{
    sngl_grp_ds_tb <- seed_ds_tb %>% dplyr::select(ds_smry_ls$id_var_nm_1L_chr, 
        ds_smry_ls$round_var_nm_1L_chr, ds_smry_ls$predr_var_nms) %>% 
        stats::na.omit()
    if ("SOFAS" %in% ds_smry_ls$predr_var_nms) 
        sngl_grp_ds_tb <- sngl_grp_ds_tb %>% dplyr::mutate(SOFAS = as.integer(round(SOFAS, 
            0)))
    sngl_grp_ds_tb <- sngl_grp_ds_tb %>% tibble::as_tibble() %>% 
        add_dates_from_dstr(bl_start_date_dtm = ds_smry_ls$bl_start_date_dtm, 
            bl_end_date_dtm = ds_smry_ls$bl_end_date_dtm, duration_args_ls = ds_smry_ls$duration_args_ls, 
            duration_fn = ds_smry_ls$duration_fn, date_var_nm_1L_chr = ds_smry_ls$date_var_nm_1L_chr, 
            id_var_nm_1L_chr = ds_smry_ls$id_var_nm_1L_chr, round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr, 
            round_bl_val_1L_chr = ds_smry_ls$round_lvls_chr[1])
    if (!is.null(ds_smry_ls$costs_mean_dbl) & !is.null(ds_smry_ls$costs_sd_dbl)) 
        sngl_grp_ds_tb <- sngl_grp_ds_tb %>% add_costs_by_tmpt(round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr, 
            round_lvls_chr = ds_smry_ls$round_lvls_chr, costs_mean_dbl = ds_smry_ls$costs_mean_dbl, 
            costs_sd_dbl = ds_smry_ls$costs_sd_dbl, extra_cost_args_ls = list(costs_var_nm_1L_chr = ds_smry_ls$costs_var_nm_1L_chr), 
            fn = add_costs_from_gamma_dstr)
    sngl_grp_ds_tb <- sngl_grp_ds_tb %>% dplyr::arrange(ds_smry_ls$id_var_nm_1L_chr)
    return(sngl_grp_ds_tb)
}
