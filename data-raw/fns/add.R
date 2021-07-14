add_change_in_ds_var <- function(ds_tb,
                                 id_var_nm_1L_chr = "fkClientID",
                                 round_var_nm_1L_chr = "round",
                                 round_bl_val_1L_chr = "Baseline",
                                 change_var_nm_1L_chr,
                                 var_nm_1L_chr,
                                 arrange_by_id_lgl = T){
  updated_ds_tb <-  ds_tb %>%
    dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) %>%
    dplyr::mutate(!!rlang::sym(change_var_nm_1L_chr) := dplyr::case_when(!!rlang::sym(round_var_nm_1L_chr) == round_bl_val_1L_chr ~ 0,
                                                                         T ~ (as.numeric(!!rlang::sym(var_nm_1L_chr)) - dplyr::lag(as.numeric(!!rlang::sym(var_nm_1L_chr)))))) %>% dplyr::ungroup()
  if(arrange_by_id_lgl)
    updated_ds_tb <-  updated_ds_tb %>%
      dplyr::arrange(!!rlang::sym(id_var_nm_1L_chr))
  return(updated_ds_tb)
}
add_costs_by_tmpt <- function(ds_tb,
                              round_var_nm_1L_chr,
                              round_lvls_chr = c("Baseline","Follow-up"),
                              costs_mean_dbl,
                              costs_sd_dbl,
                              extra_cost_args_ls = list(costs_var_nm_1L_chr = "costs_dbl"),
                              fn = add_costs_from_gamma_dstr){
  updated_ds_tb <- purrr::pmap_dfr(list(round_lvls_chr,
                                        costs_mean_dbl,
                                        costs_sd_dbl),
                                   ~ {
                                     args_ls <- list(ds_tb %>%
                                                       dplyr::filter(!!rlang::sym(round_var_nm_1L_chr) == ..1),..2,..3)
                                     if(!is.null(extra_cost_args_ls))
                                       args_ls <- append(args_ls, extra_cost_args_ls)
                                     rlang::exec(.fn = fn, !!!args_ls)
                                   })
  return(updated_ds_tb)
}
add_costs_from_gamma_dstr <- function(ds_tb,
                                      costs_mean_dbl,
                                      costs_sd_dbl,
                                      costs_var_nm_1L_chr = "costs_dbl"){

  updated_ds_tb <- dplyr::mutate(ds_tb,
                                 !!rlang::sym(costs_var_nm_1L_chr) := make_costs_vec_from_gamma_dstr(n_int = nrow(ds_tb),
                                                                                                     costs_mean_dbl = costs_mean_dbl,
                                                                                                     costs_sd_dbl = costs_sd_dbl))
  return(updated_ds_tb)

}
add_dates_from_dstr <- function(ds_tb,
                                bl_start_date_dtm,
                                bl_end_date_dtm,
                                duration_args_ls,
                                duration_fn = stats::rnorm,
                                date_var_nm_1L_chr = "date_psx",
                                id_var_nm_1L_chr = "fkClientID",
                                round_var_nm_1L_chr = "round",
                                round_bl_val_1L_chr = "Baseline",
                                origin_1L_chr = '1970-01-01'){
  args_ls <- append(list(n=nrow(ds_tb)), duration_args_ls)
  days_of_fup_int <- rlang::exec(.fn = duration_fn, !!!args_ls) %>% round(0) %>% as.integer()
  updated_ds_tb <- ds_tb %>%
    dplyr::mutate(duration_prd = dplyr::case_when(!!rlang::sym(round_var_nm_1L_chr) != round_bl_val_1L_chr ~lubridate::days(days_of_fup_int),
                                                  T ~ lubridate::days(0))) %>%
    dplyr::mutate(!!rlang::sym(date_var_nm_1L_chr) := dplyr::case_when(!!rlang::sym(round_var_nm_1L_chr) == round_bl_val_1L_chr ~ as.Date(sample(as.numeric(bl_start_date_dtm):as.numeric(bl_end_date_dtm),
                                                                                                                                                 dplyr::n(),
                                                                                                                                                 replace = T),
                                                                                                                                          origin = origin_1L_chr ))) %>%
    dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) %>%
    dplyr::mutate(!!rlang::sym(date_var_nm_1L_chr) := dplyr::case_when(!!rlang::sym(round_var_nm_1L_chr) == round_bl_val_1L_chr ~ !!rlang::sym(date_var_nm_1L_chr),
                                                                       T ~ dplyr::lag(!!rlang::sym(date_var_nm_1L_chr)) + duration_prd)) %>%
    dplyr::ungroup()  %>%
    dplyr::select(!!rlang::sym(id_var_nm_1L_chr),
                  !!rlang::sym(round_var_nm_1L_chr),
                  !!rlang::sym(date_var_nm_1L_chr),
                  duration_prd,
                  dplyr::everything())
  return(updated_ds_tb)
}
add_diffs_by_group_and_tmpt <- function(ds_tb = trial_ds_tb,
                                        cmprsn_var_nm_1L_chr = "study_arm_chr",
                                        cmprsn_group_match_val_chr = c("Intervention"),
                                        round_var_nm_1L_chr = "round",
                                        timepoint_match_val_1L_chr = "Follow-up",
                                        match_idx_var_nm_1L_chr = "match_idx_int",
                                        var_nms_chr, # From this point on should be a tibble and moved to become second arg
                                        fns_ls,
                                        abs_mean_diff_dbl,
                                        diff_sd_dbl,
                                        multiplier_dbl,
                                        min_dbl,
                                        max_dbl,
                                        integer_lgl){
  unchanged_vals_tb <- ds_tb %>%
    dplyr::filter(!(!!rlang::sym(cmprsn_var_nm_1L_chr) == cmprsn_group_match_val_chr &
                      !!rlang::sym(round_var_nm_1L_chr) == timepoint_match_val_1L_chr))
  updated_vals_tb <- ds_tb %>%
    dplyr::filter(!!rlang::sym(cmprsn_var_nm_1L_chr) == cmprsn_group_match_val_chr &
                    !!rlang::sym(round_var_nm_1L_chr) == timepoint_match_val_1L_chr) %>%
    update_multpl_cols_with_diffs(var_nms_chr = var_nms_chr,
                                  fns_ls = fns_ls,
                                  abs_mean_diff_dbl = abs_mean_diff_dbl,
                                  diff_sd_dbl = diff_sd_dbl,
                                  multiplier_dbl = multiplier_dbl,
                                  min_dbl = min_dbl,
                                  max_dbl = max_dbl,
                                  integer_lgl = integer_lgl)
  updated_ds_tb <- dplyr::bind_rows(unchanged_vals_tb,
                                    updated_vals_tb) %>%
    dplyr::arrange(!!rlang::sym(match_idx_var_nm_1L_chr))
  return(updated_ds_tb)
}
add_qalys <- function(ds_tb,
                      cmprsn_var_nm_1L_chr = "study_arm_chr",
                      duration_var_nm_1L_chr = "duration_prd",
                      id_var_nm_1L_chr = "fkClientID",
                      match_idx_var_nm_1L_chr = "match_idx_int",
                      msrmnt_date_var_nm_1L_chr = "date_dtm",
                      qalys_var_nm_1L_chr = "qalys_dbl",
                      round_var_nm_1L_chr = "round",
                      round_bl_val_1L_chr = "Baseline",
                      utl_change_var_nm_1L_chr = "utl_change_dbl",
                      utl_var_nm_1L_chr = "utility_dbl",
                      reshape_1L_lgl = T){
  if(!duration_var_nm_1L_chr %in% names(ds_tb)){
    ds_tb <- ds_tb %>%
      dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) %>%
      dplyr::mutate(!!rlang::sym(duration_var_nm_1L_chr) := dplyr::case_when(!!rlang::sym(round_var_nm_1L_chr) != round_bl_val_1L_chr ~ lubridate::as.period(!!rlang::sym(msrmnt_date_var_nm_1L_chr) - dplyr::lag(!!rlang::sym(msrmnt_date_var_nm_1L_chr))),
                                                                             T ~ lubridate::days(0))) %>%
      dplyr::ungroup()
  }
  updated_ds_tb <- ds_tb %>%
    dplyr::mutate(!!rlang::sym(qalys_var_nm_1L_chr) := (!!rlang::sym(utl_var_nm_1L_chr)-(!!rlang::sym(utl_change_var_nm_1L_chr) * 0.5)) * (!!rlang::sym(duration_var_nm_1L_chr) / lubridate::years(1)))
  if(reshape_1L_lgl){
    vars_to_spread_chr <- names(updated_ds_tb)[!names(updated_ds_tb) %in% c(cmprsn_var_nm_1L_chr,
                                                                            id_var_nm_1L_chr,
                                                                            match_idx_var_nm_1L_chr,
                                                                            round_var_nm_1L_chr)]
    updated_ds_tb <- updated_ds_tb %>% tidyr::pivot_wider(names_from = !!rlang::sym(round_var_nm_1L_chr),
                                                          values_from = tidyselect::all_of(vars_to_spread_chr))
  }
  return(updated_ds_tb)
}
add_qalys_to_ds <- function(ds_tb,
                            predn_ds_ls,
                            include_predrs_1L_lgl = T,
                            reshape_1L_lgl = T){
  if(is.null(predn_ds_ls$ds_ls$predr_vars_nms_chr))
    predn_ds_ls$ds_ls$predr_vars_nms_chr <- predn_ds_ls$mdl_ls$predictors_lup$short_name_chr
  ds_smry_ls <- predn_ds_ls$ds_ls
  if(include_predrs_1L_lgl){
    predr_vars_nms_chr <- ds_smry_ls$predr_vars_nms_chr
  }else{
    predr_vars_nms_chr <- character(0)
  }
  args_ls_ls <- purrr::map(c(predr_vars_nms_chr,
                             ds_smry_ls$utl_var_nm_1L_chr),
                           ~ list(change_var_nm_1L_chr = paste0(.x,"_change_dbl"),
                                  var_nm_1L_chr = .x))
  ds_tb <- purrr::reduce(1:length(args_ls_ls),
                         .init = ds_tb,
                         ~ add_change_in_ds_var(.x,
                                                id_var_nm_1L_chr = predn_ds_ls$ds_ls$id_var_nm_1L_chr,
                                                round_bl_val_1L_chr = predn_ds_ls$ds_ls$round_bl_val_1L_chr,
                                                round_var_nm_1L_chr = predn_ds_ls$ds_ls$round_var_nm_1L_chr,
                                                var_nm_1L_chr = args_ls_ls[[.y]]$var_nm_1L_chr,
                                                change_var_nm_1L_chr = args_ls_ls[[.y]]$change_var_nm_1L_chr)) %>%
    add_qalys(id_var_nm_1L_chr = predn_ds_ls$ds_ls$id_var_nm_1L_chr,
              msrmnt_date_var_nm_1L_chr = predn_ds_ls$ds_ls$msrmnt_date_var_nm_1L_chr,
              round_bl_val_1L_chr = predn_ds_ls$ds_ls$round_bl_val_1L_chr,
              round_var_nm_1L_chr = predn_ds_ls$ds_ls$round_var_nm_1L_chr,
              utl_change_var_nm_1L_chr = paste0(ds_smry_ls$utl_var_nm_1L_chr,"_change_dbl"),
              utl_var_nm_1L_chr = ds_smry_ls$utl_var_nm_1L_chr,
              duration_var_nm_1L_chr = "duration_prd",
              qalys_var_nm_1L_chr = "qalys_dbl",
              reshape_1L_lgl = reshape_1L_lgl)
  return(ds_tb)
}
add_utl_predn <- function(data_tb,
                          predn_ds_ls,
                          deterministic_1L_lgl = T,
                          force_min_max_1L_lgl = T,
                          key_1L_chr = NULL,
                          make_from_tbl_1L_lgl = T,
                          model_mdl = NULL,
                          new_data_is_1L_chr = "Simulated",
                          server_1L_chr = "dataverse.harvard.edu",
                          utl_cls_fn = NULL              # ,
                          ){
  id_var_nm_1L_chr = predn_ds_ls$ds_ls$id_var_nm_1L_chr
  predr_vars_nms_chr = predn_ds_ls$ds_ls$predr_vars_nms_chr
  round_var_nm_1L_chr = predn_ds_ls$ds_ls$round_var_nm_1L_chr
  round_bl_val_1L_chr = predn_ds_ls$ds_ls$round_bl_val_1L_chr
  utl_var_nm_1L_chr = predn_ds_ls$ds_ls$utl_var_nm_1L_chr
  mdl_meta_data_ls = predn_ds_ls$mdl_ls$mdl_meta_data_ls
  mdls_lup = predn_ds_ls$mdl_ls$mdls_lup
  mdl_nm_1L_chr = predn_ds_ls$mdl_ls$mdl_nm_1L_chr
  if(is.null(model_mdl))
    model_mdl <- get_model(mdls_lup,
                           mdl_nm_1L_chr = mdl_nm_1L_chr,
                           make_from_tbl_1L_lgl = make_from_tbl_1L_lgl,
                           mdl_meta_data_ls = mdl_meta_data_ls,
                           server_1L_chr = server_1L_chr,
                           key_1L_chr = key_1L_chr)
  updated_tb <- TTU::add_utl_predn_to_new_ds(data_tb,
                                             mdl_nm_1L_chr = mdl_nm_1L_chr,
                                             id_var_nm_1L_chr = id_var_nm_1L_chr,
                                             analysis_1L_chr = ready4fun::get_from_lup_obj(mdls_lup,
                                                                                           match_value_xx = mdl_nm_1L_chr,
                                                                                           match_var_nm_1L_chr = "mdl_nms_chr",
                                                                                           target_var_nm_1L_chr = "source_chr",
                                                                                           evaluate_lgl = F),
                                             ingredients_ls = get_mdl_metadata(mdls_lup,
                                                                               mdl_nm_1L_chr = mdl_nm_1L_chr),
                                             model_mdl = model_mdl,
                                             new_data_is_1L_chr = new_data_is_1L_chr,
                                             predr_vars_nms_chr = predr_vars_nms_chr,
                                             round_var_nm_1L_chr = round_var_nm_1L_chr,
                                             round_bl_val_1L_chr = round_bl_val_1L_chr,
                                             utl_var_nm_1L_chr = utl_var_nm_1L_chr)
  return(updated_tb)
}
