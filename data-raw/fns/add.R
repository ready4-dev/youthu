add_aqol6d_predn_to_ds <- function(data_tb,
                                   model_mdl,
                                   tfmn_1L_chr,
                                   predr_vars_nms_chr = NULL,
                                   utl_var_nm_1L_chr = NULL,
                                   id_var_nm_1L_chr = "fkClientID",
                                   round_var_nm_1L_chr = "round",
                                   round_bl_val_1L_chr = "Baseline",
                                   utl_cls_fn = youthvars::youthvars_aqol6d_adol,
                                   predictors_lup = NULL){
  if (is.null(predictors_lup))
    utils::data("predictors_lup", envir = environment())
  if(!is.null(names(predr_vars_nms_chr))){
    data_tb <- rename_from_nmd_vec(data_tb,
                                   nmd_vec_chr = predr_vars_nms_chr,
                                   vec_nms_as_new_1L_lgl = T)
  }
  terms_ls <- model_mdl$terms
  mdl_dep_var_1L_chr <- terms_ls[[2]] %>% as.character()
  mdl_predr_terms_chr <- terms_ls[[3]] %>% as.character()
  mdl_predr_terms_chr <- mdl_predr_terms_chr %>% strsplit(split = " +") %>% purrr::flatten_chr()
  mdl_predr_terms_chr <- mdl_predr_terms_chr[mdl_predr_terms_chr!="+"]
  mdl_predr_terms_chr <- mdl_predr_terms_chr %>% purrr::map_chr(~stringr::str_replace(.x,"_baseline","") %>%                                                         stringr::str_replace("_change","")
  ) %>% unique()
  original_ds_vars_chr <- names(data_tb)[!names(data_tb) %in% c(mdl_predr_terms_chr,
                                        ifelse(!is.null(utl_var_nm_1L_chr),
                                               utl_var_nm_1L_chr,
                                               mdl_dep_var_1L_chr))]
  updated_tb <- data_tb %>%
    transform_ds_to_predn_ds(predr_vars_nms_chr = mdl_predr_terms_chr,
                             tfmn_1L_chr = tfmn_1L_chr,
                             dep_var_nm_1L_chr = mdl_dep_var_1L_chr,
                             id_var_nm_1L_chr = id_var_nm_1L_chr,
                             round_var_nm_1L_chr = round_var_nm_1L_chr,
                             round_bl_val_1L_chr = round_bl_val_1L_chr,
                             predictors_lup = predictors_lup) %>%
    TTU::add_utility_predn_to_ds(model_mdl = model_mdl,
                                 tfmn_1L_chr = tfmn_1L_chr,
                                 dep_var_nm_1L_chr = mdl_dep_var_1L_chr,
                                 predr_vars_nms_chr = mdl_predr_terms_chr,
                                 utl_cls_fn = utl_cls_fn,
                                 rmv_tfmd_dep_var_1L_lgl = T)
  if(!is.null(utl_var_nm_1L_chr)){
    updated_tb <- updated_tb %>%
      dplyr::rename(!!rlang::sym(utl_var_nm_1L_chr):=tidyselect::all_of(mdl_dep_var_1L_chr))
  }
  if(!is.null(names(predr_vars_nms_chr))){
    updated_tb <- rename_from_nmd_vec(updated_tb,
                                      nmd_vec_chr = predr_vars_nms_chr,
                                      vec_nms_as_new_1L_lgl = F)
  }

  if("aqol6d_total_w_CLL_cloglog" %in% names(updated_tb))
    updated_tb <- updated_tb %>%
    dplyr::select(-aqol6d_total_w_CLL_cloglog)
  names_to_incl_chr <- c(names(updated_tb),
                         setdiff(names(data_tb),
                                 names(updated_tb)))
  updated_tb <- dplyr::left_join(data_tb %>% dplyr::select(tidyselect::all_of(original_ds_vars_chr)),
                                 updated_tb)
  updated_tb <- updated_tb %>%
    dplyr::select(tidyselect::all_of(names_to_incl_chr[names_to_incl_chr %in% names(updated_tb)]))
  return(updated_tb)
}
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
                              fn = add_costs_from_gamma_dist){
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
add_costs_from_gamma_dist <- function(ds_tb,
                                      costs_mean_dbl,
                                      costs_sd_dbl,
                                      costs_var_nm_1L_chr = "costs_dbl"){

  updated_ds_tb <- dplyr::mutate(ds_tb,
                                 !!rlang::sym(costs_var_nm_1L_chr) := make_costs_vec_from_gamma_dist(n_int = nrow(ds_tb),
                                                                                                     costs_mean_dbl = costs_mean_dbl,
                                                                                                     costs_sd_dbl = costs_sd_dbl))
  return(updated_ds_tb)

}
add_dates_from_dist <- function(ds_tb,
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
                      qalys_var_nm_1L_chr = "qalys_dbl",
                      round_var_nm_1L_chr = "round",
                      utl_change_var_nm_1L_chr = "utl_change_dbl",
                      utl_var_nm_1L_chr = "utility_dbl",
                      reshape_1L_lgl = T){
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
                            ds_smry_ls){
  args_ls_ls <- purrr::map(c(ds_smry_ls$predr_var_nms,
                             ds_smry_ls$utl_var_nm_1L_chr),
                           ~ list(change_var_nm_1L_chr = paste0(.x,"_change_dbl"),
                                  var_nm_1L_chr = .x))
  ds_tb <- purrr::reduce(1:length(args_ls_ls),
                         .init = ds_tb,
                         ~ add_change_in_ds_var(.x,
                                                var_nm_1L_chr = args_ls_ls[[.y]]$var_nm_1L_chr,
                                                change_var_nm_1L_chr = args_ls_ls[[.y]]$change_var_nm_1L_chr)) %>%
    add_qalys(utl_change_var_nm_1L_chr = paste0(ds_smry_ls$utl_var_nm_1L_chr,"_change_dbl"),
              utl_var_nm_1L_chr = ds_smry_ls$utl_var_nm_1L_chr,
              duration_var_nm_1L_chr = "duration_prd",
              qalys_var_nm_1L_chr = "qalys_dbl")
  return(ds_tb)
}
