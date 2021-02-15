make_balanced_fake_ds <- function(ds_tb,
                                  match_on_vars_chr,
                                  id_var_nm_1L_chr = "UID_chr",
                                  round_var_nm_1L_chr = "Timepoint_chr",
                                  timepoint_bl_val_1L_chr = "Baseline",
                                  cmprsn_var_nm_1L_chr = "study_arm_chr",
                                  cmprsn_groups_chr = c("Intervention","Control")
){
    ds_tb <- ds_tb %>%
        transform_ds_for_cmprsn(id_var_nm_1L_chr = id_var_nm_1L_chr,
                                round_var_nm_1L_chr = round_var_nm_1L_chr,
                                cmprsn_var_nm_1L_chr = cmprsn_var_nm_1L_chr,
                                cmprsn_groups_chr = cmprsn_groups_chr) %>%
        make_matched_ds_spine(round_var_nm_1L_chr = round_var_nm_1L_chr,
                        timepoint_bl_val_1L_chr = timepoint_bl_val_1L_chr,
                        cmprsn_var_nm_1L_chr = cmprsn_var_nm_1L_chr,
                        active_arm_val_1L_chr = cmprsn_groups_chr[1],
                        id_var_nm_1L_chr = id_var_nm_1L_chr,
                        match_on_vars_chr = match_on_vars_chr)
    return(ds_tb)
}
make_ce_smry <- function(ds_tb,
                         indices,
                         change_types_chr = "dbl",
                         benefits_pfx_1L_chr = "qalys_dbl",
                         benefits_var_nm_1L_chr = "qalys",
                         costs_pfx_1L_chr = "costs_dbl",
                         costs_var_nm_1L_chr = "costs",
                         change_sfx_1L_chr = "change",
                         change_vars_chr = NA_character_,
                         cmprsn_groups_chr = c("Intervention","Control"),
                         cmprsn_var_nm_1L_chr = "study_arm_chr",
                         round_fup_1L_chr = "Follow-up"){
    if(!is.na(change_vars_chr[1])){
        change_vars_with_sfx_chr <- paste0(change_vars_chr, "_", change_sfx_1L_chr,"_", change_types_chr)
        replacements_chr <- paste0(change_vars_chr, "_", change_sfx_1L_chr)
    }else{
        change_vars_with_sfx_chr <- replacements_chr <-character(0)
    }
    selected_cols_chr <- c(costs_pfx_1L_chr,
                           benefits_pfx_1L_chr,
                           change_vars_with_sfx_chr) %>% paste0("_",round_fup_1L_chr)
    rename_lup <- tibble::tibble(old_name_chr = tidyselect::all_of(selected_cols_chr),
                                 new_name_chr = c(costs_var_nm_1L_chr, benefits_var_nm_1L_chr, replacements_chr))
    new_nms_ls <- rename_lup$new_name_chr %>% purrr::map(~paste0(.x,"_",cmprsn_groups_chr))
    summary_tb <- ds_tb[indices,] %>%
        dplyr::group_by(!!rlang::sym(cmprsn_var_nm_1L_chr)) %>%
        dplyr::select(dplyr::all_of(c(selected_cols_chr, cmprsn_var_nm_1L_chr))) %>%
        dplyr::rename_with(.cols = tidyselect::all_of(selected_cols_chr),
                           ~ ready4fun::get_from_lup_obj(rename_lup,
                                                         match_var_nm_1L_chr = "old_name_chr",
                                                         match_value_xx = .x,
                                                         target_var_nm_1L_chr = "new_name_chr",
                                                         evaluate_lgl = F)) %>%
        dplyr::summarise(dplyr::across(.fns = mean)) %>%
        dplyr::ungroup()  %>%
        tidyr::pivot_wider(names_from = tidyselect::all_of(cmprsn_var_nm_1L_chr),
                           values_from = rename_lup$new_name_chr) %>%
        dplyr::mutate(!!rlang::sym(paste0("diff_",costs_var_nm_1L_chr,"_dbl")) := !!rlang::sym(new_nms_ls[[1]][1]) - !!rlang::sym(new_nms_ls[[1]][2]),
                      !!rlang::sym(paste0("diff_",benefits_var_nm_1L_chr,"_dbl")) :=  (!!rlang::sym(new_nms_ls[[2]][1]) - !!rlang::sym(new_nms_ls[[2]][2])),
                      ICER = (!!rlang::sym(paste0("diff_",costs_var_nm_1L_chr,"_dbl")) /
                                  !!rlang::sym(paste0("diff_",benefits_var_nm_1L_chr,"_dbl")))) %>%
        dplyr::rename_with(~paste0(.x,"_dbl")) %>%
        t()
    summary_dbl <- summary_tb[,1]
    return(summary_dbl)
}
make_costs_vec_from_gamma_dist <- function(n_int,
                                           costs_mean_dbl,
                                           costs_sd_dbl){
    scale_1L_dbl <- costs_sd_dbl^2/costs_mean_dbl
    shape_1L_dbl <- costs_mean_dbl / scale_1L_dbl
    costs_dbl <- rgamma(n_int,shape = shape_1L_dbl, scale = scale_1L_dbl)
    return(costs_dbl)
}
make_fake_trial_ds <- function(ds_tb,
                               id_var_nm_1L_chr = "fkClientID",
                               round_var_nm_1L_chr = "round",
                               round_lvls_chr = c("Baseline","Follow-up"),
                               match_on_vars_chr,
                               cmprsn_var_nm_1L_chr = "study_arm_chr",
                               cmprsn_groups_chr = c("Intervention","Control"),
                               fns_ls,
                               var_nms_chr,
                               abs_mean_diff_dbl,
                               diff_sd_dbl,
                               multiplier_dbl,
                               min_dbl,
                               max_dbl,
                               integer_lgl,
                               match_idx_var_nm_1L_chr = "match_idx_int"){
    updated_ds_tb <- ds_tb %>%
        make_balanced_fake_ds(id_var_nm_1L_chr = id_var_nm_1L_chr,
                              round_var_nm_1L_chr = round_var_nm_1L_chr,
                              timepoint_bl_val_1L_chr = round_lvls_chr[1],
                              match_on_vars_chr = match_on_vars_chr,
                              cmprsn_var_nm_1L_chr = cmprsn_var_nm_1L_chr ,
                              cmprsn_groups_chr = cmprsn_groups_chr) %>%
        add_diffs_by_group_and_tmpt(cmprsn_var_nm_1L_chr = cmprsn_var_nm_1L_chr ,
                                    cmprsn_group_match_val_chr = cmprsn_groups_chr[1],
                                    round_var_nm_1L_chr = round_var_nm_1L_chr,
                                    timepoint_match_val_1L_chr = round_lvls_chr[2],
                                    match_idx_var_nm_1L_chr = match_idx_var_nm_1L_chr,
                                    var_nms_chr = var_nms_chr,
                                    fns_ls = fns_ls,
                                    abs_mean_diff_dbl = abs_mean_diff_dbl,
                                    diff_sd_dbl = diff_sd_dbl,
                                    multiplier_dbl = multiplier_dbl,
                                    min_dbl = min_dbl,
                                    max_dbl = max_dbl,
                                    integer_lgl = integer_lgl
                                    )
    return(updated_ds_tb)
}
make_formula <- function(dep_var_nm_1L_chr,
                         predictors_chr,
                         environment_env = parent.frame()){
    formula_fml <- formula(paste0(dep_var_nm_1L_chr,
                                  " ~ ",
                                  paste0(predictors_chr, collapse = " + ")), env = environment_env)
    return(formula_fml)
}
make_he_smry <- function(ds_tb,
                         change_vars_chr = NA_character_,
                         wtp_dbl = 50000,
                         bootstrap_iters_1L_int = 1000,
                         change_types_chr = "dbl",
                         benefits_pfx_1L_chr = "qalys_dbl",
                         benefits_var_nm_1L_chr = "qalys",
                         costs_pfx_1L_chr = "costs_dbl",
                         costs_var_nm_1L_chr = "costs",
                         change_sfx_1L_chr = "change",
                         cmprsn_groups_chr = c("Intervention","Control"),
                         cmprsn_var_nm_1L_chr = "study_arm_chr",
                         round_fup_1L_chr = "Follow-up"
){
    bootstraps_ls <- boot::boot(ds_tb,
                                make_ce_smry,
                                R = bootstrap_iters_1L_int,
                                benefits_pfx_1L_chr = benefits_pfx_1L_chr,
                                costs_pfx_1L_chr = costs_pfx_1L_chr,
                                change_vars_chr = change_vars_chr,
                                change_sfx_1L_chr = change_sfx_1L_chr,
                                change_types_chr = change_types_chr,
                                cmprsn_groups_chr = cmprsn_groups_chr,
                                cmprsn_var_nm_1L_chr = cmprsn_var_nm_1L_chr,
                                round_fup_1L_chr = round_fup_1L_chr,
                                benefits_var_nm_1L_chr = benefits_var_nm_1L_chr,
                                costs_var_nm_1L_chr = costs_var_nm_1L_chr
    ) # round_fup_1L_chr = "Follow-up", ds_smry_ls$round_lvls_chr[2]
    costs_mat <- bootstraps_ls$t[,1:2]
    benefits_mat <- bootstraps_ls$t[,3:4]
    reordered_cmprsns_chr <- cmprsn_groups_chr[cmprsn_groups_chr %>% purrr::map_int(~which(endsWith(names(bootstraps_ls$t0)[1:2],                                                                       paste0(.x,"_dbl"))))]
    ce_res_ls <- BCEA::bcea(e = benefits_mat,
                            c = costs_mat,
                            ref = which(reordered_cmprsns_chr == cmprsn_groups_chr[1]),
                            interventions = reordered_cmprsns_chr,
                            #wtp = NULL,
                            Kmax = wtp_dbl*2)
    named_mat <- bootstraps_ls$t
    colnames(named_mat) <- names(bootstraps_ls$t0)
    he_smry_ls <- list(ce_res_ls = ce_res_ls,
                       bootstraps_ls = bootstraps_ls,
                       benefits_mat = benefits_mat,
                       costs_mat = costs_mat,
                       named_mat = named_mat
    )
    return(he_smry_ls)
}
make_matched_ds <- function(sngl_grp_ds,
                            cmprsn_smry_tb,
                            ds_smry_ls){
    matched_ds_tb <- sngl_grp_ds_tb %>%
        make_fake_trial_ds(id_var_nm_1L_chr = ds_smry_ls$id_var_nm_1L_chr,
                           round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr,
                           round_lvls_chr = ds_smry_ls$round_lvls_chr,
                           match_on_vars_chr = cmprsn_smry_tb$var_nms_chr,
                           cmprsn_var_nm_1L_chr = ds_smry_ls$cmprsn_var_nm_1L_chr ,
                           cmprsn_groups_chr = ds_smry_ls$cmprsn_groups_chr,
                           fns_ls = cmprsn_smry_tb$fns_ls,
                           var_nms_chr = cmprsn_smry_tb$var_nms_chr,
                           abs_mean_diff_dbl = cmprsn_smry_tb$abs_mean_diff_dbl,
                           diff_sd_dbl = cmprsn_smry_tb$diff_sd_dbl,
                           multiplier_dbl = cmprsn_smry_tb$multiplier_dbl,
                           min_dbl = cmprsn_smry_tb$min_dbl,
                           max_dbl = cmprsn_smry_tb$max_dbl,
                           integer_lgl = cmprsn_smry_tb$integer_lgl)
    if(any(cmprsn_smry_tb$integer_lgl))
        matched_ds_tb <-   matched_ds_tb  %>%
            dplyr::mutate(dplyr::across(cmprsn_smry_tb$var_nms_chr[cmprsn_smry_tb$integer_lgl],as.integer))
    return(matched_ds_tb)
}
make_matched_ds_spine <- function(ds_tb,
                            round_var_nm_1L_chr = "Timepoint_chr",
                            timepoint_bl_val_1L_chr = "Baseline",
                            cmprsn_var_nm_1L_chr = "study_arm_chr",
                            active_arm_val_1L_chr = "Intervention",
                            id_var_nm_1L_chr = "fkClientID",
                            match_on_vars_chr){
    match_ds <- ds_tb %>% dplyr::filter(!!rlang::sym(round_var_nm_1L_chr) == timepoint_bl_val_1L_chr) %>%
        dplyr::mutate(Intervention_lgl = dplyr::case_when(!!rlang::sym(cmprsn_var_nm_1L_chr) == active_arm_val_1L_chr ~ T,
                                                          T ~ F))
    matched_ls <- make_formula("Intervention_lgl",
                               predictors_chr = match_on_vars_chr) %>%
        MatchIt::matchit(data = match_ds, method = "nearest", ratio = 1)
    matched_ds_tb <- MatchIt::match.data(matched_ls)
    match_key_ds_tb <- c(F,T) %>% purrr::map_dfr(~matched_ds_tb %>% dplyr::filter(Intervention_lgl == .x) %>% dplyr::arrange(distance) %>% dplyr::mutate(match_idx_int = 1:dplyr::n())) %>% dplyr::arrange(match_idx_int) %>% dplyr::select(!!rlang::sym(id_var_nm_1L_chr), study_arm_chr, match_idx_int)
    matched_ds_tb <- dplyr::right_join(ds_tb,match_key_ds_tb) %>% dplyr::arrange(match_idx_int)
    return(matched_ds_tb)
}
make_sngl_grp_ds <- function(seed_ds_tb = NULL,
                             ds_smry_ls){
    sngl_grp_ds_tb <- seed_ds_tb %>%
        dplyr::select(ds_smry_ls$id_var_nm_1L_chr,
                      ds_smry_ls$round_var_nm_1L_chr,
                      ds_smry_ls$predr_var_nms) %>% #
        na.omit()
    if("SOFAS" %in% ds_smry_ls$predr_var_nms) # TEMPORARY
        sngl_grp_ds_tb <- sngl_grp_ds_tb %>%
            dplyr::mutate(SOFAS = as.integer(round(SOFAS,0)))
    sngl_grp_ds_tb <- sngl_grp_ds_tb %>%
        tibble::as_tibble() %>%
        add_dates_from_dist(bl_start_date_dtm = ds_smry_ls$bl_start_date_dtm,
                            bl_end_date_dtm = ds_smry_ls$bl_end_date_dtm,
                            duration_args_ls = ds_smry_ls$duration_args_ls,
                            duration_fn = ds_smry_ls$duration_fn,
                            date_var_nm_1L_chr = ds_smry_ls$date_var_nm_1L_chr,
                            id_var_nm_1L_chr = ds_smry_ls$id_var_nm_1L_chr,
                            round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr,
                            round_bl_val_1L_chr = ds_smry_ls$round_lvls_chr[1])
    if(!is.null(ds_smry_ls$costs_mean_dbl) & !is.null(ds_smry_ls$costs_sd_dbl))
        sngl_grp_ds_tb <- sngl_grp_ds_tb  %>%
        add_costs_by_tmpt(round_var_nm_1L_chr = ds_smry_ls$round_var_nm_1L_chr,
                          round_lvls_chr = ds_smry_ls$round_lvls_chr,
                          costs_mean_dbl = ds_smry_ls$costs_mean_dbl,
                          costs_sd_dbl = ds_smry_ls$costs_sd_dbl,
                          extra_cost_args_ls = list(costs_var_nm_1L_chr = ds_smry_ls$costs_var_nm_1L_chr),
                          fn = add_costs_from_gamma_dist)
    sngl_grp_ds_tb <- sngl_grp_ds_tb %>% dplyr::arrange(ds_smry_ls$id_var_nm_1L_chr)
    return(sngl_grp_ds_tb)
}

