make_adol_aqol6d_disv_lup <- function ()
{
    adol_aqol6d_disv_lup <- aqol6d_adult_disv_lup_tb %>% dplyr::mutate(Answer_4_dbl = dplyr::case_when(Question_chr ==
        "Q18" ~ 0.622, TRUE ~ Answer_4_dbl), Answer_5_dbl = dplyr::case_when(Question_chr ==
        "Q3" ~ 0.827, TRUE ~ Answer_5_dbl), Answer_6_dbl = dplyr::case_when(Question_chr ==
        "Q1" ~ 0.073, TRUE ~ Answer_5_dbl))
    return(adol_aqol6d_disv_lup)
}
make_aqol6d_adol_pop_tbs_ls <- function (aqol_items_props_tbs_ls, aqol_scores_pars_ls, series_names_chr,
    synth_data_spine_ls, temporal_cors_ls, id_var_nm_1L_chr = "fkClientID",
    prefix_chr = c(uid = "Participant_", aqol_item = "aqol6d_q",
        domain_unwtd_pfx_1L_chr = "aqol6d_subtotal_c_", domain_wtd_pfx_1L_chr = "aqol6d_subtotal_w_"))
{
    item_pfx_1L_chr <- prefix_chr[["aqol_item"]]
    uid_pfx_1L_chr <- prefix_chr[["uid"]]
    aqol6d_adol_pop_tbs_ls <- make_synth_series_tbs_ls(synth_data_spine_ls,
        series_names_chr = series_names_chr) %>% add_cors_and_uts_to_aqol6d_tbs_ls(aqol_scores_pars_ls = aqol_scores_pars_ls,
        aqol_items_props_tbs_ls = aqol_items_props_tbs_ls, temporal_cors_ls = temporal_cors_ls,
        prefix_chr = prefix_chr, aqol_tots_var_nms_chr = synth_data_spine_ls$aqol_tots_var_nms_chr,
        id_var_nm_1L_chr = id_var_nm_1L_chr) %>% purrr::map(~{
        domain_items_ls <- make_domain_items_ls(domain_qs_lup_tb = aqol6d_domain_qs_lup_tb,
            item_pfx_1L_chr = item_pfx_1L_chr)
        domain_items_ls %>% add_unwtd_dim_tots(items_tb = .x,
            domain_pfx_1L_chr = prefix_chr[["domain_unwtd_pfx_1L_chr"]]) %>%
            add_wtd_dim_tots(domain_items_ls = domain_items_ls,
                domain_unwtd_pfx_1L_chr = prefix_chr[["domain_unwtd_pfx_1L_chr"]],
                domain_wtd_pfx_1L_chr = prefix_chr[["domain_wtd_pfx_1L_chr"]]) %>%
            add_labels_to_aqol6d_tb()
    }) %>% purrr::map(~.x %>% dplyr::select(!!rlang::sym(id_var_nm_1L_chr),
        dplyr::starts_with(item_pfx_1L_chr), dplyr::starts_with(prefix_chr[["domain_unwtd_pfx_1L_chr"]]),
        dplyr::starts_with(prefix_chr[["domain_wtd_pfx_1L_chr"]]),
        dplyr::everything()))
    return(aqol6d_adol_pop_tbs_ls)
}
make_aqol6d_fns_ls <- function (domain_items_ls)
{
    aqol6d_disu_fn_ls <- paste0("calculate_aqol6d_dim_", 1:length(domain_items_ls),
        "_disv") %>% purrr::map(~rlang::sym(.x))
    return(aqol6d_disu_fn_ls)
}
make_aqol6d_items_tb <- function (aqol_tb, old_pfx_1L_chr, new_pfx_1L_chr)
{
    aqol6d_items_tb <- aqol_tb %>% dplyr::select(dplyr::starts_with(old_pfx_1L_chr)) %>%
        dplyr::rename_all(~{
            stringr::str_replace(., old_pfx_1L_chr, new_pfx_1L_chr)
        })
    return(aqol6d_items_tb)
}
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
make_complete_props_tbs_ls <- function (raw_props_tbs_ls, question_var_nm_1L_chr = "Question")
{
    complete_props_tbs_ls <- raw_props_tbs_ls %>% purrr::map(~{
        .x %>% dplyr::mutate(total_prop_dbl = rowSums(dplyr::select(.,
            -!!rlang::sym(question_var_nm_1L_chr)), na.rm = T) -
            100) %>% dplyr::mutate_if(is.numeric, ~purrr::map2_dbl(.,
            total_prop_dbl, ~ifelse(.x == 100, 1 - .y, .x))) %>%
            dplyr::select(-total_prop_dbl)
    })
    return(complete_props_tbs_ls)
}
make_correlated_data_tb <- function (synth_data_spine_ls, synth_data_idx_1L_dbl = 1)
{
    correlated_data_tb <- simstudy::genCorData(synth_data_spine_ls$nbr_obs_dbl[synth_data_idx_1L_dbl],
        mu = synth_data_spine_ls$means_ls[[synth_data_idx_1L_dbl]],
        sigma = synth_data_spine_ls$sds_ls[[synth_data_idx_1L_dbl]],
        corMatrix = make_pdef_cor_mat_mat(synth_data_spine_ls$cor_mat_ls[[synth_data_idx_1L_dbl]]),
        cnames = synth_data_spine_ls$var_names_chr) %>% force_min_max_and_int_cnstrs(var_names_chr = synth_data_spine_ls$var_names_chr,
        min_max_ls = synth_data_spine_ls$min_max_ls, discrete_lgl = synth_data_spine_ls$discrete_lgl)
    return(correlated_data_tb)
}
make_costs_vec_from_gamma_dist <- function(n_int,
                                           costs_mean_dbl,
                                           costs_sd_dbl){
    scale_1L_dbl <- costs_sd_dbl^2/costs_mean_dbl
    shape_1L_dbl <- costs_mean_dbl / scale_1L_dbl
    costs_dbl <- rgamma(n_int,shape = shape_1L_dbl, scale = scale_1L_dbl)
    return(costs_dbl)
}
make_dim_sclg_cons_dbl <- function (domains_chr, dim_sclg_con_lup_tb)
{
    dim_sclg_cons_dbl <- purrr::map_dbl(domains_chr, ~ready4fun::get_from_lup_obj(dim_sclg_con_lup_tb,
        match_var_nm_1L_chr = "Dimension_chr", match_value_xx = .x,
        target_var_nm_1L_chr = "Constant_dbl", evaluate_lgl = F))
    return(dim_sclg_cons_dbl)
}
make_domain_items_ls <- function (domain_qs_lup_tb, item_pfx_1L_chr)
{
    domains_chr <- domain_qs_lup_tb$Domain_chr %>% unique()
    q_nbrs_ls <- purrr::map(domains_chr, ~domain_qs_lup_tb %>%
        dplyr::filter(Domain_chr == .x) %>% dplyr::pull(Question_dbl))
    domain_items_ls <- purrr::map(q_nbrs_ls, ~paste0(item_pfx_1L_chr,
        .x)) %>% stats::setNames(domains_chr)
    return(domain_items_ls)
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
make_item_wrst_wghts_ls_ls <- function (domain_items_ls, itm_wrst_wghts_lup_tb)
{
    item_wrst_wghts_ls_ls <- domain_items_ls %>% purrr::map(~{
        purrr::map_dbl(.x, ~{
            ready4fun::get_from_lup_obj(itm_wrst_wghts_lup_tb,
                                        match_var_nm_1L_chr = "Question_chr", match_value_xx = .x,
                                        target_var_nm_1L_chr = "Worst_Weight_dbl", evaluate_lgl = F)
        })
    })
    return(item_wrst_wghts_ls_ls)
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
make_pdef_cor_mat_mat <- function (lower_diag_mat)
{
    pdef_cor_mat <- lower_diag_mat %>% Matrix::forceSymmetric(uplo = "L") %>%
        as.matrix()
    if (!matrixcalc::is.positive.definite(pdef_cor_mat)) {
        pdef_cor_mat <- psych::cor.smooth(pdef_cor_mat)
    }
    return(pdef_cor_mat)
}
make_sngl_grp_ds <- function(seed_ds_tb = NULL,
                             ds_smry_ls){
    sngl_grp_ds_tb <- seed_ds_tb %>%
        dplyr::select(ds_smry_ls$id_var_nm_1L_chr,
                      ds_smry_ls$round_var_nm_1L_chr,
                      ds_smry_ls$predr_var_nms) %>% #
        na.omit()
    if("SOFAS" %in% ds_smry_ls$round_var_nm_1L_chr) # TEMPORARY
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
make_synth_series_tbs_ls <- function (synth_data_spine_ls, series_names_chr)
{
    synth_series_tbs_ls <- 1:length(series_names_chr) %>% purrr::map(~make_correlated_data_tb(synth_data_spine_ls = synth_data_spine_ls,
        synth_data_idx_1L_dbl = .x) %>% replace_with_missing_vals(synth_data_spine_ls = synth_data_spine_ls,
        idx_int = .x)) %>% stats::setNames(series_names_chr)
    return(synth_series_tbs_ls)
}
make_vec_with_sum_of_int <- function (target_int, start_int, end_int, length_int)
{
    vec_int <- Surrogate::RandVec(a = start_int, b = end_int,
        s = target_int, n = length_int, m = 1) %>% purrr::pluck("RandVecOutput") %>%
        as.vector() %>% round() %>% as.integer() %>% force_vec_to_sum_to_int(target_1L_int = target_int)
    return(vec_int)
}
