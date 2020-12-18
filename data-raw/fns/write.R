write_box_cox_tfmn <- function (data_tb, predr_var_nm_1L_chr, path_to_write_to_1L_chr,
    dep_var_nm_1L_chr = "aqol6d_total_w", covar_var_nms_chr = NA_character_,
    fl_nm_pfx_1L_chr = "A_RT", height_1L_dbl = 6, width_1L_dbl = 6,
    start_1L_chr = NULL, mdl_types_lup = NULL)
{
    if (is.null(mdl_types_lup))
        utils::data("mdl_types_lup", envir = environment())
    mdl <- make_mdl(data_tb, dep_var_nm_1L_chr = dep_var_nm_1L_chr,
        predr_var_nm_1L_chr = predr_var_nm_1L_chr, covar_var_nms_chr = covar_var_nms_chr,
        mdl_type_1L_chr = "OLS_NTF", mdl_types_lup = mdl_types_lup,
        start_1L_chr = start_1L_chr)
    path_to_plot_1L_chr <- write_brm_mdl_plt_fl(plt_fn = MASS::boxcox,
        fn_args_ls = list(mdl, plotit = T), path_to_write_to_1L_chr = path_to_write_to_1L_chr,
        plt_nm_1L_chr = paste0(fl_nm_pfx_1L_chr, "_", predr_var_nm_1L_chr,
            "_", "BOXCOX"), height_1L_dbl = height_1L_dbl, width_1L_dbl = width_1L_dbl)
    return(path_to_plot_1L_chr)
}
write_results_to_csv <- function (synth_data_spine_ls, output_dir_1L_chr = ".")
{
    measurements_tb <- tibble::tibble(timepoint_nms_chr = synth_data_spine_ls$timepoint_nms_chr,
        nbr_obs_dbl = synth_data_spine_ls$nbr_obs_dbl)
    var_summ_res_tb <- suppressMessages(purrr::map_dfr(1:length(synth_data_spine_ls$timepoint_nms_chr),
        ~{
            idx_dbl <- .x
            suppressWarnings({
                synth_data_spine_ls[c(4:6)] %>% purrr::map_dfc(~.x[idx_dbl])
            }) %>% stats::setNames(c("Mean", "SD", "N_Missing")) %>%
                dplyr::mutate(var_names_chr = synth_data_spine_ls$var_names_chr,
                  timepoint_nms_chr = synth_data_spine_ls$timepoint_nms_chr[idx_dbl]) %>%
                dplyr::select(timepoint_nms_chr, var_names_chr,
                  dplyr::everything())
        }))
    cor_tb_ls <- synth_data_spine_ls$cor_mat_ls %>% purrr::map(~tibble::as_tibble(.x) %>%
        stats::setNames(synth_data_spine_ls$var_names_chr) %>%
        dplyr::mutate(var_names_chr = synth_data_spine_ls$var_names_chr) %>%
        dplyr::select(var_names_chr, dplyr::everything())) %>%
        stats::setNames(paste0(synth_data_spine_ls$timepoint_nms_chr,
            "_correlations_tb"))
    var_class_pars_tb <- synth_data_spine_ls[7:9] %>% tibble::as_tibble() %>%
        dplyr::mutate(min_dbl = purrr::map_dbl(min_max_ls, ~.x[1]),
            max_dbl = purrr::map_dbl(min_max_ls, ~.x[2])) %>%
        dplyr::select(var_names_chr, dplyr::everything(), -min_max_ls)
    output_ls <- list(measurements_tb = measurements_tb, var_summ_res_tb = var_summ_res_tb,
        var_class_pars_tb = var_class_pars_tb) %>% append(cor_tb_ls)
    dss_tb <- tibble::tibble(ds_obj_nm_chr = names(output_ls),
        title_chr = c("Brief summary table of the number of observations for which data was collected at each study timepoint.",
            "Summary statistics (Mean, SD and Number Missing) for AQoL6D health utility and six mental health outcome measures for each study timepoint.",
            "Brief information about the data structure (whether discrete and allowable range) of AQoL6D health utility and six mental health outcome variables.",
            paste0("Correlation matrix for AQoL6D health utility and six mental health outcome measures at the ",
                synth_data_spine_ls$timepoint_nms_chr, " study timepoint.")))
    purrr::walk2(output_ls, names(output_ls), ~write.csv(.x,
        file = paste0(output_dir_1L_chr, "/", .y, ".csv"), row.names = F))
    return(dss_tb)
}

