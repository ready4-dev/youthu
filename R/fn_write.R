#' Write results to comma separated variables file
#' @description write_results_to_csv() is a Write function that writes a file to a specified local directory. Specifically, this function implements an algorithm to write results to comma separated variables file. The function returns Datasets (a tibble).
#' @param synth_data_spine_ls Synth data spine (a list)
#' @param output_dir_1L_chr Output directory (a character vector of length one), Default: '.'
#' @return Datasets (a tibble)
#' @rdname write_results_to_csv
#' @export 
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map_dfr map_dfc map map_dbl walk2
#' @importFrom stats setNames
#' @importFrom dplyr mutate select everything
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
