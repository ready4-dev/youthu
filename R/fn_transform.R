#' Transform dataset for tstng
#' @description transform_ds_for_tstng() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform dataset for tstng. Function argument data_tb specifies the object to be updated. Argument dep_var_nm_1L_chr provides the object to be updated. The function returns Transformed data (a tibble).
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param dep_var_max_val_1L_dbl Dep var max value (a double vector of length one), Default: 0.999
#' @param candidate_predrs_chr Candidate predrs (a character vector), Default: 'NA'
#' @param covar_var_nms_chr Covar var names (a character vector), Default: 'NA'
#' @param round_var_nm_1L_chr Round var name (a character vector of length one), Default: 'round'
#' @param round_val_1L_chr Round value (a character vector of length one), Default: 'Baseline'
#' @param remove_all_mssng_1L_lgl Remove all mssng (a logical vector of length one), Default: F
#' @return Transformed data (a tibble)
#' @rdname transform_ds_for_tstng
#' @export 
#' @importFrom purrr discard
#' @importFrom dplyr filter select mutate
#' @importFrom rlang sym syms
#' @importFrom stats na.omit
transform_ds_for_tstng <- function (data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", dep_var_max_val_1L_dbl = 0.999, 
    candidate_predrs_chr = NA_character_, covar_var_nms_chr = NA_character_, 
    round_var_nm_1L_chr = "round", round_val_1L_chr = "Baseline", 
    remove_all_mssng_1L_lgl = F) 
{
    vars_to_keep_chr <- c(dep_var_nm_1L_chr, candidate_predrs_chr, 
        covar_var_nms_chr) %>% purrr::discard(is.na)
    tfd_data_tb <- data_tb %>% dplyr::filter(!!rlang::sym(round_var_nm_1L_chr) == 
        round_val_1L_chr) %>% dplyr::select(!!!rlang::syms(vars_to_keep_chr)) %>% 
        dplyr::mutate(`:=`(!!rlang::sym(dep_var_nm_1L_chr), ifelse(!!rlang::sym(dep_var_nm_1L_chr) == 
            1, 0.999, !!rlang::sym(dep_var_nm_1L_chr))))
    if (remove_all_mssng_1L_lgl) 
        tfd_data_tb <- tfd_data_tb %>% stats::na.omit()
    return(tfd_data_tb)
}
#' Transform raw Assessment of Quality of Life tibble to Assessment of Quality of Life Six Dimension
#' @description transform_raw_aqol_tb_to_aqol6d_tb() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform raw assessment of quality of life tibble to assessment of quality of life six dimension tibble. Function argument raw_aqol_tb specifies the object to be updated. The function returns Assessment of Quality of Life Six Dimension (a tibble).
#' @param raw_aqol_tb Raw Assessment of Quality of Life (a tibble)
#' @return Assessment of Quality of Life Six Dimension (a tibble)
#' @rdname transform_raw_aqol_tb_to_aqol6d_tb
#' @export 
#' @importFrom dplyr mutate select contains rename
transform_raw_aqol_tb_to_aqol6d_tb <- function (raw_aqol_tb) 
{
    aqol6d_tb <- raw_aqol_tb %>% dplyr::mutate(d_agegroup = cut(d_age, 
        breaks = c(11, 17, 30), labels = c("Age 12-17", "Age 18-26"))) %>% 
        dplyr::mutate(round = factor(round, labels = c("Baseline", 
            "Follow-up"))) %>% dplyr::select(fkClientID, c_p_diag_s, 
        s_centre, c_clinical_staging_s, d_age, d_agegroup, d_gender, 
        d_sex_birth_s, d_sexual_ori_s, d_country_bir_s, d_ATSI, 
        d_english_home, d_english_native, d_relation_s, d_studying_working, 
        k6_total, phq9_total, bads_total, gad7_total, oasis_total, 
        scared_total, dplyr::contains("aqol6d"), c_sofas, round) %>% 
        dplyr::mutate(Gender = factor(ifelse(d_gender == "Genderqueer/gender nonconforming/agender" | 
            d_gender == "Transgender", "Other", as.character(d_gender)))) %>% 
        dplyr::mutate(Region = as.factor(ifelse(s_centre == "Canberra" | 
            s_centre == "Southport" | s_centre == "Knox", "Metro", 
            "Regional"))) %>% dplyr::mutate(CALD = factor(ifelse(d_country_bir_s == 
        "Other" | d_english_home == "No" | d_english_native == 
        "No" | d_ATSI == "Yes", "Yes", "No"))) %>% dplyr::rename(PHQ9 = phq9_total, 
        BADS = bads_total, GAD7 = gad7_total, OASIS = oasis_total, 
        SCARED = scared_total, K6 = k6_total, SOFAS = c_sofas)
    return(aqol6d_tb)
}
#' Transform ts mdl data
#' @description transform_ts_mdl_data() is a Transform function that edits an object in such a way that core object attributes - e.g. shape, dimensions, elements, type - are altered. Specifically, this function implements an algorithm to transform ts mdl data. Function argument mdl_ls specifies the object to be updated. Argument data_tb provides the object to be updated. The function returns Cnfdl mdl (a list).
#' @param mdl_ls Mdl (a list)
#' @param data_tb Data (a tibble)
#' @param dep_var_nm_1L_chr Dep var name (a character vector of length one), Default: 'aqol6d_total_w'
#' @param predr_vars_nms_chr Predr vars names (a character vector)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param mdl_nm_1L_chr Mdl name (a character vector of length one)
#' @return Cnfdl mdl (a list)
#' @rdname transform_ts_mdl_data
#' @export 
#' @importFrom dplyr select all_of summarise across everything
#' @importFrom purrr map flatten_chr
transform_ts_mdl_data <- function (mdl_ls, data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", 
    predr_vars_nms_chr, id_var_nm_1L_chr = "fkClientID", mdl_nm_1L_chr) 
{
    old_data_tb <- data_tb %>% dplyr::select(c(dplyr::all_of(id_var_nm_1L_chr), 
        dplyr::all_of(dep_var_nm_1L_chr), predr_vars_nms_chr %>% 
            purrr::map(~paste0(.x, c("", "_baseline", "_change"))) %>% 
            purrr::flatten_chr()))
    cnfdl_mdl_ls <- mdl_ls
    cnfdl_mdl_ls$data <- old_data_tb %>% as.data.frame() %>% 
        dplyr::summarise(dplyr::across(dplyr::everything(), ~sample(.x, 
            1)))
    return(cnfdl_mdl_ls)
}
