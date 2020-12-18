#' Make adolescent Assessment of Quality of Life Six Dimension disvalue
#' @description make_adol_aqol6d_disv_lup() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make adolescent assessment of quality of life six dimension disvalue lookup table. The function returns Adolescent Assessment of Quality of Life Six Dimension disvalue (a lookup table).

#' @return Adolescent Assessment of Quality of Life Six Dimension disvalue (a lookup table)
#' @rdname make_adol_aqol6d_disv_lup
#' @export 
#' @importFrom dplyr mutate case_when
#' @keywords internal
make_adol_aqol6d_disv_lup <- function () 
{
    adol_aqol6d_disv_lup <- aqol6d_adult_disv_lup_tb %>% dplyr::mutate(Answer_4_dbl = dplyr::case_when(Question_chr == 
        "Q18" ~ 0.622, TRUE ~ Answer_4_dbl), Answer_5_dbl = dplyr::case_when(Question_chr == 
        "Q3" ~ 0.827, TRUE ~ Answer_5_dbl), Answer_6_dbl = dplyr::case_when(Question_chr == 
        "Q1" ~ 0.073, TRUE ~ Answer_5_dbl))
    return(adol_aqol6d_disv_lup)
}
#' Make Assessment of Quality of Life Six Dimension adolescent pop tibbles
#' @description make_aqol6d_adol_pop_tbs_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make assessment of quality of life six dimension adolescent pop tibbles list. The function returns Assessment of Quality of Life Six Dimension adolescent pop tibbles (a list).
#' @param aqol_items_props_tbs_ls Assessment of Quality of Life items props tibbles (a list)
#' @param aqol_scores_pars_ls Assessment of Quality of Life scores parameters (a list)
#' @param series_names_chr Series names (a character vector)
#' @param synth_data_spine_ls Synth data spine (a list)
#' @param temporal_cors_ls Temporal cors (a list)
#' @param id_var_nm_1L_chr Id var name (a character vector of length one), Default: 'fkClientID'
#' @param prefix_chr Prefix (a character vector), Default: c(uid = "Participant_", aqol_item = "aqol6d_q", domain_unwtd_pfx_1L_chr = "aqol6d_subtotal_c_", 
#'    domain_wtd_pfx_1L_chr = "aqol6d_subtotal_w_")
#' @return Assessment of Quality of Life Six Dimension adolescent pop tibbles (a list)
#' @rdname make_aqol6d_adol_pop_tbs_ls
#' @export 
#' @importFrom purrr map
#' @importFrom dplyr select starts_with everything
#' @importFrom rlang sym
#' @keywords internal
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
#' Make Assessment of Quality of Life Six Dimension functions
#' @description make_aqol6d_fns_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make assessment of quality of life six dimension functions list. The function returns Assessment of Quality of Life Six Dimension disu (a list of functions).
#' @param domain_items_ls Domain items (a list)
#' @return Assessment of Quality of Life Six Dimension disu (a list of functions)
#' @rdname make_aqol6d_fns_ls
#' @export 
#' @importFrom purrr map
#' @importFrom rlang sym
make_aqol6d_fns_ls <- function (domain_items_ls) 
{
    aqol6d_disu_fn_ls <- paste0("calculate_aqol6d_dim_", 1:length(domain_items_ls), 
        "_disv") %>% purrr::map(~rlang::sym(.x))
    return(aqol6d_disu_fn_ls)
}
#' Make Assessment of Quality of Life Six Dimension items
#' @description make_aqol6d_items_tb() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make assessment of quality of life six dimension items tibble. The function returns Assessment of Quality of Life Six Dimension items (a tibble).
#' @param aqol_tb Assessment of Quality of Life (a tibble)
#' @param old_pfx_1L_chr Old prefix (a character vector of length one)
#' @param new_pfx_1L_chr New prefix (a character vector of length one)
#' @return Assessment of Quality of Life Six Dimension items (a tibble)
#' @rdname make_aqol6d_items_tb
#' @export 
#' @importFrom dplyr select starts_with rename_all
#' @importFrom stringr str_replace
make_aqol6d_items_tb <- function (aqol_tb, old_pfx_1L_chr, new_pfx_1L_chr) 
{
    aqol6d_items_tb <- aqol_tb %>% dplyr::select(dplyr::starts_with(old_pfx_1L_chr)) %>% 
        dplyr::rename_all(~{
            stringr::str_replace(., old_pfx_1L_chr, new_pfx_1L_chr)
        })
    return(aqol6d_items_tb)
}
#' Make complete props tibbles
#' @description make_complete_props_tbs_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make complete props tibbles list. The function returns Complete props tibbles (a list).
#' @param raw_props_tbs_ls Raw props tibbles (a list)
#' @param question_var_nm_1L_chr Question var name (a character vector of length one), Default: 'Question'
#' @return Complete props tibbles (a list)
#' @rdname make_complete_props_tbs_ls
#' @export 
#' @importFrom purrr map map2_dbl
#' @importFrom dplyr mutate select mutate_if
#' @importFrom rlang sym
#' @keywords internal
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
#' Make correlated data
#' @description make_correlated_data_tb() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make correlated data tibble. The function returns Correlated data (a tibble).
#' @param synth_data_spine_ls Synth data spine (a list)
#' @param synth_data_idx_1L_dbl Synth data index (a double vector of length one), Default: 1
#' @return Correlated data (a tibble)
#' @rdname make_correlated_data_tb
#' @export 
#' @importFrom simstudy genCorData
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
#' Make dimension sclg cons
#' @description make_dim_sclg_cons_dbl() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make dimension sclg cons double vector. The function returns Dimension sclg cons (a double vector).
#' @param domains_chr Domains (a character vector)
#' @param dim_sclg_con_lup_tb Dimension sclg con lookup table (a tibble)
#' @return Dimension sclg cons (a double vector)
#' @rdname make_dim_sclg_cons_dbl
#' @export 
#' @importFrom purrr map_dbl
#' @importFrom ready4fun get_from_lup_obj
make_dim_sclg_cons_dbl <- function (domains_chr, dim_sclg_con_lup_tb) 
{
    dim_sclg_cons_dbl <- purrr::map_dbl(domains_chr, ~ready4fun::get_from_lup_obj(dim_sclg_con_lup_tb, 
        match_var_nm_1L_chr = "Dimension_chr", match_value_xx = .x, 
        target_var_nm_1L_chr = "Constant_dbl", evaluate_lgl = F))
    return(dim_sclg_cons_dbl)
}
#' Make domain items
#' @description make_domain_items_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make domain items list. The function returns Domain items (a list).
#' @param domain_qs_lup_tb Domain questions lookup table (a tibble)
#' @param item_pfx_1L_chr Item prefix (a character vector of length one)
#' @return Domain items (a list)
#' @rdname make_domain_items_ls
#' @export 
#' @importFrom purrr map
#' @importFrom dplyr filter pull
#' @importFrom stats setNames
make_domain_items_ls <- function (domain_qs_lup_tb, item_pfx_1L_chr) 
{
    domains_chr <- domain_qs_lup_tb$Domain_chr %>% unique()
    q_nbrs_ls <- purrr::map(domains_chr, ~domain_qs_lup_tb %>% 
        dplyr::filter(Domain_chr == .x) %>% dplyr::pull(Question_dbl))
    domain_items_ls <- purrr::map(q_nbrs_ls, ~paste0(item_pfx_1L_chr, 
        .x)) %>% stats::setNames(domains_chr)
    return(domain_items_ls)
}
#' Make item wrst wghts
#' @description make_item_wrst_wghts_ls_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make item wrst wghts list list. The function returns Item wrst wghts (a list of lists).
#' @param domain_items_ls Domain items (a list)
#' @param itm_wrst_wghts_lup_tb Itm wrst wghts lookup table (a tibble)
#' @return Item wrst wghts (a list of lists)
#' @rdname make_item_wrst_wghts_ls_ls
#' @export 
#' @importFrom purrr map map_dbl
#' @importFrom ready4fun get_from_lup_obj
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
#' Make pdef cor matrix
#' @description make_pdef_cor_mat_mat() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make pdef cor matrix matrix. The function returns Pdef cor (a matrix).
#' @param lower_diag_mat Lower diag (a matrix)
#' @return Pdef cor (a matrix)
#' @rdname make_pdef_cor_mat_mat
#' @export 
#' @importFrom Matrix forceSymmetric
#' @importFrom matrixcalc is.positive.definite
#' @importFrom psych cor.smooth
make_pdef_cor_mat_mat <- function (lower_diag_mat) 
{
    pdef_cor_mat <- lower_diag_mat %>% Matrix::forceSymmetric(uplo = "L") %>% 
        as.matrix()
    if (!matrixcalc::is.positive.definite(pdef_cor_mat)) {
        pdef_cor_mat <- psych::cor.smooth(pdef_cor_mat)
    }
    return(pdef_cor_mat)
}
#' Make synth series tibbles
#' @description make_synth_series_tbs_ls() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make synth series tibbles list. The function returns Synth series tibbles (a list).
#' @param synth_data_spine_ls Synth data spine (a list)
#' @param series_names_chr Series names (a character vector)
#' @return Synth series tibbles (a list)
#' @rdname make_synth_series_tbs_ls
#' @export 
#' @importFrom purrr map
#' @importFrom stats setNames
make_synth_series_tbs_ls <- function (synth_data_spine_ls, series_names_chr) 
{
    synth_series_tbs_ls <- 1:length(series_names_chr) %>% purrr::map(~make_correlated_data_tb(synth_data_spine_ls = synth_data_spine_ls, 
        synth_data_idx_1L_dbl = .x) %>% replace_with_missing_vals(synth_data_spine_ls = synth_data_spine_ls, 
        idx_int = .x)) %>% stats::setNames(series_names_chr)
    return(synth_series_tbs_ls)
}
#' Make vec with sum of
#' @description make_vec_with_sum_of_int() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make vec with sum of integer vector. The function returns Vec (an integer vector).
#' @param target_int Target (an integer vector)
#' @param start_int Start (an integer vector)
#' @param end_int End (an integer vector)
#' @param length_int Length (an integer vector)
#' @return Vec (an integer vector)
#' @rdname make_vec_with_sum_of_int
#' @export 
#' @importFrom Surrogate RandVec
#' @importFrom purrr pluck
make_vec_with_sum_of_int <- function (target_int, start_int, end_int, length_int) 
{
    vec_int <- Surrogate::RandVec(a = start_int, b = end_int, 
        s = target_int, n = length_int, m = 1) %>% purrr::pluck("RandVecOutput") %>% 
        as.vector() %>% round() %>% as.integer() %>% force_vec_to_sum_to_int(target_1L_int = target_int)
    return(vec_int)
}
