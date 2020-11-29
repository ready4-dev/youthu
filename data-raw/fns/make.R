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
make_brms_mdl_print_ls <- function (mdl_ls, label_stub_1L_chr, caption_1L_chr, output_type_1L_chr = "PDF",
    digits_1L_dbl = 2, big_mark_1L_chr = " ")
{
    smry_mdl_ls <- summary(mdl_ls, digits = 4)
    mdl_smry_chr <- smry_mdl_ls %>% utils::capture.output()
    idx_dbl <- c("Formula: ", "Samples: ", "Group-Level Effects: ",
        "Population-Level Effects: ", "Family Specific Parameters: ",
        "Samples were drawn using ") %>% purrr::map_dbl(~mdl_smry_chr %>%
        startsWith(.x) %>% which())
    data_tb <- make_brms_mdl_smry_tbl(smry_mdl_ls, grp_1L_chr = mdl_smry_chr[idx_dbl[3]],
        pop_1L_chr = mdl_smry_chr[idx_dbl[4]], fam_1L_chr = mdl_smry_chr[idx_dbl[5]])
    bold_lgl <- data_tb$Parameter %in% c(mdl_smry_chr[idx_dbl[3]],
        mdl_smry_chr[idx_dbl[4]], mdl_smry_chr[idx_dbl[5]])
    if (output_type_1L_chr == "PDF") {
        data_tb <- data_tb %>% dplyr::mutate(Parameter = purrr::map_chr(Parameter,
            ~.x %>% Hmisc::latexTranslate()))
    }
    data_tb <- data_tb %>% dplyr::mutate(Parameter = Parameter %>%
        purrr::map2_chr(dplyr::all_of(bold_lgl), ~ifelse(.y &
            output_type_1L_chr == "PDF", paste0("\\textbf{",
            .x, "}"), .x)))
    if (output_type_1L_chr != "PDF") {
        data_tb <- data_tb %>% dplyr::mutate(dplyr::across(c(Bulk_ESS,
            Tail_ESS), ~format(., big.mark = big_mark_1L_chr)))
    }
    if (output_type_1L_chr == "HTML") {
        data_tb <- data_tb %>% dplyr::mutate(dplyr::across(where(is.numeric),
            ~format(round(., digits = digits_1L_dbl), digits = digits_1L_dbl,
                nsmall = digits_1L_dbl)))
    }
    data_tb <- data_tb %>% dplyr::mutate(dplyr::across(where(is.character),
        ~dplyr::case_when(is.na(.) ~ "", . == "NA" ~ "", endsWith(.,
            " NA") ~ "", TRUE ~ .)))
    end_matter_1L_chr <- trimws(mdl_smry_chr[idx_dbl[6]:length(mdl_smry_chr)]) %>%
        paste0(collapse = " ")
    brms_mdl_print_ls <- list(part_1 = mdl_smry_chr[idx_dbl[1]],
        part_2 = "\n\n", part_3 = c(trimws(mdl_smry_chr[1:(idx_dbl[2] -
            1)][-idx_dbl[1]]), paste0(trimws(mdl_smry_chr[idx_dbl[2]]),
            " ", trimws(mdl_smry_chr[idx_dbl[2] + 1]), collapse = " ")) %>%
            paste0(collapse = ifelse(output_type_1L_chr == "PDF",
                "\n\n", "\n")), part_4 = "\n\n", part_5 = list(data_tb = data_tb,
            output_type_1L_chr = output_type_1L_chr, caption_1L_chr = caption_1L_chr,
            label = paste0("tab:", label_stub_1L_chr), merge_row_idx_int = as.integer(which(bold_lgl)),
            digits_dbl = c(ifelse(output_type_1L_chr == "PDF",
                0, NA_real_) %>% purrr::discard(is.na), names(data_tb) %>%
                purrr::map_dbl(~ifelse(.x %in% c("Bulk_ESS",
                  "Tail_ESS"), 0, digits_1L_dbl))), big_mark_1L_chr = big_mark_1L_chr,
            hline.after = c(-1, 0), sanitize_fn = force, footnotes_chr = NA_character_),
        part_6 = end_matter_1L_chr)
    if (output_type_1L_chr != "PDF") {
        brms_mdl_print_ls$part_5$footnotes_chr <- c(paste0(brms_mdl_print_ls$part_1,
            ifelse(output_type_1L_chr == "Word", "", "\n")),
            brms_mdl_print_ls$part_3, brms_mdl_print_ls$part_6)
        brms_mdl_print_ls$part_6 <- NULL
    }
    else {
        footnotes_chr <- c(mdl_smry_chr[idx_dbl[1]], trimws(mdl_smry_chr[1:(idx_dbl[2] -
            1)][-idx_dbl[1]]), trimws(mdl_smry_chr[idx_dbl[2]]),
            trimws(mdl_smry_chr[idx_dbl[2] + 1]), trimws(mdl_smry_chr[idx_dbl[6]:length(mdl_smry_chr)])) %>%
            Hmisc::latexTranslate()
        footnotes_chr[1] <- footnotes_chr[1] %>% stringr::str_replace("~",
            "\\\\textasciitilde")
        brms_mdl_print_ls$part_5$addtorow <- list(pos = purrr::map(c(0,
            rep(nrow(data_tb), length(footnotes_chr) + 1)), ~.x),
            command = c(names(data_tb) %>% Hmisc::latexTranslate() %>%
                paste0(collapse = " & ") %>% paste0("\\\\\n"),
                c("\\toprule\n"), footnotes_chr %>% purrr::map_chr(~paste0("\\multicolumn{",
                  ncol(data_tb), "}{l}{", paste0("{\\footnotesize ",
                    .x, "}\n", collapse = ","), "}\\\\\n"))))
    }
    return(brms_mdl_print_ls)
}
make_brms_mdl_smry_tbl <- function (smry_mdl_ls, grp_1L_chr, pop_1L_chr, fam_1L_chr)
{
    brms_mdl_smry_tb <- purrr::map(1:length(smry_mdl_ls$random),
        ~make_mdl_smry_elmt_tbl(cat_chr = c(ifelse(.x == 1, grp_1L_chr,
            character(0)), paste0(names(smry_mdl_ls$ngrps)[.x],
            " (Number of levels: ", smry_mdl_ls$ngrps[.x][[1]],
            ")")), mat = smry_mdl_ls$random[.x][[1]])) %>% dplyr::bind_rows(make_mdl_smry_elmt_tbl(mat = smry_mdl_ls$fixed,
        cat_chr = pop_1L_chr), make_mdl_smry_elmt_tbl(mat = smry_mdl_ls$spec_pars,
        cat_chr = fam_1L_chr))
    return(brms_mdl_smry_tb)
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
make_corstars_tbl_xx <- function (x, method = c("pearson", "spearman"), removeTriangle = c("upper",
    "lower"), result = c("none", "html", "latex"))
{
    x <- as.matrix(x)
    correlation_matrix <- Hmisc::rcorr(x, type = method[1])
    R <- correlation_matrix$r
    p <- correlation_matrix$P
    mystars <- ifelse(p < 1e-04, "****", ifelse(p < 0.001, "*** ",
        ifelse(p < 0.01, "**  ", ifelse(p < 0.05, "*   ", "    "))))
    R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[, -1]
    Rnew <- matrix(paste(R, mystars, sep = ""), ncol = ncol(x))
    diag(Rnew) <- paste(diag(R), " ", sep = "")
    rownames(Rnew) <- colnames(x)
    colnames(Rnew) <- paste(colnames(x), "", sep = "")
    if (removeTriangle[1] == "upper") {
        Rnew <- as.matrix(Rnew)
        Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
        Rnew <- as.data.frame(Rnew)
    }
    else if (removeTriangle[1] == "lower") {
        Rnew <- as.matrix(Rnew)
        Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
        Rnew <- as.data.frame(Rnew)
    }
    Rnew <- cbind(Rnew[1:length(Rnew) - 1])
    if (result[1] == "none")
        return(Rnew)
    else {
        if (result[1] == "html")
            print(xtable(Rnew), type = "html")
        else print(xtable(Rnew), type = "latex")
    }
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
make_folds_ls <- function (data_tb, dep_var_nm_1L_chr = "aqol6d_total_w", n_folds_1L_int = 10L)
{
    folds_ls <- caret::createFolds(data_tb %>% dplyr::pull(!!rlang::sym(dep_var_nm_1L_chr)),
        k = n_folds_1L_int, list = TRUE, returnTrain = FALSE)
    return(folds_ls)
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
make_knit_pars_ls <- function (mdl_smry_dir_1L_chr, mdl_types_chr, predr_vars_nms_ls,
    output_type_1L_chr = "HTML", mdl_types_lup = NULL, plt_types_lup = NULL,
    plt_types_chr = c("coefs", "hetg", "dnst", "sctr_plt"), section_type_1L_chr = "#")
{
    if (is.null(mdl_types_lup))
        utils::data(mdl_types_lup, envir = environment())
    if (is.null(plt_types_lup))
        utils::data(plt_types_lup, envir = environment())
    lab_idx_dbl <- 1:(length(mdl_types_chr) * length(predr_vars_nms_ls))
    knit_pars_ls <- purrr::pmap(list(predr_vars_nms_ls, split(lab_idx_dbl,
        ceiling(seq_along(lab_idx_dbl)/length(mdl_types_chr))),
        make_unique_ls_elmt_idx_int(predr_vars_nms_ls), make_mdl_nms_ls(predr_vars_nms_ls,
            mdl_types_chr = mdl_types_chr)), ~{
        mdl_nms_chr <- ..4
        path_to_mdl_stubs_chr <- paste0(mdl_smry_dir_1L_chr,
            "/", mdl_nms_chr)
        paths_to_mdls_chr <- paste0(path_to_mdl_stubs_chr, ".RDS")
        paths_to_mdl_plts_ls <- purrr::map(path_to_mdl_stubs_chr,
            ~paste0(..1, paste0("_", plt_types_chr, ".png")))
        mdl_ttls_chr <- paste0(..1[1], ifelse(is.na(..1[2]),
            "", paste(" with ", ..1[2])), " ", ready4fun::get_from_lup_obj(mdl_types_lup,
            match_var_nm_1L_chr = "short_name_chr", match_value_xx = mdl_types_chr,
            target_var_nm_1L_chr = "long_name_chr", evaluate_lgl = F))
        section_ttls_chr <- paste0(section_type_1L_chr, " ",
            mdl_ttls_chr)
        list(plt_nms_ls = purrr::map(mdl_ttls_chr, ~{
            paste0(.x, " ", ready4fun::get_from_lup_obj(plt_types_lup,
                match_var_nm_1L_chr = "short_name_chr", match_value_xx = plt_types_chr,
                target_var_nm_1L_chr = "long_name_chr", evaluate_lgl = F))
        }), paths_to_mdls_chr = paths_to_mdls_chr, tbl_captions_chr = mdl_ttls_chr,
            label_stubs_chr = paste0("lab", ..2), output_type_1L_chr = output_type_1L_chr,
            section_ttls_chr = section_ttls_chr, ls_elmt_idx_1L_int = ..3)
    })
    return(knit_pars_ls)
}
make_mdl_smry_elmt_tbl <- function (mat, cat_chr)
{
    tb <- mat %>% tibble::as_tibble() %>% dplyr::mutate(Parameter = rownames(mat)) %>%
        dplyr::select(Parameter, dplyr::everything())
    mdl_elmt_sum_tb <- tb %>% dplyr::filter(F) %>% tibble::add_case(Parameter = cat_chr) %>%
        dplyr::bind_rows(tb)
    return(mdl_elmt_sum_tb)
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
make_predr_vals <- function (predr_var_nm_1L_chr, candidate_predrs_lup = NULL)
{
    if (is.null(candidate_predrs_lup)) {
        utils::data("candidate_predrs_lup", envir = environment())
    }
    args_ls <- purrr::map_dbl(names(candidate_predrs_lup)[3:5],
        ~candidate_predrs_lup %>% ready4fun::get_from_lup_obj(match_value_xx = predr_var_nm_1L_chr,
            match_var_nm_1L_chr = "short_name_chr", target_var_nm_1L_chr = .x,
            evaluate_lgl = F)) %>% as.list()
    predr_vals_dbl <- rlang::exec(.fn = seq, !!!args_ls)
    return(predr_vals_dbl)
}
make_predr_vars_nms_ls <- function (main_predrs_chr, covars_ls)
{
    predr_vars_nms_ls <- covars_ls %>% purrr::map(~{
        covars_chr <- .x
        purrr::map(main_predrs_chr, ~list(c(.x), c(.x, covars_chr))) %>%
            purrr::flatten()
    }) %>% purrr::flatten()
    predr_vars_nms_ls <- predr_vars_nms_ls[order(sapply(predr_vars_nms_ls,
        length))]
    return(predr_vars_nms_ls)
}
make_prefd_mdls_vec <- function (smry_of_sngl_predr_mdls_tb, choose_from_pfx_chr = c("GLM",
    "OLS"), mdl_types_lup = NULL)
{
    if (is.null(mdl_types_lup))
        utils::data("mdl_types_lup", envir = environment())
    ordered_mdl_types_chr <- dplyr::inner_join(smry_of_sngl_predr_mdls_tb %>%
        dplyr::select(Model) %>% dplyr::rename(long_name_chr = Model),
        mdl_types_lup) %>% dplyr::pull(short_name_chr)
    prefd_mdls_chr <- purrr::map_chr(choose_from_pfx_chr, ~ordered_mdl_types_chr[startsWith(ordered_mdl_types_chr,
        .x)][1])
    return(prefd_mdls_chr)
}
make_synth_series_tbs_ls <- function (synth_data_spine_ls, series_names_chr)
{
    synth_series_tbs_ls <- 1:length(series_names_chr) %>% purrr::map(~make_correlated_data_tb(synth_data_spine_ls = synth_data_spine_ls,
        synth_data_idx_1L_dbl = .x) %>% replace_with_missing_vals(synth_data_spine_ls = synth_data_spine_ls,
        idx_int = .x)) %>% stats::setNames(series_names_chr)
    return(synth_series_tbs_ls)
}
make_unique_ls_elmt_idx_int <- function (data_ls)
{
    combos_tb <- tibble::as_tibble(data_ls, .name_repair = ~paste0("r_",
        1:length(data_ls))) %>% t() %>% as.data.frame()
    combos_tb <- combos_tb %>% tibble::as_tibble()
    combos_tb <- combos_tb %>% dplyr::mutate(V2 = dplyr::case_when(V1 ==
        V2 ~ NA_character_, T ~ V2)) %>% dplyr::mutate(combo_chr = purrr::map2_chr(V1,
        V2, ~ifelse(ncol(combos_tb) == 1 | is.na(.y), .x, paste0(.x,
            "_", .y))))
    combos_tb <- combos_tb %>% dplyr::group_by(combo_chr) %>%
        dplyr::mutate(combo_id = dplyr::row_number())
    unique_ls_elmt_idx_int <- purrr::map(data_ls %>% unique(),
        ~ready4fun::get_from_lup_obj(combos_tb %>% dplyr::ungroup(),
            match_var_nm_1L_chr = "combo_chr", match_value_xx = paste0(.x[1],
                ifelse(is.na(.x[2]), "", paste0("_", .x[2]))),
            target_var_nm_1L_chr = "combo_id", evaluate_lgl = F)) %>%
        purrr::flatten_int()
    return(unique_ls_elmt_idx_int)
}
make_vec_with_sum_of_int <- function (target_int, start_int, end_int, length_int)
{
    vec_int <- Surrogate::RandVec(a = start_int, b = end_int,
        s = target_int, n = length_int, m = 1) %>% purrr::pluck("RandVecOutput") %>%
        as.vector() %>% round() %>% as.integer() %>% force_vec_to_sum_to_int(target_1L_int = target_int)
    return(vec_int)
}
