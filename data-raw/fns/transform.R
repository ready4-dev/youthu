transform_ds_for_cmprsn <- function(ds_tb,
                                    cmprsn_var_nm_1L_chr,
                                    id_var_nm_1L_chr = "UID_chr",
                                    round_var_nm_1L_chr = "Timepoint_chr",
                                    cmprsn_groups_chr = c("Intervention","Control")){
  ds_tb <- ds_tb %>%
    dplyr::group_by(!!rlang::sym(id_var_nm_1L_chr)) %>%
    dplyr::mutate(n_measures_int = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n_measures_int == 2)
  UIDs_with_FUP_data_chr <- ds_tb  %>% dplyr::pull(!!rlang::sym(id_var_nm_1L_chr)) %>% unique()
  group_1_IDs <- sample(UIDs_with_FUP_data_chr, floor(length(UIDs_with_FUP_data_chr)/2))
  group_2_IDs <- setdiff(UIDs_with_FUP_data_chr,group_1_IDs)[1:length(group_1_IDs)]
  ds_tb <- ds_tb %>% dplyr::mutate(!!rlang::sym(cmprsn_var_nm_1L_chr) := dplyr::case_when(!!rlang::sym(id_var_nm_1L_chr) %in% group_1_IDs ~ cmprsn_groups_chr[1],
                                                                                          !!rlang::sym(id_var_nm_1L_chr) %in% group_2_IDs ~ cmprsn_groups_chr[2],
                                                                                          T ~ NA_character_
  )) %>%
    dplyr::select(-n_measures_int) %>%
    dplyr::filter(!!rlang::sym(cmprsn_var_nm_1L_chr) %in% cmprsn_groups_chr)
  return(ds_tb)
}
transform_ds_to_drop_msng <- function(ds_tb,
                                       predictors_chr,
                                       uid_var_nm_1L_chr = "UID_chr"){
  drop_chr <- ds_tb[rowSums(is.na(ds_tb[predictors_chr])) >0,] %>% dplyr::pull(!!rlang::sym(uid_var_nm_1L_chr))
  ds_tb <- dplyr::filter(ds_tb, !(!!rlang::sym(uid_var_nm_1L_chr) %in% drop_chr))
  return(ds_tb)
}
transform_ds_to_long <- function(ds_tb,
                              predictors_chr,
                              drop_underscore_1L_lgl = T,
                              msrmnt_date_var_nm_1L_chr = "date_dtm",
                              round_var_nm_1L_chr = "Timepoint_chr",
                              row_id_nm_1L_chr = "case_id",
                              time_is_sfx_1L_lgl = T#,
                              #uid_var_nm_1L_chr = "UID_chr",
                              ){
  names_chr <- names(ds_tb)
  if(time_is_sfx_1L_lgl){
    crop_fn <- stringi::stri_replace_first_fixed
    select_fn <- startsWith
  }else{
    crop_fn <- stringi::stri_replace_last_fixed
    select_fn <- endsWith
  }
  predictors_ls <- predictors_chr %>% purrr::map(~{
    predictor_1L_chr <- .x
    predictor_vars_chr <- names_chr[names_chr %>% purrr::map_lgl(~ {
      name_1L_chr <- .x
      #predictors_chr %>% # purrr::map_lgl(~
      select_fn(name_1L_chr, predictor_1L_chr)
                                       # )
      })]
  }) %>% stats::setNames(predictors_chr)
  prefixes_chr <- suffixes_chr <- character(0)
  extensions_chr <- predictors_ls[[1]] %>% purrr::map_chr(~crop_fn(.x, pattern = names(predictors_ls)[1], replacement = ""))
  if(time_is_sfx_1L_lgl){
    suffixes_chr <- extensions_chr
  }else{
    prefixes_chr <- extensions_chr
  }
  if(drop_underscore_1L_lgl){
    tfmn_fn <- function(x){stringr::str_remove_all(x,"_")}
  }else{
    tfmn_fn <- identity
  }
  predictor_vars_chr <- predictors_ls %>% purrr::flatten_chr()
  other_vars_chr <- setdiff(names_chr, c(paste0(prefixes_chr,msrmnt_date_var_nm_1L_chr,suffixes_chr), predictor_vars_chr))
  ds_tb <- c(msrmnt_date_var_nm_1L_chr, predictors_chr) %>%
    purrr::map(~ ds_tb %>%
                 tidyr::pivot_longer(cols = paste0(prefixes_chr,.x,suffixes_chr),
                                     names_to = round_var_nm_1L_chr, values_to = .x) %>%
                 tibble::rowid_to_column(row_id_nm_1L_chr) %>%
                 dplyr::mutate(!!rlang::sym(round_var_nm_1L_chr) := !!rlang::sym(round_var_nm_1L_chr) %>%
                                 crop_fn(pattern = .x, replacement = "") %>%
                                 tfmn_fn() %>%
                                 #stringr::str_sub(start = -2L) %>%
                                 factor()) %>%
                 dplyr::select(tidyselect::all_of(c(row_id_nm_1L_chr, other_vars_chr, round_var_nm_1L_chr, .x)))) %>%
    purrr::reduce(~dplyr::left_join(.x,.y)) %>%
    dplyr::select(-tidyselect::all_of(row_id_nm_1L_chr))
  return(ds_tb)
  }
  # data_tb <- c("d_interview_date", "gad7", "phq9") %>%
  #   purrr::map(~ data_tb %>%
  #                tidyr::pivot_longer(cols = paste0(.x,"_t",1:2),
  #                                    names_to = "When", values_to = .x) %>%
  #                tibble::rowid_to_column("case_id") %>%
  #                dplyr::mutate(When = When %>% stringr::str_sub(start = -2L) %>% factor()) %>%
  #                dplyr::select(tidyselect::all_of(c("case_id", "fkClientID","When",.x)))) %>%
  #   purrr::reduce(~dplyr::left_join(.x,.y)) %>%
  #   dplyr::select(-case_id)



