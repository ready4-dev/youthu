get_mdl_from_dv <- function(mdl_nm_1L_chr,
                            dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/JC6PTV",
                            server_1L_chr = "dataverse.harvard.edu",
                            key_1L_chr = NULL){
  ds_ls <- dataverse::dataset_files(dv_ds_nm_1L_chr,
                                    server = server_1L_chr,
                                    key = key_1L_chr)
  all_mdls_chr <- purrr::map_chr(ds_ls,~.x$label)
  idx_1L_int <- which(all_mdls_chr == paste0(mdl_nm_1L_chr,".RDS"))
  model_mdl <- readRDS(url(paste0("https://dataverse.harvard.edu/api/access/datafile/",ds_ls[[idx_1L_int]]$dataFile$id)))
  return(model_mdl)
}
get_mdls_using_predrs <- function(mdl_predrs_in_ds_chr,
                                  mdls_lup = NULL){
  if (is.null(mdls_lup))
    utils::data("mdls_lup", envir = environment())
  args_ls <- mdl_predrs_in_ds_chr %>% purrr::map(~c(NA_character_,.x)) %>% stats::setNames(paste0("var",1:length(mdl_predrs_in_ds_chr)))
  tb <- rlang::exec(.fn = tidyr::crossing, !!!args_ls)
  include_lgl <- tb %>%
    dplyr::filter(tb %>% is.na() %>% rowSums()<length(mdl_predrs_in_ds_chr)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),~ifelse(is.na(.),"",.))) %>%
    tidyr::unite("combinations_chr", colnames(tb), sep = " ",remove = T) %>%
    dplyr::pull(combinations_chr) %>%
    stringr::str_trim() %>%
    purrr::map(~strsplit(.x,split = " ")) %>%
    purrr::flatten() %>%
    purrr::map(~{
      terms_to_match_chr <- .x
      mdls_lup$predrs_ls %>% purrr::map_lgl(~ {
        setdiff(.x,terms_to_match_chr) %>% identical(character(0))}
      )
    }) %>%
    tibble::as_tibble(.name_repair = "unique")  %>%
    rowSums() > 0
  filtered_mdls_lup <- mdls_lup %>% dplyr::filter(include_lgl)
  return(filtered_mdls_lup)
}
get_signft_covars <- function (mdls_with_covars_smry_tb, covar_var_nms_chr)
{
    signif_vars_chr <- mdls_with_covars_smry_tb$Significant %>%
        purrr::map(~strsplit(.x, " ")) %>% purrr::flatten() %>%
        purrr::flatten_chr() %>% unique()
    signt_covars_chr <- covar_var_nms_chr[covar_var_nms_chr %in%
        signif_vars_chr]
    return(signt_covars_chr)
}
get_tfmn_from_lup <- function(mdl_nm_1L_chr, mdls_lup = NULL){
  if (is.null(mdls_lup))
    utils::data("mdls_lup", envir = environment())
  tfmn_1L_chr <- ready4fun::get_from_lup_obj(mdls_lup,
                                             target_var_nm_1L_chr = "tfmn_chr",
                                             match_value_xx = mdl_nm_1L_chr,
                                             match_var_nm_1L_chr = "mdl_nms_chr",
                                             evaluate_lgl = F)

  return(tfmn_1L_chr)
}
