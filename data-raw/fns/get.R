get_mdl_from_dv <- function(mdl_nm_1L_chr,
                            dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/JC6PTV"){
  ds_ls <- dataverse::dataset_files(dv_ds_nm_1L_chr)
  all_mdls_chr <- purrr::map_chr(ds_ls,~.x$label)
  idx_1L_int <- which(all_mdls_chr == paste0(mdl_nm_1L_chr,".RDS"))
  model_mdl <- readRDS(url(paste0("https://dataverse.harvard.edu/api/access/datafile/",ds_ls[[idx_1L_int]]$dataFile$id)))
  return(model_mdl)
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
