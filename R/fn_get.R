#' Get mdl from dataverse
#' @description get_mdl_from_dv() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get mdl from dataverse. Function argument mdl_nm_1L_chr specifies the where to look for the required object. The function is called for its side effects and does not return a value.
#' @param mdl_nm_1L_chr Mdl name (a character vector of length one)
#' @param dv_ds_nm_1L_chr Dataverse dataset name (a character vector of length one), Default: 'https://doi.org/10.7910/DVN/JC6PTV'
#' @return NA ()
#' @rdname get_mdl_from_dv
#' @export 
#' @importFrom dataverse dataset_files
#' @importFrom purrr map_chr
#' @keywords internal
get_mdl_from_dv <- function (mdl_nm_1L_chr, dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/JC6PTV") 
{
    ds_ls <- dataverse::dataset_files(dv_ds_nm_1L_chr)
    all_mdls_chr <- purrr::map_chr(ds_ls, ~.x$label)
    idx_1L_int <- which(all_mdls_chr == paste0(mdl_nm_1L_chr, 
        ".RDS"))
    model_mdl <- readRDS(url(paste0("https://dataverse.harvard.edu/api/access/datafile/", 
        ds_ls[[idx_1L_int]]$dataFile$id)))
    return(model_mdl)
}
#' Get signft covars
#' @description get_signft_covars() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get signft covars. Function argument mdls_with_covars_smry_tb specifies the where to look for the required object. The function returns Signt covars (a character vector).
#' @param mdls_with_covars_smry_tb Mdls with covars smry (a tibble)
#' @param covar_var_nms_chr Covar var names (a character vector)
#' @return Signt covars (a character vector)
#' @rdname get_signft_covars
#' @export 
#' @importFrom purrr map flatten flatten_chr
#' @keywords internal
get_signft_covars <- function (mdls_with_covars_smry_tb, covar_var_nms_chr) 
{
    signif_vars_chr <- mdls_with_covars_smry_tb$Significant %>% 
        purrr::map(~strsplit(.x, " ")) %>% purrr::flatten() %>% 
        purrr::flatten_chr() %>% unique()
    signt_covars_chr <- covar_var_nms_chr[covar_var_nms_chr %in% 
        signif_vars_chr]
    return(signt_covars_chr)
}
#' Get tfmn from
#' @description get_tfmn_from_lup() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get tfmn from lookup table. Function argument mdl_nm_1L_chr specifies the where to look for the required object. The function returns Tfmn (a character vector of length one).
#' @param mdl_nm_1L_chr Mdl name (a character vector of length one)
#' @param mdls_lup Mdls (a lookup table), Default: NULL
#' @return Tfmn (a character vector of length one)
#' @rdname get_tfmn_from_lup
#' @export 
#' @importFrom utils data
#' @importFrom ready4fun get_from_lup_obj
#' @keywords internal
get_tfmn_from_lup <- function (mdl_nm_1L_chr, mdls_lup = NULL) 
{
    if (is.null(mdls_lup)) 
        utils::data("mdls_lup", envir = environment())
    tfmn_1L_chr <- ready4fun::get_from_lup_obj(mdls_lup, target_var_nm_1L_chr = "tfmn_chr", 
        match_value_xx = mdl_nm_1L_chr, match_var_nm_1L_chr = "mdl_nms_chr", 
        evaluate_lgl = F)
    return(tfmn_1L_chr)
}
