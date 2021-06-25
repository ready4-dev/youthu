#' Get model catalogue references
#' @description get_mdl_catalogue_refs() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get model catalogue references. Function argument predictors_chr specifies the where to look for the required object. The function returns Catalogue references (a character vector).
#' @param predictors_chr Predictors (a character vector)
#' @param ingredients_ls Ingredients (a list)
#' @return Catalogue references (a character vector)
#' @rdname get_mdl_catalogue_refs
#' @export 
#' @importFrom dplyr pull
#' @keywords internal
get_mdl_catalogue_refs <- function (predictors_chr, ingredients_ls) 
{
    catalogue_refs_chr <- get_mdls_using_predrs("k10", mdls_lup = ingredients_ls$mdls_lup) %>% 
        dplyr::pull(mdl_nms_chr)
    return(catalogue_refs_chr)
}
#' Get model from dataverse
#' @description get_mdl_from_dv() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get model from dataverse. Function argument mdl_nm_1L_chr specifies the where to look for the required object. The function returns Model (a model).
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @param dv_ds_nm_1L_chr Dataverse dataset name (a character vector of length one), Default: 'https://doi.org/10.7910/DVN/JC6PTV'
#' @param server_1L_chr Server (a character vector of length one), Default: 'dataverse.harvard.edu'
#' @param key_1L_chr Key (a character vector of length one), Default: NULL
#' @return Model (a model)
#' @rdname get_mdl_from_dv
#' @export 
#' @importFrom dataverse dataset_files
#' @importFrom purrr map_chr
get_mdl_from_dv <- function (mdl_nm_1L_chr, dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/JC6PTV", 
    server_1L_chr = "dataverse.harvard.edu", key_1L_chr = NULL) 
{
    ds_ls <- dataverse::dataset_files(dv_ds_nm_1L_chr, server = server_1L_chr, 
        key = key_1L_chr)
    all_mdls_chr <- purrr::map_chr(ds_ls, ~.x$label)
    idx_1L_int <- which(all_mdls_chr == paste0(mdl_nm_1L_chr, 
        ".RDS"))
    model_mdl <- readRDS(url(paste0("https://dataverse.harvard.edu/api/access/datafile/", 
        ds_ls[[idx_1L_int]]$dataFile$id)))
    return(model_mdl)
}
#' Get model summarys
#' @description get_mdl_smrys() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get model summarys. Function argument ingredients_ls specifies the where to look for the required object. The function returns Models summary (a list).
#' @param ingredients_ls Ingredients (a list)
#' @param mdl_nms_chr Model names (a character vector), Default: NULL
#' @return Models summary (a list)
#' @rdname get_mdl_smrys
#' @export 
#' @importFrom purrr map
#' @importFrom dplyr filter
#' @importFrom stats setNames
#' @keywords internal
get_mdl_smrys <- function (ingredients_ls, mdl_nms_chr = NULL) 
{
    if (is.null(mdl_nms_chr)) 
        mdl_nms_chr <- ingredients_ls$mdls_smry_tb$Model %>% 
            unique()
    mdls_smry_ls <- mdl_nms_chr %>% purrr::map(~ingredients_ls$mdls_smry_tb %>% 
        dplyr::filter(Model == .x)) %>% stats::setNames(mdl_nms_chr)
    return(mdls_smry_ls)
}
#' Get models using predictors
#' @description get_mdls_using_predrs() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get models using predictors. Function argument mdl_predrs_in_ds_chr specifies the where to look for the required object. The function returns Filtered models (a lookup table).
#' @param mdl_predrs_in_ds_chr Model predictors in dataset (a character vector)
#' @param mdls_lup Models (a lookup table), Default: NULL
#' @return Filtered models (a lookup table)
#' @rdname get_mdls_using_predrs
#' @export 
#' @importFrom utils data
#' @importFrom purrr map flatten map_lgl
#' @importFrom stats setNames
#' @importFrom rlang exec
#' @importFrom tidyr crossing unite
#' @importFrom dplyr filter mutate across everything pull
#' @importFrom stringr str_trim
#' @importFrom tibble as_tibble
get_mdls_using_predrs <- function (mdl_predrs_in_ds_chr, mdls_lup = NULL) 
{
    if (is.null(mdls_lup)) 
        utils::data("mdls_lup", envir = environment())
    args_ls <- mdl_predrs_in_ds_chr %>% purrr::map(~c(NA_character_, 
        .x)) %>% stats::setNames(paste0("var", 1:length(mdl_predrs_in_ds_chr)))
    tb <- rlang::exec(.fn = tidyr::crossing, !!!args_ls)
    include_lgl <- tb %>% dplyr::filter(tb %>% is.na() %>% rowSums() < 
        length(mdl_predrs_in_ds_chr)) %>% dplyr::mutate(dplyr::across(dplyr::everything(), 
        ~ifelse(is.na(.), "", .))) %>% tidyr::unite("combinations_chr", 
        colnames(tb), sep = " ", remove = T) %>% dplyr::pull(combinations_chr) %>% 
        stringr::str_trim() %>% purrr::map(~strsplit(.x, split = " ")) %>% 
        purrr::flatten() %>% purrr::map(~{
        terms_to_match_chr <- .x
        mdls_lup$predrs_ls %>% purrr::map_lgl(~{
            setdiff(.x, terms_to_match_chr) %>% identical(character(0))
        })
    }) %>% tibble::as_tibble(.name_repair = "unique") %>% rowSums() > 
        0
    filtered_mdls_lup <- mdls_lup %>% dplyr::filter(include_lgl)
    return(filtered_mdls_lup)
}
#' Get predictors
#' @description get_predictors() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get predictors. Function argument ingredients_ls specifies the where to look for the required object. The function returns Predictors (a tibble).
#' @param ingredients_ls Ingredients (a list)
#' @return Predictors (a tibble)
#' @rdname get_predictors
#' @export 
#' @importFrom dplyr select rename
#' @keywords internal
get_predictors <- function (ingredients_ls) 
{
    predictors_tb <- ingredients_ls$predictors_lup %>% dplyr::select(short_name_chr, 
        long_name_chr) %>% dplyr::rename(Variable = short_name_chr, 
        Description = long_name_chr)
    return(predictors_tb)
}
#' Get transformation from
#' @description get_tfmn_from_lup() is a Get function that retrieves a pre-existing data object from memory, local file system or online repository. Specifically, this function implements an algorithm to get transformation from lookup table. Function argument mdl_nm_1L_chr specifies the where to look for the required object. The function returns Transformation (a character vector of length one).
#' @param mdl_nm_1L_chr Model name (a character vector of length one)
#' @param mdls_lup Models (a lookup table), Default: NULL
#' @return Transformation (a character vector of length one)
#' @rdname get_tfmn_from_lup
#' @export 
#' @importFrom utils data
#' @importFrom ready4fun get_from_lup_obj
get_tfmn_from_lup <- function (mdl_nm_1L_chr, mdls_lup = NULL) 
{
    if (is.null(mdls_lup)) 
        utils::data("mdls_lup", envir = environment())
    tfmn_1L_chr <- ready4fun::get_from_lup_obj(mdls_lup, target_var_nm_1L_chr = "tfmn_chr", 
        match_value_xx = mdl_nm_1L_chr, match_var_nm_1L_chr = "mdl_nms_chr", 
        evaluate_lgl = F)
    return(tfmn_1L_chr)
}
