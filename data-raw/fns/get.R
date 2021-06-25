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
get_mdl_catalogue_refs <- function(predictors_chr,
                                   ingredients_ls){
  catalogue_refs_chr <- get_mdls_using_predrs("k10",
                                              mdls_lup = ingredients_ls$mdls_lup) %>%
    dplyr::pull(mdl_nms_chr)
  return(catalogue_refs_chr)
}
get_mdl_smrys <- function(ingredients_ls,
                          mdl_nms_chr = NULL){
  if(is.null(mdl_nms_chr))
    mdl_nms_chr <- ingredients_ls$mdls_smry_tb$Model %>% unique()
  mdls_smry_ls <- mdl_nms_chr %>%
    purrr::map(~ingredients_ls$mdls_smry_tb %>%
                 dplyr::filter(Model == .x)) %>%
    stats::setNames(mdl_nms_chr)
  return(mdls_smry_ls)
}
get_predictors <- function(ingredients_ls){
  predictors_tb <- ingredients_ls$predictors_lup %>%
    dplyr::select(short_name_chr, long_name_chr) %>%
    dplyr::rename(Variable = short_name_chr,
                  Description = long_name_chr)
  return(predictors_tb)
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
