get_dv_ds_publication <- function(ds_url_1L_chr,
                                  server_1L_chr = "dataverse.harvard.edu",
                                  key_1L_chr = NULL){
  ds_md_ls <- dataverse::dataset_metadata(ds_url_1L_chr,
                                          key = key_1L_chr,
                                          server = server_1L_chr)
  doi_url_1L_chr <- ds_md_ls$fields %>%
    dplyr::filter(typeName == "publication") %>%
    dplyr::pull(value) %>% purrr::pluck(1)
  doi_url_1L_chr <- ifelse(!is.null(doi_url_1L_chr),doi_url_1L_chr %>%
                             dplyr::pull(publicationURL) %>%
                             dplyr::pull(value), "")
  return(doi_url_1L_chr)
}
get_dv_dss_mdl_smrys <- function(ids_chr,
                                 server_1L_chr = "dataverse.harvard.edu",
                                 key_1L_chr = NULL){
  dv_dss_mdl_smrys_ls <- ids_chr %>% purrr::map(~get_mdl_from_dv("mdl_ingredients",
                                                                 dv_ds_nm_1L_chr = .x,
                                                                 server_1L_chr = server_1L_chr,
                                                                 key_1L_chr = key_1L_chr)) %>%
    stats::setNames(ids_chr)
  return(dv_dss_mdl_smrys_ls)
}
get_dv_mdl_smrys <- function(mdls_lup,
                             mdl_nms_chr = NULL){
  ds_urls_chr <- mdls_lup$ds_url %>% unique()
  dss_mdl_smrys_ls <- get_dv_dss_mdl_smrys(ds_urls_chr)
  dv_mdl_smrys <- dss_mdl_smrys_ls %>% purrr::map2(ds_urls_chr,
                                                   ~get_mdl_smrys(.x,
                                                                  mdl_nms_chr = mdls_lup %>%
                                                                    dplyr::filter(ds_url ==.y) %>%
                                                                    dplyr::pull(mdl_nms_chr))) %>%
    purrr::flatten()
  if(!is.null(mdl_nms_chr)){
    dv_mdl_smrys <- dv_mdl_smrys[names(dv_mdl_smrys) %in% mdl_nms_chr]

  }
  return(dv_mdl_smrys)
}
get_filtered_ttu_dss <- function(ttu_dv_dss_tb = NULL,
                                 mdl_predrs_in_ds_chr = NULL,
                                 utility_type_chr = NULL,
                                 ttu_dv_nms_chr = "TTU",
                                 server_1L_chr = "dataverse.harvard.edu",
                                 key_1L_chr = NULL){
  if(is.null(ttu_dv_dss_tb))
    ttu_dv_dss_tb <- get_ttu_dv_dss(ttu_dv_nms_chr = ttu_dv_nms_chr,
                                    server_1L_chr = server_1L_chr,
                                    key_1L_chr = NULL)
  if(is.null(ttu_dv_dss_tb)){
    if(is.null(mdl_predrs_in_ds_chr))
      mdl_predrs_in_ds_chr <- get_ttu_dv_predrs(ttu_dv_dss_tb)
    ttu_dv_dss_tb <- ttu_dv_dss_tb %>%
      dplyr::filter(predrs_ls %>% purrr::map_lgl(~!identical(intersect(.x,mdl_predrs_in_ds_chr),
                                                             character(0))))
    if(!is.null(utility_type_chr))
      ttu_dv_dss_tb <- ttu_dv_dss_tb %>%
      dplyr::filter(utility %in% utility_type_chr)
  }
  return(ttu_dv_dss_tb)
}
get_mdl_ctlg_url <- function(mdls_lup,
                             mdl_nm_1L_chr,
                             server_1L_chr = "dataverse.harvard.edu",
                             key_1L_chr = NULL){
  dv_ds_nm_1L_chr <- get_mdl_ds_url(mdls_lup, mdl_nm_1L_chr = mdl_nm_1L_chr)
  ds_ls <- dataverse::dataset_files(dv_ds_nm_1L_chr,
                                    server = server_1L_chr,
                                    key = key_1L_chr)
  all_lbls_chr <- purrr::map_chr(ds_ls,~.x$label)
  include_lgl <- all_lbls_chr %>% purrr::map_lgl(~startsWith(.x,"AAA_TTU_MDL_CTG"))
  all_descs_chr <- purrr::map_chr(ds_ls,~.x$description)
  include_lgl <- include_lgl & (all_descs_chr %>% purrr::map_lgl(~stringr::str_detect(.x,
                                                                                     ready4::get_from_lup_obj(mdls_lup,
                                                                                                                 match_value_xx = mdl_nm_1L_chr,
                                                                                                                 match_var_nm_1L_chr = "mdl_nms_chr",
                                                                                                                 target_var_nm_1L_chr = "source_chr",
                                                                                                                 evaluate_1L_lgl = F))))
  idx_1L_int <- which(include_lgl)
  if(identical(idx_1L_int,integer(0))){
    ctlg_url <- NULL
  }else{
    ctlg_url <- paste0("https://dataverse.harvard.edu/api/access/datafile/",ds_ls[[idx_1L_int]]$dataFile$id)
  }
  return(ctlg_url)
}
get_mdl_ds_url <- function(mdls_lup,
                           mdl_nm_1L_chr){
  mdl_ds_url <- ready4::get_from_lup_obj(mdls_lup,
                                                 match_value_xx = mdl_nm_1L_chr,
                                                 match_var_nm_1L_chr = "mdl_nms_chr",
                                                 target_var_nm_1L_chr = "ds_url",
                                                 evaluate_1L_lgl = F)
  return(mdl_ds_url)
}
get_mdl_from_dv <- function(mdl_nm_1L_chr,
                            dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/JC6PTV",
                            server_1L_chr = "dataverse.harvard.edu",
                            key_1L_chr = NULL){
  ds_ls <- dataverse::dataset_files(dv_ds_nm_1L_chr,
                                    server = server_1L_chr,
                                    key = key_1L_chr)
  all_mdls_chr <- purrr::map_chr(ds_ls,~.x$label)
  idx_1L_int <- which(all_mdls_chr == paste0(mdl_nm_1L_chr,".RDS"))
  if(identical(idx_1L_int,integer(0))){
    model_mdl <- NULL
  }else{
    model_mdl <- readRDS(url(paste0("https://dataverse.harvard.edu/api/access/datafile/",ds_ls[[idx_1L_int]]$dataFile$id)))
  }
  return(model_mdl)
}
get_mdls_lup <- function(ttu_dv_dss_tb = NULL,
                         mdl_predrs_in_ds_chr = NULL,
                         utility_type_chr = NULL,
                         ttu_dv_nms_chr = "TTU",
                         server_1L_chr = "dataverse.harvard.edu",
                         key_1L_chr = NULL){
  if(is.null(ttu_dv_dss_tb))
    ttu_dv_dss_tb <- get_ttu_dv_dss(ttu_dv_nms_chr = ttu_dv_nms_chr,
                                    server_1L_chr = server_1L_chr,
                                    key_1L_chr = NULL)
  if(!is.null(ttu_dv_dss_tb)){
    if(is.null(mdl_predrs_in_ds_chr))
      mdl_predrs_in_ds_chr <- get_ttu_dv_predrs(ttu_dv_dss_tb)
    ttu_dv_dss_tb <- get_filtered_ttu_dss(ttu_dv_dss_tb,
                                          mdl_predrs_in_ds_chr = mdl_predrs_in_ds_chr,
                                          utility_type_chr = utility_type_chr)
  }
  if(!is.null(ttu_dv_dss_tb)){
    ds_smrys_ls <- get_ttu_ds_smrys("TTU", reference_int = ttu_dv_dss_tb$reference_int)
    mdls_lup <- ds_smrys_ls %>%
      purrr::map2_dfr(names(ds_smrys_ls),
                      ~{
                        predictors_lup <- .x$predictors_lup
                        urls_chr <- .y
                        predrs_short_nms_chr <- mdl_predrs_in_ds_chr %>%
                          purrr::map_chr(~ready4::get_from_lup_obj(predictors_lup,
                                                                      match_value_xx = .x,
                                                                      match_var_nm_1L_chr = "long_name_chr",
                                                                      target_var_nm_1L_chr = "short_name_chr",
                                                                      evaluate_1L_lgl = F))
                        .x$mdls_lup %>%
                          dplyr::filter(predrs_ls %>%
                                          purrr::map_lgl(~{!identical(intersect(.x,
                                                                                predrs_short_nms_chr),
                                                                      character(0))})) %>%
                          dplyr::mutate(ds_url = urls_chr)
                      })
  }else{
    mdls_lup <- NULL
  }
  return(mdls_lup)
}
# get_mdl_catalogue_refs <- function(predictors_chr,
#                                    ingredients_ls,
#                                    mdl_predrs_in_ds_chr){
#   catalogue_refs_chr <- get_mdls_lup(mdl_predrs_in_ds_chr = mdl_predrs_in_ds_chr,
#                                               mdls_lup = ingredients_ls$mdls_lup) %>%
#     dplyr::pull(mdl_nms_chr)
#   return(catalogue_refs_chr)
# }
get_mdl_smrys <- function(ingredients_ls,
                          mdl_nms_chr = NULL){
  if(is.null(mdl_nms_chr))
    mdl_nms_chr <- ingredients_ls$mdls_smry_tb$Model %>% unique()
  mdls_smry_ls <- mdl_nms_chr %>%
    purrr::map(~{
      df <- ingredients_ls$mdls_smry_tb %>%
                 dplyr::filter(Model == .x) %>%
                 dplyr::select(-Model)
      rownames(df) <- NULL
      df
      }) %>%
    stats::setNames(mdl_nms_chr)
  return(mdls_smry_ls)
}
get_mdl_metadata <- function(mdls_lup,
                             mdl_nm_1L_chr,
                             server_1L_chr = "dataverse.harvard.edu",
                             key_1L_chr = NULL){
  dv_ds_nm_1L_chr <- get_mdl_ds_url(mdls_lup, mdl_nm_1L_chr = mdl_nm_1L_chr)
  ingredients_ls <- get_mdl_from_dv("mdl_ingredients",
                                    dv_ds_nm_1L_chr = dv_ds_nm_1L_chr,
                                    server_1L_chr = server_1L_chr,
                                    key_1L_chr = key_1L_chr)
  return(ingredients_ls)
}
get_model <- function(mdls_lup,
                      mdl_nm_1L_chr,
                      make_from_tbl_1L_lgl = T,
                      mdl_meta_data_ls = NULL,
                      server_1L_chr = "dataverse.harvard.edu",
                      key_1L_chr = NULL){
  if(make_from_tbl_1L_lgl){
    if(is.null(mdl_meta_data_ls)){
      mdl_meta_data_ls <- get_mdl_metadata(mdls_lup,
                                         mdl_nm_1L_chr = mdl_nm_1L_chr,
                                         server_1L_chr = server_1L_chr,
                                         key_1L_chr = key_1L_chr)
    }
    model_mdl <- TTU::get_table_predn_mdl(mdl_nm_1L_chr,
                                          ingredients_ls = mdl_meta_data_ls,
                                          analysis_1L_chr = ready4::get_from_lup_obj(mdls_lup,
                                                                                   match_value_xx = mdl_nm_1L_chr,
                                                                                   match_var_nm_1L_chr = "mdl_nms_chr",
                                                                                   target_var_nm_1L_chr = "source_chr",
                                                                                   evaluate_1L_lgl = F))
  }else{
    model_mdl <- get_mdl_from_dv(mdl_nm_1L_chr,
                                 dv_ds_nm_1L_chr = get_mdl_ds_url(mdls_lup,
                                                                  mdl_nm_1L_chr = mdl_nm_1L_chr),
                                 server_1L_chr = server_1L_chr,
                                 key_1L_chr = key_1L_chr)
  }
  return(model_mdl)
}
get_predictors_lup <- function(mdl_meta_data_ls = NULL,
                               mdls_lup = NULL,
                               mdl_nm_1L_chr = NULL,
                               outp_is_abbrs_tb = F,
                               server_1L_chr = "dataverse.harvard.edu",
                               key_1L_chr = NULL){
  if(is.null(mdl_meta_data_ls))
    mdl_meta_data_ls <- get_mdl_metadata(mdls_lup = mdls_lup,
                                         mdl_nm_1L_chr = mdl_nm_1L_chr,
                                         server_1L_chr = server_1L_chr,
                                         key_1L_chr = key_1L_chr)
  predictors_tb <- mdl_meta_data_ls$predictors_lup
  if(!is.null(mdl_nm_1L_chr)){
    predictors_tb <- predictors_tb %>%
      dplyr::filter(short_name_chr %in% (ready4::get_from_lup_obj(mdls_lup,
                                match_value_xx = mdl_nm_1L_chr,
                                match_var_nm_1L_chr = "mdl_nms_chr",
                                target_var_nm_1L_chr = "predrs_ls",
                                evaluate_1L_lgl = F) %>% purrr::flatten_chr()))
  }
  if(outp_is_abbrs_tb){
    predictors_tb <- predictors_tb %>%
      dplyr::select(short_name_chr, long_name_chr) %>%
      dplyr::rename(Variable = short_name_chr,
                    Description = long_name_chr)
  }

  return(predictors_tb)
}
get_ttu_dv_predrs <- function(ttu_dv_dss_tb = NULL,
                              ttu_dv_nms_chr = "TTU",
                              server_1L_chr = "dataverse.harvard.edu",
                              key_1L_chr = NULL){
  if(is.null(ttu_dv_dss_tb))
    ttu_dv_dss_tb <- get_ttu_dv_dss(ttu_dv_nms_chr = ttu_dv_nms_chr,
                                    server_1L_chr = server_1L_chr,
                                    key_1L_chr = NULL)
  predrs_chr <- ttu_dv_dss_tb$predrs_ls %>%
    purrr::flatten_chr() %>%
    unique() %>% sort()
  return(predrs_chr)
}
get_tfmn_from_lup <- function(mdl_nm_1L_chr, mdls_lup = NULL){
  if (is.null(mdls_lup))
    utils::data("mdls_lup", envir = environment())
  tfmn_1L_chr <- ready4::get_from_lup_obj(mdls_lup,
                                             target_var_nm_1L_chr = "tfmn_chr",
                                             match_value_xx = mdl_nm_1L_chr,
                                             match_var_nm_1L_chr = "mdl_nms_chr",
                                             evaluate_1L_lgl = F)

  return(tfmn_1L_chr)
}
get_ttu_ds_smrys <- function(ttu_dv_nm_1L_chr = "TTU",
                             server_1L_chr = "dataverse.harvard.edu",
                             key_1L_chr = NULL,
                             reference_int = NULL){

  ds_ls <- dataverse::dataverse_contents(ttu_dv_nm_1L_chr,
                                         key = key_1L_chr,
                                         server = server_1L_chr)
  ids_chr <- ds_ls %>%
    purrr::map_chr(~.x$persistentUrl)
  dv_dss_mdl_smrys_ls <- get_dv_dss_mdl_smrys(ids_chr,
                                        server_1L_chr = server_1L_chr,
                                        key_1L_chr = key_1L_chr)
  dv_dss_mdl_smrys_ls <- dv_dss_mdl_smrys_ls %>%
    purrr::compact()
  if(!is.null(reference_int) & length(dv_dss_mdl_smrys_ls)>0)
    dv_dss_mdl_smrys_ls <- reference_int %>%
    purrr::map(~dv_dss_mdl_smrys_ls %>%
    purrr::pluck(.x)) %>%
    stats::setNames(names(dv_dss_mdl_smrys_ls)[reference_int])
  return(dv_dss_mdl_smrys_ls)
}
get_ttu_dv_dss <- function(ttu_dv_nms_chr = "TTU",
                           server_1L_chr = "dataverse.harvard.edu",
                           key_1L_chr = NULL){
  dv_dss_mdl_smrys_ls <- ttu_dv_nms_chr %>%
    purrr::map(~get_ttu_ds_smrys(.x,
                                 server_1L_chr = server_1L_chr,
                                 key_1L_chr = key_1L_chr)) %>%
    purrr::flatten()

  if(length(dv_dss_mdl_smrys_ls) > 0){
    ttu_dss_ls <- names(dv_dss_mdl_smrys_ls) %>%
      purrr::map(~ dataverse::get_dataset(.x,
                                          key = key_1L_chr,
                                          server = server_1L_chr)) %>%
      stats::setNames(names(dv_dss_mdl_smrys_ls))
    ttu_dv_dss_tb <- purrr::map_dfr(1:length(ttu_dss_ls),
                 ~ tibble::tibble(reference_int = .x,
                                  utility_chr = ifelse((dv_dss_mdl_smrys_ls %>%
                                    purrr::pluck(.x))$depnt_var_nm_1L_chr =="EQ5D_total_dbl",
                                    "EQ-5D",
                                    ifelse((dv_dss_mdl_smrys_ls %>%
                                              purrr::pluck(.x))$depnt_var_nm_1L_chr =="aqol6d_total_w",
                                           "AQoL-6D",
                                           (dv_dss_mdl_smrys_ls %>%
                                             purrr::pluck(.x))$depnt_var_nm_1L_chr)),
                                  predrs_ls = list(get_predictors_lup(dv_dss_mdl_smrys_ls %>%
                                                                purrr::pluck(.x)) %>%
                                                                dplyr::pull(long_name_chr)),
                                  ds_url = names(ttu_dss_ls)[.x])) %>%
      dplyr::mutate(publication_url = purrr::map_chr(ds_url,
                                                     ~ get_dv_ds_publication(.x,
                                                                             key_1L_chr = key_1L_chr,
                                                                             server_1L_chr = server_1L_chr)))

    ttu_dv_dss_tb <- ready4use::add_labels_from_dictionary(ttu_dv_dss_tb,
                                            dictionary_tb = tibble::tibble(var_nm_chr = names(ttu_dv_dss_tb),
                                                                           var_desc_chr = c("ID","Utility","Predictors","Dataset","Article")))

  }else{
    ttu_dv_dss_tb <- NULL
  }
  return(ttu_dv_dss_tb)
}
