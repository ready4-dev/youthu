get_dss_using_predrs <- function(mdl_predrs_in_ds_chr,
                                 ttu_dv_dss_tb = NULL,
                                 ttu_dv_nm_1L_chr = "firstbounce",
                                 server_1L_chr = "dataverse.harvard.edu",
                                 key_1L_chr = NULL){
  if(is.null(ttu_dv_dss_tb))
    ttu_dv_dss_tb <- get_ttu_dv_dss(ttu_dv_nm_1L_chr = ttu_dv_nm_1L_chr,
                                    server_1L_chr = server_1L_chr,
                                    key_1L_chr = NULL)
  if(is.null(ttu_dv_dss_tb))
    ttu_dv_dss_tb <- ttu_dv_dss_tb %>%
      dplyr::filter(predrs_ls %>% purrr::map_lgl(~!identical(intersect(.x,mdl_predrs_in_ds_chr),
                                                             character(0))))
  return(ttu_dv_dss_tb)
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
get_mdl_ctlg_url <- function(mdls_lup,
                             mdl_nm_1L_chr,
                             server_1L_chr = "dataverse.harvard.edu",
                             key_1L_chr = NULL){
  dv_ds_nm_1L_chr <- ready4fun::get_from_lup_obj(mdls_lup,
                              match_value_xx = mdl_nm_1L_chr,
                              match_var_nm_1L_chr = "mdl_nms_chr",
                              target_var_nm_1L_chr = "ds_url",
                              evaluate_lgl = F)
  ds_ls <- dataverse::dataset_files(dv_ds_nm_1L_chr,
                                    server = server_1L_chr,
                                    key = key_1L_chr)
  all_lbls_chr <- purrr::map_chr(ds_ls,~.x$label)
  include_lgl <- all_lbls_chr %>% purrr::map_lgl(~startsWith(.x,"TS_TTU_Mdls_Smry"))
  all_descs_chr <- purrr::map_chr(ds_ls,~.x$description)
  include_lgl <- include_lgl & (all_descs_chr %>% purrr::map_lgl(~stringr::str_detect(.x,
                                                                                     ready4fun::get_from_lup_obj(mdls_lup,
                                                                                                                 match_value_xx = mdl_nm_1L_chr,
                                                                                                                 match_var_nm_1L_chr = "mdl_nms_chr",
                                                                                                                 target_var_nm_1L_chr = "source_chr",
                                                                                                                 evaluate_lgl = F))))
  idx_1L_int <- which(include_lgl)
  if(identical(idx_1L_int,integer(0))){
    ctlg_url <- NULL
  }else{
    ctlg_url <- paste0("https://dataverse.harvard.edu/api/access/datafile/",ds_ls[[idx_1L_int]]$dataFile$id)
  }
  return(ctlg_url)
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
get_mdls_using_predrs <- function(mdl_predrs_in_ds_chr,
                                  ttu_dv_dss_tb = NULL,
                                  ttu_dv_nm_1L_chr = "firstbounce",
                                  server_1L_chr = "dataverse.harvard.edu",
                                  key_1L_chr = NULL){
  if(is.null(ttu_dv_dss_tb))
    ttu_dv_dss_tb <- get_ttu_dv_dss(ttu_dv_nm_1L_chr = ttu_dv_nm_1L_chr,
                                    server_1L_chr = server_1L_chr,
                                    key_1L_chr = NULL)
  if(!is.null(ttu_dv_dss_tb))
    ttu_dv_dss_tb <- get_dss_using_predrs(ttu_dv_dss_tb,
                                          mdl_predrs_in_ds_chr = mdl_predrs_in_ds_chr)
  if(!is.null(ttu_dv_dss_tb)){
    ds_smrys_ls <- get_ttu_ds_smrys("firstbounce", reference_int = ttu_dv_dss_tb$reference)
    mdls_lup <- ds_smrys_ls %>%
      purrr::map2_dfr(names(ds_smrys_ls),
                      ~{
                        predictors_lup <- .x$predictors_lup
                        urls_chr <- .y
                        predrs_short_nms_chr <- mdl_predrs_in_ds_chr %>%
                          purrr::map_chr(~ready4fun::get_from_lup_obj(predictors_lup,
                                                                      match_value_xx = .x,
                                                                      match_var_nm_1L_chr = "long_name_chr",
                                                                      target_var_nm_1L_chr = "short_name_chr",
                                                                      evaluate_lgl = F))
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
get_predictors <- function(ingredients_ls){
  predictors_tb <- ingredients_ls$predictors_lup %>%
    dplyr::select(short_name_chr, long_name_chr) %>%
    dplyr::rename(Variable = short_name_chr,
                  Description = long_name_chr)
  return(predictors_tb)
}
get_ttu_dv_predrs <- function(ttu_dv_dss_tb = NULL,
                              ttu_dv_nm_1L_chr = "firstbounce",
                              server_1L_chr = "dataverse.harvard.edu",
                              key_1L_chr = NULL){
  if(is.null(ttu_dv_dss_tb))
    ttu_dv_dss_tb <- get_ttu_dv_dss(ttu_dv_nm_1L_chr = ttu_dv_nm_1L_chr,
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
  tfmn_1L_chr <- ready4fun::get_from_lup_obj(mdls_lup,
                                             target_var_nm_1L_chr = "tfmn_chr",
                                             match_value_xx = mdl_nm_1L_chr,
                                             match_var_nm_1L_chr = "mdl_nms_chr",
                                             evaluate_lgl = F)

  return(tfmn_1L_chr)
}
get_ttu_ds_smrys <- function(ttu_dv_nm_1L_chr = "firstbounce",
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
  if(!is.null(reference_1L_int) & length(dv_dss_mdl_smrys_ls)>0)
    dv_dss_mdl_smrys_ls <- reference_int %>%
    purrr::map(~dv_dss_mdl_smrys_ls %>%
    purrr::pluck(.x)) %>%
    stats::setNames(names(dv_dss_mdl_smrys_ls)[reference_int])
  return(dv_dss_mdl_smrys_ls)
}
get_ttu_dv_dss <- function(ttu_dv_nm_1L_chr = "firstbounce",
                           server_1L_chr = "dataverse.harvard.edu",
                           key_1L_chr = NULL){
  dv_dss_mdl_smrys_ls <- get_ttu_ds_smrys(ttu_dv_nm_1L_chr,
                                    server_1L_chr = server_1L_chr,
                                    key_1L_chr = key_1L_chr)

  if(length(dv_dss_mdl_smrys_ls) > 0){
    ttu_dss_ls <- names(dv_dss_mdl_smrys_ls) %>%
      purrr::map(~ dataverse::get_dataset(.x,
                                          key = key_1L_chr,
                                          server = server_1L_chr)) %>%
      stats::setNames(names(dv_dss_mdl_smrys_ls))
    ttu_dv_dss_tb <- purrr::map_dfr(1:length(ttu_dss_ls),
                 ~ tibble::tibble(reference = .x,
                                  title = (ttu_dss_ls %>%
                                             purrr::pluck(.x))$metadataBlocks$citation$fields %>%
                                  ready4fun::get_from_lup_obj(match_value_xx = "title",
                                                              match_var_nm_1L_chr = "typeName",
                                                              target_var_nm_1L_chr = "value",
                                                              evaluate_lgl = F) %>%
                                  purrr::flatten_chr(),
                                  predrs_ls = list(get_predictors(dv_dss_mdl_smrys_ls %>%
                                                                purrr::pluck(.x)) %>%
                                                                dplyr::pull(Description)),
                                  id = names(ttu_dss_ls)[.x],),
        )
  }else{
    ttu_dv_dss_tb <- NULL
  }
  return(ttu_dv_dss_tb)
}
