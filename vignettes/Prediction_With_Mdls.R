params <-
list(output_type_1L_chr = "HTML")

## ----message=FALSE------------------------------------------------------------
library(ready4)
library(ready4use)
library(youthu)

## -----------------------------------------------------------------------------
ttu_dv_dss_tb <- get_ttu_dv_dss("TTU")

## ----ttudss, echo = F, tab.cap='Transfer to Utility Datasets', tab.id = 'ttudss', results="asis"----
ttu_dv_dss_tb[,c(1:3,5)] %>% 
  ready4show::print_table(output_type_1L_chr = params$output_type_1L_chr,
                          caption_1L_chr = knitr::opts_current$get("tab.cap"),
                          use_lbls_as_col_nms_1L_lgl = T,
                          mkdn_tbl_ref_1L_chr = paste0("tab:",knitr::opts_current$get("tab.id")),
                          add_to_row_ls = NULL,
                          scroll_box_args_ls = list(width = "100%")) 

## -----------------------------------------------------------------------------
mdls_lup <- get_mdls_lup(ttu_dv_dss_tb = ttu_dv_dss_tb,
                         utility_type_chr = "AQoL-6D",
                         mdl_predrs_in_ds_chr = c("PHQ9 total score",
                                                  "SOFAS total score"))

## ----mdlslupsel, echo = F, tab.cap='Selected elements from Models Look-Up Table', tab.id = 'mdlslupsel', results="asis"----
mdls_lup[,c(1,2,5)] %>% 
  ready4use::add_labels_from_dictionary(dictionary_tb = tibble::tibble(var_nm_chr = names(mdls_lup)[c(1,2,5)],
                                                                       var_desc_chr = c("Catalogue reference", "Predictors", "Analysis"))) %>%
  ready4show::print_table(output_type_1L_chr = params$output_type_1L_chr,
                          caption_1L_chr = knitr::opts_current$get("tab.cap"),
                          use_lbls_as_col_nms_1L_lgl = T,
                          mkdn_tbl_ref_1L_chr = paste0("tab:",knitr::opts_current$get("tab.id")),
                          add_to_row_ls = NULL,
                          scroll_box_args_ls = list(width = "100%")) 

## -----------------------------------------------------------------------------
get_dv_mdl_smrys(mdls_lup,
                 mdl_nms_chr = "PHQ9_SOFAS_1_OLS_CLL")


## ----results='asis'-----------------------------------------------------------
get_mdl_ctlg_url(mdls_lup,
                 mdl_nm_1L_chr = "PHQ9_SOFAS_1_OLS_CLL")

## -----------------------------------------------------------------------------
data_tb <- make_fake_ds_one()

## ----dstb, echo = F, tab.cap='Illustrative example of a prediction dataset', tab.id = 'dstb', results="asis"----
data_tb %>% 
  head() %>%
  ready4show::print_table(output_type_1L_chr = params$output_type_1L_chr,
                          caption_1L_chr = knitr::opts_current$get("tab.cap"),
                          mkdn_tbl_ref_1L_chr = paste0("tab:",knitr::opts_current$get("tab.id")),
                          add_to_row_ls = NULL,
                          scroll_box_args_ls = list(width = "100%")) 

## -----------------------------------------------------------------------------
predictors_lup <- get_predictors_lup(mdls_lup = mdls_lup,
                                     mdl_nm_1L_chr = "PHQ9_SOFAS_1_OLS_CLL")

## ----predluptb, echo = F, tab.cap='Model predictors lookup table', tab.id = 'predluptb', results="asis"----
predictors_lup %>% 
  ready4show::print_table(output_type_1L_chr = params$output_type_1L_chr,
                          caption_1L_chr = knitr::opts_current$get("tab.cap"),
                          mkdn_tbl_ref_1L_chr = paste0("tab:",knitr::opts_current$get("tab.id")),
                          add_to_row_ls = NULL,
                          scroll_box_args_ls = list(width = "100%")) 

## -----------------------------------------------------------------------------
mdl_meta_data_ls <- ingest(Ready4useRepos(gh_repo_1L_chr = "ready4-dev/youthu", gh_tag_1L_chr = "v0.0.0.91125"), fls_to_ingest_chr = c("mdl_meta_data_ls"), metadata_1L_lgl = F)

## -----------------------------------------------------------------------------
predn_ds_ls <- make_predn_metadata_ls(data_tb,
                                      id_var_nm_1L_chr = "UID",
                                      msrmnt_date_var_nm_1L_chr = "Date",
                                      predr_vars_nms_chr = c(PHQ9 = "PHQ_total",SOFAS = "SOFAS_total"),
                                      round_var_nm_1L_chr = "Timepoint",
                                      round_bl_val_1L_chr = "Baseline",
                                      utl_var_nm_1L_chr = "AQoL6D_HU",
                                      mdls_lup = mdls_lup,
                                      mdl_nm_1L_chr = "PHQ9_SOFAS_1_OLS_CLL")

## -----------------------------------------------------------------------------
data_tb <- add_utl_predn(data_tb,
                         predn_ds_ls = predn_ds_ls)


## ----updtb, echo = F, tab.cap='Prediction dataset with predicted utilities', tab.id = 'updtb', results="asis"----
data_tb %>% 
  head() %>%
  ready4show::print_table(output_type_1L_chr = params$output_type_1L_chr,
                          caption_1L_chr = knitr::opts_current$get("tab.cap"),
                          mkdn_tbl_ref_1L_chr = paste0("tab:",knitr::opts_current$get("tab.id")),
                          add_to_row_ls = NULL,
                          scroll_box_args_ls = list(width = "100%")) 

## -----------------------------------------------------------------------------
summary(data_tb$AQoL6D_HU)

## -----------------------------------------------------------------------------
data_tb <- data_tb %>% add_qalys_to_ds(predn_ds_ls = predn_ds_ls,
                                       include_predrs_1L_lgl = F,
                                       reshape_1L_lgl = F)

## ----qalytb, echo = F, tab.cap='Prediction dataset with QALYs', tab.id = 'qalytb', results="asis"----
data_tb %>% 
  head() %>%
  ready4show::print_table(output_type_1L_chr = params$output_type_1L_chr,
                          caption_1L_chr = knitr::opts_current$get("tab.cap"),
                          mkdn_tbl_ref_1L_chr = paste0("tab:",knitr::opts_current$get("tab.id")),
                          add_to_row_ls = NULL,
                          scroll_box_args_ls = list(width = "100%")) 


