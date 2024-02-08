params <-
list(output_type_1L_chr = "HTML")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE------------------------------------------------------------
library(youthu)
library(ggplot2)
library(ready4use)
set.seed(1234)

## -----------------------------------------------------------------------------
ds_tb <- make_fake_ds_two()

## ----echo = F-----------------------------------------------------------------
ds_tb %>% 
  head() %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "First few records from input dataset",
                          scroll_box_args_ls = list(width = "100%"))

## -----------------------------------------------------------------------------
summary((ds_tb %>% dplyr::filter(study_arm_chr == "Control" & round == "Baseline"))[5:6])

## -----------------------------------------------------------------------------
summary((ds_tb %>% dplyr::filter(study_arm_chr == "Control" & round == "Follow-up"))[5:7])

## -----------------------------------------------------------------------------
summary((ds_tb %>% dplyr::filter(study_arm_chr == "Intervention" & round == "Baseline"))[5:6])

## -----------------------------------------------------------------------------
summary((ds_tb %>% dplyr::filter(study_arm_chr == "Intervention" & round == "Follow-up"))[5:7])

## -----------------------------------------------------------------------------
mdl_meta_data_ls <- ingest(Ready4useRepos(gh_repo_1L_chr = "ready4-dev/youthu", gh_tag_1L_chr = "v0.0.0.91125"), fls_to_ingest_chr = c("mdl_meta_data_ls"), metadata_1L_lgl = F)

## -----------------------------------------------------------------------------
predn_ds_ls <- make_predn_metadata_ls(ds_tb,
                                      cmprsn_groups_chr = c("Intervention", "Control"),
                                      cmprsn_var_nm_1L_chr = "study_arm_chr",
                                      costs_var_nm_1L_chr = "costs_dbl",
                                      id_var_nm_1L_chr = "fkClientID",
                                      mdl_meta_data_ls = mdl_meta_data_ls,
                                      msrmnt_date_var_nm_1L_chr = "date_psx",
                                      round_var_nm_1L_chr = "round",
                                      round_bl_val_1L_chr = "Baseline",
                                      utl_var_nm_1L_chr = "AQoL6D_HU",
                                      mdls_lup = get_mdls_lup(utility_type_chr = "AQoL-6D",
                                                              mdl_predrs_in_ds_chr = c("PHQ9 total score",
                                                                                       "SOFAS total score"),
                                                              ttu_dv_nms_chr = "TTU"),
                                      mdl_nm_1L_chr =  "PHQ9_SOFAS_1_OLS_CLL")

## -----------------------------------------------------------------------------
ds_tb <- add_utl_predn(ds_tb,
                       predn_ds_ls = predn_ds_ls) %>%
  dplyr::select(fkClientID, round, study_arm_chr, date_psx, duration_prd, dplyr::everything())

## -----------------------------------------------------------------------------
ds_tb  <- ds_tb %>% add_qalys_to_ds(predn_ds_ls = predn_ds_ls,
                                    include_predrs_1L_lgl = T,
                                    reshape_1L_lgl = T)

## ----echo = F-----------------------------------------------------------------
ds_tb %>% head() %>%
    ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "First few records from updated dataset with QALYs", 
                          use_lbls_as_col_nms_1L_lgl = T,
                          scroll_box_args_ls = list(width = "100%"))

## -----------------------------------------------------------------------------
he_smry_ls <- ds_tb %>% make_hlth_ec_smry(predn_ds_ls = predn_ds_ls,
                                                 wtp_dbl = 50000,
                                                 bootstrap_iters_1L_int = 1000L)

## ----fig.width=6--------------------------------------------------------------
BCEA::ceplane.plot(he_smry_ls$ce_res_ls, wtp =50000,  graph = "ggplot2", theme = ggplot2::theme_light())

