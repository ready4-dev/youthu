params <-
list(output_type_1L_chr = "HTML")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE------------------------------------------------------------
library(youthu)
set.seed(1234)

## -----------------------------------------------------------------------------
data("replication_popl_tb", package = "youthvars")
seed_ds_tb <- replication_popl_tb %>% youthvars::transform_raw_ds_for_analysis() %>% dplyr::filter(fkClientID %in% (replication_popl_tb %>% dplyr::filter(round=="Baseline" & PHQ9<20) %>% dplyr::pull(fkClientID)))

## -----------------------------------------------------------------------------
ds_smry_ls <- list(bl_start_date_dtm = Sys.Date() - lubridate::days(300),
                   bl_end_date_dtm = Sys.Date() - lubridate::days(120),
                   cmprsn_var_nm_1L_chr = "study_arm_chr",
                   cmprsn_groups_chr = c("Intervention","Control"),
                   costs_mean_dbl = c(400,1500),
                   costs_sd_dbl = c(100,220),
                   costs_var_nm_1L_chr = "costs_dbl",
                   date_var_nm_1L_chr = "date_psx",
                   duration_args_ls = list(a = 160, b = 2200, mean = 180, sd = 7),
                   duration_fn = truncnorm::rtruncnorm,
                   id_var_nm_1L_chr = "fkClientID",
                   predr_var_nms = c("PHQ9", "SOFAS"),
                   round_var_nm_1L_chr = "round",
                   round_lvls_chr = c("Baseline","Follow-up"), 
                   utl_var_nm_1L_chr = "AQoL6D_HU")


## -----------------------------------------------------------------------------
sngl_grp_ds_tb <- make_sngl_grp_ds(seed_ds_tb,ds_smry_ls = ds_smry_ls)

## -----------------------------------------------------------------------------
matched_ds_tb <- make_matched_ds(sngl_grp_ds_tb, 
                                 cmprsn_smry_tb = tibble::tibble(var_nms_chr = c(ds_smry_ls$predr_var_nms, ds_smry_ls$costs_var_nm_1L_chr),
                                 fns_ls = list(stats::rnorm,stats::rnorm,stats::rnorm),
                                 abs_mean_diff_dbl = c(2,2,300),
                                 diff_sd_dbl = c(2,2,200),
                                 multiplier_dbl = c(-1,-1,1),
                                 min_dbl = c(0,0,0),
                                 max_dbl = c(27,100,Inf),
                                 integer_lgl = c(T,T,F)), 
                                 ds_smry_ls = ds_smry_ls)

## -----------------------------------------------------------------------------
matched_ds_tb %>% 
  head() %>%
  ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = "Input dataset")

## -----------------------------------------------------------------------------
summary((matched_ds_tb %>% dplyr::filter(study_arm_chr == "Control" & round == "Baseline"))[5:6])

## -----------------------------------------------------------------------------
summary((matched_ds_tb %>% dplyr::filter(study_arm_chr == "Control" & round == "Follow-up"))[5:7])

## -----------------------------------------------------------------------------
summary((matched_ds_tb %>% dplyr::filter(study_arm_chr == "Intervention" & round == "Baseline"))[5:6])

## -----------------------------------------------------------------------------
summary((matched_ds_tb %>% dplyr::filter(study_arm_chr == "Intervention" & round == "Follow-up"))[5:7])

## -----------------------------------------------------------------------------
data("predictors_lup", package = "youthvars")
candidate_mdls_tb <- get_mdls_using_predrs(predictors_lup$short_name_chr[c(5,7)])
model_mdl <- get_mdl_from_dv(candidate_mdls_tb$mdl_nms_chr[4])

