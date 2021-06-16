params <-
list(output_type_1L_chr = "HTML")

## ----message=FALSE------------------------------------------------------------
library(magrittr)
library(youthu)

## ----results='hide'-----------------------------------------------------------
data("predictors_lup", package = "youthvars")
predictors_lup[,1:5]

## ----predrstb, echo = F, tab.cap='Predictors used in youthu models', tab.id = 'predrstb', results="asis"----
predictors_lup[,1:5] %>%
  ready4show::print_table(output_type_1L_chr = params$output_type_1L_chr,
                          caption_1L_chr = knitr::opts_current$get("tab.cap"),
                          mkdn_tbl_ref_1L_chr = paste0("tab:",knitr::opts_current$get("tab.id")),
                          add_to_row_ls = NULL) 

## -----------------------------------------------------------------------------
data_tb <- youthvars::replication_popl_tb %>%
  youthvars::transform_raw_ds_for_analysis() %>%
  dplyr::select(fkClientID,round,PHQ9, SOFAS) %>% #
  dplyr::arrange(fkClientID) %>%
  na.omit() %>%
  dplyr::rename(UID = fkClientID,
                Timepoint = round,
                PHQ_total = PHQ9,
                SOFAS_total = SOFAS) %>%
  dplyr::mutate(SOFAS_total = as.integer(round(SOFAS_total,0))) %>%
  tibble::as_tibble()

## ----results='hide'-----------------------------------------------------------
data_tb %>% head()

## ----dstb, echo = F, tab.cap='Illustrative example of a prediction dataset', tab.id = 'dstb', results="asis"----
data_tb %>% 
  head() %>%
  ready4show::print_table(output_type_1L_chr = params$output_type_1L_chr,
                          caption_1L_chr = knitr::opts_current$get("tab.cap"),
                          mkdn_tbl_ref_1L_chr = paste0("tab:",knitr::opts_current$get("tab.id")),
                          add_to_row_ls = NULL) 

## -----------------------------------------------------------------------------
mdl_predrs_in_ds_chr <- predictors_lup$short_name_chr[c(5,7)]
mdl_predrs_in_ds_chr

## -----------------------------------------------------------------------------
summary(data_tb$PHQ_total)

## -----------------------------------------------------------------------------
class(data_tb$PHQ_total)

## -----------------------------------------------------------------------------
summary(data_tb$SOFAS_total)

## -----------------------------------------------------------------------------
class(data_tb$SOFAS_total)

## -----------------------------------------------------------------------------
phq9_test <- youthvars::youthvars_phq9(data_tb$PHQ_total) 
phq9_test %>% class()

## -----------------------------------------------------------------------------
sofas_test <- youthvars::youthvars_sofas(data_tb$SOFAS_total) 
sofas_test %>% class()

## -----------------------------------------------------------------------------
data_tb$UID %>% class()

## -----------------------------------------------------------------------------
data_tb$Timepoint %>% class()
levels(data_tb$Timepoint)

## ----message=FALSE------------------------------------------------------------
data("mdls_lup", package ="youthu")
candidate_mdls_tb <- get_mdls_using_predrs(mdl_predrs_in_ds_chr)

## ----results='hide'-----------------------------------------------------------
candidate_mdls_tb

## ----candmdlstb, echo = F, tab.cap='youthu models that can be applied to the prediction dataset', tab.id = 'candmdlstb', results="asis"----
candidate_mdls_tb %>%
  ready4show::print_table(output_type_1L_chr = params$output_type_1L_chr,
                          caption_1L_chr = knitr::opts_current$get("tab.cap"),
                          mkdn_tbl_ref_1L_chr = paste0("tab:",knitr::opts_current$get("tab.id")),
                          add_to_row_ls = NULL)

## ----warning=FALSE------------------------------------------------------------
model_mdl <- get_mdl_from_dv(candidate_mdls_tb$mdl_nms_chr[4])

## -----------------------------------------------------------------------------
model_mdl 

## ----include=F----------------------------------------------------------------
# TTU::make_brms_mdl_print_ls(mdl_ls = model_mdl,
#                                             label_stub_1L_chr = 5,
#                                             caption_1L_chr = "ANOTHER TEST",
#                                             output_type_1L_chr = "HTML")

