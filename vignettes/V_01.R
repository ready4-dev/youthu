## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  # Uncomment and run if installation is required.
#  # utils::install.packages("devtools")
#  # devtools::install_github("ready4-dev/youthu")

## ----message=FALSE------------------------------------------------------------
library(ready4)
library(ready4use)
library(specific)
library(youthu)

## -----------------------------------------------------------------------------
X <- Ready4useRepos(dv_nm_1L_chr = "fakes", dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/HJXYKQ", 
                    dv_server_1L_chr = "dataverse.harvard.edu",
                    gh_repo_1L_chr = "ready4-dev/youthu", gh_tag_1L_chr = "v0.0.0.91125")

## -----------------------------------------------------------------------------
data_tb <- ingest(X, fls_to_ingest_chr = c("ymh_clinical_tb"), metadata_1L_lgl = F)

## -----------------------------------------------------------------------------
data_tb <- data_tb %>% dplyr:: select(c("fkClientID", "round", 
                                        "d_interview_date", 
                                        "gad7_total", "phq9_total")) %>%
  tidyr::pivot_wider(names_from = c("round"), 
                     values_from = c("d_interview_date", "gad7_total", "phq9_total")) %>%
  dplyr::rename_with(~stringr::str_replace(.x,"_Baseline","_t1") %>% 
                       stringr::str_replace("_Follow-up","_t2") %>% 
                       stringr::str_replace("_total",""))

## ----eval=FALSE---------------------------------------------------------------
#  data_tb %>% head()

## ----echo=FALSE---------------------------------------------------------------
data_tb %>% head() %>% ready4show::print_table(output_type_1L_chr = "HTML",
                          caption_1L_chr = knitr::opts_current$get("tab.cap"),
                          use_lbls_as_col_nms_1L_lgl = T,
                          mkdn_tbl_ref_1L_chr = paste0("tab:",knitr::opts_current$get("tab.id")),
                          add_to_row_ls = NULL,
                          scroll_box_args_ls = list(width = "100%")) 

