## This script creates the data files embedded with this package.
## As it involves interaction, it must be run  in four separate parts.
##
## PART ONE
##
# 1. Load magrittr package to that the pipe operator ("%>%") can be used in this script.
library(magrittr)
# 2. Create "fns", "gnrcs" and "mthds" sub-directories.
ready4fun::write_fn_type_dirs()
# 3. MANUAL STEP. Write all your functions to R files in the new "fns" directory.
##
## PART TWO
##
# 4. Set-up package structure
ready4fun::make_pkg_desc_ls(pkg_title_1L_chr = "Youth Outcomes to Health Utility",
                            pkg_desc_1L_chr = "Tools for mapping measures routinely collected in youth mental health services to AQOL 6D Health Utility. Part of the First Bounce model of primary youth mental health services.
  This development version of the youthu package has been made available as part of the process of testing and documenting the package. The tools contained in this development release are designed for use in conjunction with model objects stored in data repositories. The real model objects will be publicly released once the associated scientific manuscript is published. In the mean time, we have included links to placeholder model objects derived from synthetic data.
  For this reason, this release is for demonstration purposes only and this package should not yet be used in analyses deigned to inform policy decisions. If you have any questions, please contact the authors (matthew.hamilton@orygen.org.au).
  The documentation for this package has been automatically generated by the ready4fun package and is therefore quite rudimentary. Human authored documentation will follow in 2021.",
                            authors_prsn = c(utils::person(given = "Matthew",family = "Hamilton",email = "matthew.hamilton@orygen.org.au", role = c("aut", "cre"),comment = c(ORCID = "0000-0001-7407-9194")),
                                              utils::person(given = "Caroline",family = "Gao",email = "caroline.gao@orygen.org.au", role = c("aut"),comment = c(ORCID = "0000-0002-0987-2759")),
                                              utils::person("Orygen", role = c("cph", "fnd")),
                                              utils::person("Headspace", role = c( "fnd")),
                                              utils::person("National Health and Medical Research Council", role = c( "fnd"))),
                            urls_chr = c("https://ready4-dev.github.io/youthu/",
                                         "https://github.com/ready4-dev/youthu",
                                         "https://www.ready4-dev.com/")) %>%
  ready4fun::write_pkg_setup_fls(incr_ver_1L_lgl = F,
                                 delete_r_dir_cnts_1L_lgl = T,
                                 copyright_holders_chr = "Orygen",
                                 check_type_1L_chr = "gh",
                                 path_to_pkg_logo_1L_chr = "../../../../../Documentation/Images/youthu-logo/default.png",
                                 github_repo = "ready4-dev/youthu",
                                 lifecycle_stage_1L_chr = "experimental",
                                 badges_lup = ready4fun::badges_lup,
                                 addl_badges_ls = list(ready4 = "prediction"))
## PAUSE FOR INTERACTION
##
## PART THREE
##
object_type_lup <- ready4fun::get_rds_from_dv("object_type_lup")
pkg_dss_tb <- ready4fun::get_rds_from_dv("abbreviations_lup") %>%
  ready4fun::write_abbr_lup(object_type_lup = object_type_lup)
utils::data("abbreviations_lup")
# 5. Create function types and generics look-up tables
# 5.1 Create a lookup table of function types used in this package and save it as a package dataset (data gets saved in the data directory, documentation script is created in R directory).
pkg_dss_tb <- ready4fun::get_rds_from_dv("fn_type_lup_tb") %>%
  ready4fun::write_dmtd_fn_type_lup(abbreviations_lup = abbreviations_lup,
                                    object_type_lup = object_type_lup,
                                    pkg_dss_tb = pkg_dss_tb)
utils::data("fn_type_lup_tb")
#
# 6. Create a table of all functions to document
fns_dmt_tb <- ready4fun::make_dmt_for_all_fns(paths_ls = ready4fun::make_fn_nms()[1],
                                   undocumented_fns_dir_chr = ready4fun::make_undmtd_fns_dir_chr()[1],
                                   custom_dmt_ls = list(details_ls = NULL,
                                                        inc_for_main_user_lgl_ls = list(force_true_chr = c("add_utl_predn",
                                                                                                           "add_qalys_to_ds",
                                                                                                           "get_mdl_from_dv",
                                                                                                           "get_mdls_using_predrs",
                                                                                                           "get_tfmn_from_lup",
                                                                                                           "make_hlth_ec_smry",
                                                                                                           "make_sngl_grp_ds",
                                                                                                           "make_matched_ds",
                                                                                                           "rename_from_nmd_vec"),
                                                                                                   force_false_chr = NA_character_),
                                                                   args_ls_ls = NULL),
                                   fn_type_lup_tb = fn_type_lup_tb,
                                   abbreviations_lup = abbreviations_lup,
                                   object_type_lup = object_type_lup)
pkg_dss_tb <- fns_dmt_tb %>%
  ready4fun::write_and_doc_ds(db_1L_chr = "fns_dmt_tb",
                              title_1L_chr = "youthu function documentation table",
                              desc_1L_chr = "Meta-data on each youthu function used to create package documentation",
                              url_1L_chr = "https://ready4-dev.github.io/youthu/",
                              abbreviations_lup = abbreviations_lup,
                              object_type_lup = object_type_lup,
                              pkg_dss_tb = pkg_dss_tb)
##
mdls_smry_tb <- ready4use::ready4_dv_import_lup() %>%
  tibble::add_case(data_repo_db_ui_chr = "https://doi.org/10.7910/DVN/JC6PTV", # NOT UP TO DATE
                   file_name_chr = "mdls_smry_tb",
                   file_type_chr = ".csv",
                   data_repo_file_ext_chr = ".tab") %>%
  ready4use::get_data()
utils::data("predictors_lup", package = "youthvars")
pkg_dss_tb <- tibble::tibble(mdl_nms_chr = mdls_smry_tb$Model %>% unique()) %>%
  dplyr::mutate(predrs_ls = mdl_nms_chr %>% strsplit("_") %>% purrr::map(~ .x[.x %in% c(predictors_lup$short_name_chr)]),
                mdl_type_chr = mdl_nms_chr %>% strsplit("_") %>% purrr::map(~ .x[.x %in% c("GLM", "OLS")]) %>% purrr::flatten_chr(),
                tfmn_chr = mdl_nms_chr %>% purrr::map_chr(~stringr::str_sub(.x,start = 1 +  stringi::stri_locate_last(.x,fixed = "_")[1,1] %>% as.vector())))  %>%
  ready4fun::write_and_doc_ds(db_1L_chr = "mdls_lup",
                              title_1L_chr = "Lookup table of prediction models",
                              desc_1L_chr = "A summary of the key descriptive features of the prediction models included in the youthu package.",
                              abbreviations_lup = abbreviations_lup,
                              object_type_lup = object_type_lup,
                              pkg_dss_tb = pkg_dss_tb)
# 7. Save copy of package documentation to online data repo.
# ds_ls <- ready4use::write_pkg_dss_to_dv_ds_csvs(pkg_dss_tb,
#                                                 dv_nm_1L_chr = "ready4models",
#                                                 ds_url_1L_chr = "https://doi.org/10.7910/DVN/RXGPAT",
#                                                 parent_dv_dir_1L_chr = "../../../../../Data/Dataverse",
#                                                 wait_time_in_secs_int = 5L)
# NOTE: NEED TO UPDATE DIR PATH FOR MODELS
## Note files to be rewritten cannot be open in RStudio.
## 8. Document functions.
usethis::use_build_ignore("initial_setup.R")
usethis::use_package("knitrBootstrap")
usethis::use_package("truncnorm")
usethis::use_package("rmarkdown",type = "Suggests")
usethis::use_dev_package("ready4show")
usethis::use_dev_package("ready4use")
usethis::use_dev_package("youthvars")
usethis::use_dev_package("TTU")
readLines(".github/workflows/R-CMD-check.yaml")[-28] %>%
  writeLines(".github/workflows/R-CMD-check.yaml")
ready4fun::write_and_doc_fn_fls(fns_dmt_tb,
                                r_dir_1L_chr = "R",
                                dev_pkgs_chr = c("ready4fun","ready4show","ready4use","youthvars","TTU","dataverse"),
                                update_pkgdown_1L_lgl = T,
                                path_to_dvpr_dmt_dir_1L_chr = "../../../../../Documentation/Code/Developer",
                                path_to_user_dmt_dir_1L_chr = "../../../../../Documentation/Code/User")
# PAUSE FOR INTERACTIVE
# PART FOUR
# data("prototype_lup")
# if(!identical(prototype_lup,ready4fun::get_rds_from_dv("prototype_lup"))){
#   prototype_lup %>%
#     ready4use::write_paired_ds_fls_to_dv(fl_nm_1L_chr = "prototype_lup",
#                                          desc_1L_chr = "Prototypes lookup table")
# }
# devtools::build_vignettes()
#
ready4fun::write_links_for_website(user_manual_url_1L_chr = "https://github.com/ready4-dev/youthu/releases/download/v0.0.0.9065/youthu_user_0.0.0.9065.pdf",
                                   developer_manual_url_1L_chr = "https://github.com/ready4-dev/youthu/releases/download/v0.0.0.9065/youthu_developer_0.0.0.9065.pdf",
                                   project_website_url_1L_chr = "https://www.ready4-dev.com/")
#


