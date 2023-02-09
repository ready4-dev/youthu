library(generics)
library(ready4)
library(ready4show)
library(ready4use)
library(youthvars)
library(scorz)
library(specific)
ready4fun::write_fn_type_dirs()
# MANUAL STEP. Write all your functions to R files in the new "fns" directory.
fns_env_ls <- ready4fun::read_fns(c("data-raw/fns/","data-raw/mthds/"),
                                  fns_env = new.env(parent = globalenv()))
x <- ready4fun::make_pkg_desc_ls(pkg_title_1L_chr = "Transform Youth Outcomes to Health Utility Predictions With Ready4",
                                 pkg_desc_1L_chr = "Tools for mapping measures routinely collected in youth mental health services to Quality Adjusted Life Years (QALYs). Part of the ready4 youth mental health systems model (https://www.ready4-dev.com/).
  This development version of the youthu package has been made available as part of the process of testing and documenting the package. If you have any questions, please contact the authors (matthew.hamilton@orygen.org.au).",
                                 authors_prsn = c(utils::person(given = "Matthew",family = "Hamilton",email = "matthew.hamilton@orygen.org.au", role = c("aut", "cre"),comment = c(ORCID = "0000-0001-7407-9194")),
                                                  utils::person(given = "Caroline",family = "Gao",email = "caroline.gao@orygen.org.au", role = c("aut"),comment = c(ORCID = "0000-0002-0987-2759")),
                                                  utils::person("Orygen", role = c("cph", "fnd")),
                                                  utils::person("Headspace", role = c( "fnd")),
                                                  utils::person("National Health and Medical Research Council", role = c( "fnd"))),
                                 urls_chr = c("https://ready4-dev.github.io/youthu/",
                                              "https://github.com/ready4-dev/youthu",
                                              "https://www.ready4-dev.com/")) %>%
  ready4fun::make_manifest(addl_pkgs_ls = ready4fun::make_addl_pkgs_ls(#depends_chr = "TTU",#c("eq5d","ggfortify"),
                                                                       suggests_chr = c("knitr","rmarkdown")),
  build_ignore_ls = ready4fun::make_build_ignore_ls(file_nms_chr = c("initial_setup.R")),
  check_type_1L_chr = "ready4",
  copyright_holders_chr = "Orygen",
  custom_dmt_ls = ready4fun::make_custom_dmt_ls(user_manual_fns_chr = c("add_utl_predn",
                                                                        "add_qalys_to_ds",
                                                                        "get_mdl_from_dv",
                                                                        "get_mdls_using_predrs",
                                                                        "get_tfmn_from_lup",
                                                                        "make_hlth_ec_smry",
                                                                        "make_sngl_grp_ds",
                                                                        "make_matched_ds",
                                                                        "rename_from_nmd_vec")),##
  dev_pkgs_chr = c("ready4",#"ready4fun",
                   "ready4use","ready4show",
                   "youthvars", "specific"),
  lifecycle_stage_1L_chr = "experimental",
  path_to_pkg_logo_1L_chr = "../../../../../Documentation/Images/youthu-logo/default.png",
  piggyback_to_1L_chr = "ready4-dev/ready4",
  ready4_type_1L_chr = "prediction",
  zenodo_badge_1L_chr = "[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5646668.svg)](https://doi.org/10.5281/zenodo.5646668)")
mdl_ingredients_ls <- ready4::get_rds_from_dv("mdl_ingredients",
                                              dv_ds_nm_1L_chr = "https://doi.org/10.7910/DVN/DKDIB0")
datasets_ls <- list(
  mdl_ingredients_ls$mdls_lup %>%
    ready4fun::make_pkg_ds_ls(db_1L_chr = "mdls_lup",
                              title_1L_chr = "Lookup table of prediction models",
                              desc_1L_chr = "A summary of the key descriptive features of the prediction models included in the youthu package."))
z <- ready4pack::make_pt_ready4pack_manifest(x,
                                             #constructor_r3 = x_ready4class_constructor,
                                             pkg_ds_ls_ls = datasets_ls) %>%
  ready4pack::ready4pack_manifest()
z <- ready4::author(z)
ready4::write_extra_pkgs_to_actions()
devtools::build_vignettes()
# ready4::write_citation_cff(packageDescription("youthu"),
#                            citation_chr = readLines("inst/CITATION"))
# usethis::use_dev_package("TTU",
#                          type = "Depends",
#                          remote = "ready4-dev/TTU")
# usethis::use_dev_package("specific",
#                          remote = "ready4-dev/specific")
# usethis::use_package("truncnorm")

