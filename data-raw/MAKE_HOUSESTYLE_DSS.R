library(magrittr)
ready4fun::read_fns("data-raw/fns")
# source("data-raw/MAKE_CLASSES.R")
abbreviations_lup <- ready4fun::get_rds_from_dv("abbreviations_lup")
fn_type_lup_tb <- ready4fun::get_rds_from_dv("fn_type_lup_tb")
# Edits go here
# Push updates to dataverse
abbreviations_lup %>%
  ready4use::write_paired_ds_fls_to_dv(fl_nm_1L_chr = "abbreviations_lup",
                            desc_1L_chr = "Abbreviations lookup table")
fn_type_lup_tb %>%
  ready4use::write_paired_ds_fls_to_dv(fl_nm_1L_chr = "fn_type_lup_tb",
                            desc_1L_chr = "Function type lookup table")
# NOTE: Don't forget to review and publish the updated dataset.
# Previous edits
# fn_type_lup_tb <- fn_type_lup_tb %>%
#   ready4fun::add_rows_to_fn_type_lup(fn_type_nm_chr = ready4fun::get_new_fn_types(abbreviations_lup = abbreviations_lup,
#                                                                                   fn_type_lup_tb = fn_type_lup_tb),
#                                      fn_type_desc_chr = c("Extracts data from an object.",
#                                                           "Renames elements of an object based on a pre-speccified schema."),
#                                      is_generic_lgl = F,
#                                      is_method_lgl = F)

