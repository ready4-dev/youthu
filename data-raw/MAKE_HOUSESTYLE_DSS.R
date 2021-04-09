library(magrittr)
ready4fun::read_fns("data-raw/fns")
# source("data-raw/MAKE_CLASSES.R")
abbreviations_lup <- ready4fun::get_rds_from_dv("abbreviations_lup")
fn_type_lup_tb <- ready4fun::get_rds_from_dv("fn_type_lup_tb")
object_type_lup <- ready4fun::get_rds_from_dv("object_type_lup")
seed_obj_lup_tb <- ready4fun::get_rds_from_dv("seed_obj_lup_tb")
# Edits go here
seed_obj_lup_tb <- make_obj_lup_spine(seed_obj_lup_tb,
                                      new_entries_tb = tibble::tibble(short_name_chr = c("env","plt"),
                                                                      long_name_chr = c("environment","plot"),
                                                                      atomic_element_lgl = c(rep(F,2)),
                                                                      r3_element_lgl = c(rep(F,2))))
object_type_lup <- make_obj_lup(seed_obj_lup_tb)
obj_type_lup_new_cses_tb <- get_obj_type_lup_new_cses_tb(updated_obj_type_lup_tb = object_type_lup,
                                                         old_obj_type_lup_tb = get_rds_from_dv("object_type_lup"),
                                                         excluded_chr = c("env","plt"))

abbreviations_lup <- abbreviations_lup %>%
  add_lups(new_lup = obj_type_lup_new_cses_tb,
           key_var_nm_1L_chr = "short_name_chr") %>%
  ready4fun::update_abbr_lup(short_name_chr = c("abs",
                                                "cse",
                                                "cst",
                                                "diff",
                                                "dstr",
                                                "efc",
                                                "efcn",
                                                "env",
                                                "fup",
                                                "lvl",
                                                "multpl",
                                                "nmd",
                                                "qaly",
                                                "sd",
                                                "tmpt",
                                                "wtp",
                                                "wtpt"),
                             long_name_chr = c("absolute",
                                               "case",
                                               "cost",
                                               "difference",
                                               "distribution",
                                               "effective",
                                               "effectiveness",
                                               "environment",
                                               "follow-up",
                                               "level",
                                               "multiplier",
                                               "named",
                                               "Quality Adjusted Life Year",
                                               "standard deviation",
                                               "time point",
                                               "willingness to pay",
                                               "willingness to pay threshold"),
                             no_plural_chr = c("effective","effectiveness",
                                               "named","willingness to pay"))
# Push updates to dataverse
seed_obj_lup_tb %>%
  ready4use::write_paired_ds_fls_to_dv(fl_nm_1L_chr = "seed_obj_lup_tb",
                                       desc_1L_chr = "Seed object type lookup table")
object_type_lup %>%
  ready4use::write_paired_ds_fls_to_dv(fl_nm_1L_chr = "object_type_lup",
                                       desc_1L_chr = "Object type lookup table")
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

