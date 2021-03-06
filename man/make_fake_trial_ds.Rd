% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_make.R
\name{make_fake_trial_ds}
\alias{make_fake_trial_ds}
\title{Make fake trial dataset}
\usage{
make_fake_trial_ds(
  ds_tb,
  id_var_nm_1L_chr = "fkClientID",
  round_var_nm_1L_chr = "round",
  round_lvls_chr = c("Baseline", "Follow-up"),
  match_on_vars_chr,
  cmprsn_var_nm_1L_chr = "study_arm_chr",
  cmprsn_groups_chr = c("Intervention", "Control"),
  fns_ls,
  var_nms_chr,
  abs_mean_diff_dbl,
  diff_sd_dbl,
  multiplier_dbl,
  min_dbl,
  max_dbl,
  integer_lgl,
  match_idx_var_nm_1L_chr = "match_idx_int"
)
}
\arguments{
\item{ds_tb}{Dataset (a tibble)}

\item{id_var_nm_1L_chr}{Identity variable name (a character vector of length one), Default: 'fkClientID'}

\item{round_var_nm_1L_chr}{Round variable name (a character vector of length one), Default: 'round'}

\item{round_lvls_chr}{Round levels (a character vector), Default: c("Baseline", "Follow-up")}

\item{match_on_vars_chr}{Match on variables (a character vector)}

\item{cmprsn_var_nm_1L_chr}{Comparison variable name (a character vector of length one), Default: 'study_arm_chr'}

\item{cmprsn_groups_chr}{Comparison groups (a character vector), Default: c("Intervention", "Control")}

\item{fns_ls}{Functions (a list)}

\item{var_nms_chr}{Variable names (a character vector)}

\item{abs_mean_diff_dbl}{Absolute mean difference (a double vector)}

\item{diff_sd_dbl}{Difference standard deviation (a double vector)}

\item{multiplier_dbl}{Multiplier (a double vector)}

\item{min_dbl}{Minimum (a double vector)}

\item{max_dbl}{Maximum (a double vector)}

\item{integer_lgl}{Integer (a logical vector)}

\item{match_idx_var_nm_1L_chr}{Match index variable name (a character vector of length one), Default: 'match_idx_int'}
}
\value{
Updated dataset (a tibble)
}
\description{
make_fake_trial_ds() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make fake trial dataset. The function returns Updated dataset (a tibble).
}
\keyword{internal}
