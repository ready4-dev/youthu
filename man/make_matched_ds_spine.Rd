% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_make.R
\name{make_matched_ds_spine}
\alias{make_matched_ds_spine}
\title{Make matched dataset spine}
\usage{
make_matched_ds_spine(
  ds_tb,
  round_var_nm_1L_chr = "Timepoint_chr",
  timepoint_bl_val_1L_chr = "Baseline",
  cmprsn_var_nm_1L_chr = "study_arm_chr",
  active_arm_val_1L_chr = "Intervention",
  id_var_nm_1L_chr = "fkClientID",
  match_on_vars_chr
)
}
\arguments{
\item{ds_tb}{Dataset (a tibble)}

\item{round_var_nm_1L_chr}{Round var name (a character vector of length one), Default: 'Timepoint_chr'}

\item{timepoint_bl_val_1L_chr}{Timepoint bl value (a character vector of length one), Default: 'Baseline'}

\item{cmprsn_var_nm_1L_chr}{Cmprsn var name (a character vector of length one), Default: 'study_arm_chr'}

\item{active_arm_val_1L_chr}{Active arm value (a character vector of length one), Default: 'Intervention'}

\item{id_var_nm_1L_chr}{Id var name (a character vector of length one), Default: 'fkClientID'}

\item{match_on_vars_chr}{Match on vars (a character vector)}
}
\value{
Matched dataset (a tibble)
}
\description{
make_matched_ds_spine() is a Make function that creates a new R object. Specifically, this function implements an algorithm to make matched dataset spine. The function returns Matched dataset (a tibble).
}
\keyword{internal}