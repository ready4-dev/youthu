% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fn_add.R
\name{add_qalys}
\alias{add_qalys}
\title{Add Quality Adjusted Life Years}
\usage{
add_qalys(
  ds_tb,
  cmprsn_var_nm_1L_chr = "study_arm_chr",
  duration_var_nm_1L_chr = "duration_prd",
  id_var_nm_1L_chr = "fkClientID",
  match_idx_var_nm_1L_chr = "match_idx_int",
  msrmnt_date_var_nm_1L_chr = "date_dtm",
  qalys_var_nm_1L_chr = "qalys_dbl",
  round_var_nm_1L_chr = "round",
  round_bl_val_1L_chr = "Baseline",
  utl_change_var_nm_1L_chr = "utl_change_dbl",
  utl_var_nm_1L_chr = "utility_dbl",
  reshape_1L_lgl = T
)
}
\arguments{
\item{ds_tb}{Dataset (a tibble)}

\item{cmprsn_var_nm_1L_chr}{Comparison variable name (a character vector of length one), Default: 'study_arm_chr'}

\item{duration_var_nm_1L_chr}{Duration variable name (a character vector of length one), Default: 'duration_prd'}

\item{id_var_nm_1L_chr}{Identity variable name (a character vector of length one), Default: 'fkClientID'}

\item{match_idx_var_nm_1L_chr}{Match index variable name (a character vector of length one), Default: 'match_idx_int'}

\item{msrmnt_date_var_nm_1L_chr}{Measurement date variable name (a character vector of length one), Default: 'date_dtm'}

\item{qalys_var_nm_1L_chr}{Quality Adjusted Life Years variable name (a character vector of length one), Default: 'qalys_dbl'}

\item{round_var_nm_1L_chr}{Round variable name (a character vector of length one), Default: 'round'}

\item{round_bl_val_1L_chr}{Round baseline value (a character vector of length one), Default: 'Baseline'}

\item{utl_change_var_nm_1L_chr}{Utility change variable name (a character vector of length one), Default: 'utl_change_dbl'}

\item{utl_var_nm_1L_chr}{Utility variable name (a character vector of length one), Default: 'utility_dbl'}

\item{reshape_1L_lgl}{Reshape (a logical vector of length one), Default: T}
}
\value{
Updated dataset (a tibble)
}
\description{
add_qalys() is an Add function that updates an object by adding new values to new or empty fields. Specifically, this function implements an algorithm to add quality adjusted life years. The function returns Updated dataset (a tibble).
}
\keyword{internal}
