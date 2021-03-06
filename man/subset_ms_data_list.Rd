% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standard_db_construction_utils.R
\name{subset_ms_data_list}
\alias{subset_ms_data_list}
\title{Subset MS Data List}
\usage{
subset_ms_data_list(ms_data_list, table, filter_call, join_type = "semi_join")
}
\arguments{
\item{ms_data_list}{a list containg mass spectrometry experimental data.
  Three entires are required: ions, samples and measurements:
\itemize{
  \item{ions}
  \item{samples}
  \item{elutions [optional]}
  \item{fragmentations [optional]}
  \item{fragmentationData [optional]}
}}

\item{table}{the table to filter based upon}

\item{filter_call}{the conditions on which to filter (contained within a quosure). e.g., filter_call = quo(stringr::str_detect(compoundName, '^NAD'))}

\item{join_type}{either semi_join or anti_join}
}
\description{
Find the .entry_ids associated with a filtering condition and reduce all tables to that subset
}
