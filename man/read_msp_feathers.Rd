% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msp.R
\name{read_msp_feathers}
\alias{read_msp_feathers}
\title{Read .msp Feathers}
\usage{
read_msp_feathers(savepath, file_id, feather_type = "all")
}
\arguments{
\item{savepath}{Path to directory containing msp files}

\item{file_id}{String name of previously parsed file}

\item{feather_type}{which type of feather to read: both, attr, or frag}
}
\value{
A list containing compound_attributes and
  compound_fragmentationData tibbles.
}
\description{
Read attributes and fragmentation tables produced by
  \code{\link{run_msp_parse}}.
}
