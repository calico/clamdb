% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standard_db_construction.R
\name{standard_db_summarize_compounds}
\alias{standard_db_summarize_compounds}
\title{Summarize Compounds}
\usage{
standard_db_summarize_compounds(
  ms_data_ions,
  mass_spec_standards_con,
  systematic_compounds_con,
  match_var_order = c("inchikey", "name", "alias", "inchikey_noproto",
    "inchikey_connectivity"),
  annotationSource
)
}
\arguments{
\item{ms_data_ions}{ions table from \code{ms_data_list}.}

\item{mass_spec_standards_con}{A connection to an SQL database which stores
experimental mass spectrometry data.}

\item{systematic_compounds_con}{A connection to an SQL databasing housing
names, structure, aliases and systematic identifiers of known compounds.}

\item{match_var_order}{identifier priority order for matching.}

\item{annotationSource}{A unique character string to label the data source.}
}
\value{
\describe{
  \item{entry_compound_dict}{entry_id to compoundId}
  \item{compound_attr}{compoundId attributes - e.g., systematic id,
    structure, name, formula}
  \item{compound_annotation_sources}{compound_annotation_sources -
    annotations associated with a compoundId for a given annotationSource}
  \item{compoundId_conversions}{rewrite databased compoundIds (if they
    fused with another w/ further annotation)}
  \item{standard_inconsistencies}{multiple compoundIds fusing based on
    annotations}
  \item{systematic_inconsistencies}{multiple systematc compounds matching
    to a single compoundId}
  }
}
\description{
Reduce a set of compounds associated with data to unique entries by name,
  structures, and database aliases.
}
\details{
Compounds are split into two categories - those with identifying fields
  (name and/or structure) and those without any identifying fields.
  ions without identifying fields are each stored as a unique entry.
  ions with identifying fields are grouped by name + structure and
  compared to standards and systematic compounds.
}
