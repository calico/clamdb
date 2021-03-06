% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standard_db_construction_utils.R
\name{add_sql_table}
\alias{add_sql_table}
\title{Add SQL Table}
\usage{
add_sql_table(
  src_sql_object,
  Rsql_class_maps,
  table_name,
  indicies = NA,
  variable_Rtypes,
  overwrite = FALSE
)
}
\arguments{
\item{src_sql_object}{An SQL database connection.}

\item{Rsql_class_maps}{A mapping between R and SQL classes generated from \link{sql_class_lookup_table}.}

\item{table_name}{The table name to generate.}

\item{indicies}{NA or a character vector of indicies to create.}

\item{variable_Rtypes}{A named vector where names are variables and values are their R class.}

\item{overwrite}{If TRUE overwrite existing tables, if FALSE only create tables which are absent.}
}
\value{
Adds a table to src_sql_object and silently returns 1
}
\description{
Add SQL Table
}
\examples{
\dontrun{
src_sql_object <- DBI::dbConnect(RMySQL::MySQL(), user = "admin", password = .rs.askForPassword("Database Password:"), dbname = "mass_spec_standards", host = "104.196.252.153")
Rsql_class_maps <- sql_class_lookup_table(RMySQL::MySQL())
table_name <- "compounds"
variable_Rtypes <- c("compoundId" = "integer", "compoundName" = "character", "mass" = "numeric", "formula" = "character", "inchikey" = "character")
}
}
