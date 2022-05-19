#' Subset MS Data List
#'
#' Find the .entry_ids associated with a filtering condition and reduce all tables to that subset
#'
#' @inheritParams standard_db_populate_data
#' @param table the table to filter based upon
#' @param filter_call the conditions on which to filter (contained within a quosure). e.g., filter_call = quo(stringr::str_detect(compoundName, '^NAD'))
#' @param join_type either semi_join or anti_join
#'
#' @export
subset_ms_data_list <- function(ms_data_list, table, filter_call, join_type = "semi_join") {
  stopifnot(class(ms_data_list) == "list")
  stopifnot(length(table) == 1, class(table) == "character")
  stopifnot("quosure" %in% class(filter_call))
  stopifnot(length(join_type) == 1, class(join_type) == "character", join_type %in% c("semi_join", "anti_join"))

  filter_matches <- ms_data_list[[table]] %>%
    dplyr::filter(!!filter_call)

  purrr::map(ms_data_list, function(x) {
    do.call(semi_join, list(
      x = x,
      y = filter_matches,
      by = ".entry_id"
    ))
  })
}

#' SQL class lookup table
#'
#' @inheritParams DBI::dbDataType
#'
#' @return a named vector mapping from R classes (names) to SQL classses (values).
#' @examples
#' dbObj <- RMySQL::MySQL()
sql_class_lookup_table <- function(dbObj) {
  R_classes <- list(
    "character" = "bar",
    "integer" = 1L,
    "numeric" = 1.2,
    "logical" = TRUE,
    "Date" = lubridate::mdy("09262017"),
    "POSIXct" = lubridate::mdy_hms("09262017 10:55:00")
  )

  R_classes %>%
    purrr::map(., DBI::dbDataType, dbObj = dbObj) %>%
    unlist()
}

#' Add SQL Table
#'
#' @param src_sql_object An SQL database connection.
#' @param Rsql_class_maps A mapping between R and SQL classes generated from \link{sql_class_lookup_table}.
#' @param table_name The table name to generate.
#' @param indicies NA or a character vector of indicies to create.
#' @param variable_Rtypes A named vector where names are variables and values are their R class.
#' @param overwrite If TRUE overwrite existing tables, if FALSE only create tables which are absent.
#'
#' @return Adds a table to src_sql_object and silently returns 1
#'
#' @examples
#' \dontrun{
#' src_sql_object <- DBI::dbConnect(RMySQL::MySQL(), user = "admin", password = .rs.askForPassword("Database Password:"), dbname = "mass_spec_standards", host = "104.196.252.153")
#' Rsql_class_maps <- sql_class_lookup_table(RMySQL::MySQL())
#' table_name <- "compounds"
#' variable_Rtypes <- c("compoundId" = "integer", "compoundName" = "character", "mass" = "numeric", "formula" = "character", "inchikey" = "character")
#' }
add_sql_table <- function(src_sql_object, Rsql_class_maps, table_name, indicies = NA, variable_Rtypes, overwrite = FALSE) {

  # test inputs
  stopifnot(class(table_name) == "character", length(table_name) == 1)
  if (any(!(variable_Rtypes %in% names(Rsql_class_maps)))) {
    stop(paste(setdiff(variable_Rtypes, names(Rsql_class_maps)), collapse = " & "), " classes are not are not defined in Rsql_class_maps; add these classes to clamdb::sql_class_lookup_table()\n")
  }
  stopifnot(class(overwrite) == "logical", length(overwrite) == 1, overwrite %in% c(TRUE, FALSE))
  stopifnot(length(indicies) == 1 & is.na(indicies) | class(indicies) == "character" & !any(is.na(indicies)))
  if (!any(is.na(indicies))) {
    stopifnot(length(indicies) == length(unique(indicies)))
    if (!all(indicies %in% names(variable_Rtypes))) {
      stop(paste(setdiff(indicies, names(variable_Rtypes)), collapse = " & "), " are provided indicies which are not present as variables; add these variables or remove the indicies")
    }
  }

  # reinitialize existing tables if overwrite is TRUE

  table_exists <- table_name %in% DBI::dbListTables(src_sql_object)
  if (table_exists & overwrite == FALSE) {
    warning(table_name, " already found in connected SQL database, this table will not be reinitialized - set overwrite to TRUE to force initialization\n")
    return(invisible(0))
  } else if (table_exists & overwrite == TRUE) {
    DBI::dbRemoveTable(src_sql_object, table_name)
  }

  table_vars <- plyr::revalue(variable_Rtypes, Rsql_class_maps, warn_missing = FALSE)

  # generate table

  DBI::dbExecute(src_sql_object, DBI::sqlCreateTable(src_sql_object, table_name, table_vars))

  # add indicies

  if (!any(is.na(indicies))) {
    for (an_index in indicies) {
      DBI::dbSendQuery(src_sql_object, sprintf("CREATE INDEX %s ON %s(%s)", paste(table_name, an_index, sep = "_"), table_name, an_index))
    }
  }

  return(invisible(0))
}

standard_db_test_ms_data_list <- function(ms_data_list, strict = TRUE) {

  # test for required tables
  required_tables <- c("ions", "samples")
  ms_data_list_tables <- names(ms_data_list)

  if (length(setdiff(required_tables, ms_data_list_tables)) != 0) {
    stop("missing required tables from ms_data_list: ", paste(setdiff(required_tables, ms_data_list_tables), collapse = " & "))
  }

  possible_tables <- c(required_tables, "elutions", "fragmentations", "fragmentationData")
  if (length(setdiff(ms_data_list_tables, possible_tables)) != 0) {
    stop(paste(setdiff(ms_data_list_tables, possible_tables), collapse = " & "), " are not recognized entries of ms_data_list. Valid entries are: ", paste(possible_tables, collapse = ", "))
  }

  # test that either an elutions table or fragmentation + fragmentationData tables are present

  elution_data_present <- ifelse("elutions" %in% ms_data_list_tables, TRUE, FALSE)
  fragmentation_data_present <- ifelse(all(c("fragmentations", "fragmentationData") %in% ms_data_list_tables), TRUE, FALSE)

  if (!elution_data_present && !fragmentation_data_present) {
    stop("measurement tables were not found: either supply elution data (elutions table), fragmentation data (fragmentations and fragmentationData), or both")
  }

  # Test ms_data_list formatting
  entry_id_key_found <- purrr::map_lgl(ms_data_list, function(x) {
    ".entry_id" %in% colnames(x)
  })
  if (any(!entry_id_key_found)) {
    stop(paste(names(entry_id_key_found)[!entry_id_key_found], collapse = " & "), " tables do not contain a .entry variable:\n  this variable is required in all tables since it is used to match samples and ions to their measurements.")
  }

  data_frames_found <- purrr::map_lgl(ms_data_list, function(x) {
    "data.frame" %in% class(x)
  })
  if (any(!data_frames_found)) {
    stop(paste(names(data_frames_found)[!data_frames_found], collapse = " & "), " tables are not data.frames (or tibbles); all entries in the ms_data_list must be tables.")
  }

  # Test variable classes
  standard_db_test_ms_data_list_classes(ms_data_list)

  if (!all(ms_data_list$samples$mode %in% c("negative", "positive"))) {
    stop("some samples had a \"mode\" that was not negative or postiive")
  }

  # Check that ions are
  ms_data_ions <- ms_data_list$ions

  if (!all(ms_data_ions$mode %in% c("negative", "positive"))) {
    stop("some ions had a \"mode\" that was not negative or postiive")
  }

  invalid_charges <- ms_data_ions %>%
    dplyr::count(precursorCharge, mode) %>%
    dplyr::filter(
      mode == "negative" & precursorCharge > 0 |
        mode == "positive" & precursorCharge < 0
    )

  if (nrow(invalid_charges) > 0) {
    msg <- glue::glue(
      "{sum(invalid_charges$n)} entries had invalid charges:
        a charge > 0 in negative mode or a charge < 0 in positive mode
      "
    )

    if (strict) {
      stop(msg)
    } else {
      warning(msg)
    }
  }
}

standard_db_test_ms_data_list_classes <- function(ms_data_list) {
  class_specifications <- ms_data_list_class_specification()

  specification_errors <- tibble::tibble(tbl = names(ms_data_list), ms_tbls = ms_data_list) %>%
    dplyr::left_join(class_specifications, by = "tbl") %>%
    dplyr::mutate(errors = purrr::map2(ms_tbls, specification_tbls, ms_data_list_class_specification_evaluation)) %>%
    tidyr::unnest(errors) %>%
    dplyr::mutate(message = glue::glue("{tbl} : {errors}"))

  if (nrow(specification_errors) != 0) {
    stop(glue::glue_collapse(specification_errors$message, sep = "\n"))
  }

  return(invisible(0))
}

ms_data_list_class_specification_evaluation <- function(ms_tbl, specification_tbl) {
  ms_tbl_format <- tibble::tibble(
    variable = colnames(ms_tbl),
    ms_class = purrr::map_chr(ms_tbl, class)
  )

  # check classes for provided vars
  wrong_classes <- ms_tbl_format %>%
    dplyr::inner_join(specification_tbl, by = "variable") %>%
    dplyr::filter(ms_class != required_class) %>%
    {
      glue::glue("{.$variable} is {.$ms_class} and should be {.$required_class}")
    }

  missing_required_vars <- specification_tbl %>%
    dplyr::filter(required) %>%
    dplyr::anti_join(ms_tbl_format, by = "variable") %>%
    {
      glue::glue("{.$variable} is missing")
    }

  return(c(wrong_classes, missing_required_vars))
}

ms_data_list_class_specification <- function() {
  ions_specification <- tibble::tribble(
    ~variable, ~required_class, ~required,
    ".entry_id", "integer", TRUE,
    "compoundName", "character", TRUE,
    "inchi", "character", FALSE,
    "inchikey", "character", TRUE,
    "smiles", "character", FALSE,
    "formula", "character", FALSE,
    "adductName", "character", TRUE,
    "mass", "numeric", TRUE,
    "precursorMz", "numeric", TRUE,
    "mode", "character", TRUE,
    "precursorCharge", "integer", TRUE
  )

  samples_specification <- tibble::tribble(
    ~variable, ~required_class, ~required,
    ".entry_id", "integer", TRUE,
    "chemicalClass", "character", TRUE,
    "sampleType", "character", TRUE,
    "mode", "character", TRUE,
    "analyticalMethod", "character", FALSE,
    "chromatographicMethod", "character", FALSE,
    "date", "Date", FALSE,
    "path", "character", FALSE
  )

  fragmentations_specification <- tibble::tribble(
    ~variable, ~required_class, ~required,
    ".entry_id", "integer", TRUE,
    "energy", "character", FALSE,
    "msLevel", "integer", FALSE
  )


  fragmentationData_specification <- tibble::tribble(
    ~variable, ~required_class, ~required,
    ".entry_id", "integer", TRUE,
    "mz", "numeric", TRUE,
    "ic", "numeric", TRUE,
    "label", "character", FALSE
  )

  elutions_specification <- tibble::tribble(
    ~variable, ~required_class, ~required,
    ".entry_id", "integer", TRUE,
    "rtObserved", "numeric", TRUE,
    "rtAligned", "numeric", FALSE
  )

  specification_nest <- ms_data_specifications <- tibble::tribble(
    ~tbl, ~specification_tbls,
    "ions", ions_specification,
    "samples", samples_specification,
    "fragmentations", fragmentations_specification,
    "fragmentationData", fragmentationData_specification,
    "elutions", elutions_specification
  )

  return(specification_nest)
}

#' Check Standards Database
#'
#' Check the standards for pathologies
#'
#' @inheritParams standard_db_init
#'
#' @return 0 invisibly
check_standards_database <- function(mass_spec_standards_con) {
  compounds <- tbl(mass_spec_standards_con, "compounds") %>%
    dplyr::collect()

  check_standards_database_compounds(compounds)

  return(invisible(0))
}

#' Check Standards Database Compounds
#'
#' @param compounds a compounds table from the standards database
#'
#' @return 0 invisibly
check_standards_database_compounds <- function(compounds) {

  # check for uniqueness of compoundId, compoundName, inchikey, smiles,
  #   and systematicCompoundId

  required_unique_felds <- tibble::tribble(
    ~variable, ~na_okay,
    "compoundId", FALSE,
    "compoundName", FALSE,
    "inchikey", TRUE,
    "smiles", TRUE,
    "systematicCompoundId", TRUE
  )

  purrr::walk2(
    required_unique_felds$variable,
    required_unique_felds$na_okay,
    function(x, y) {
      var_duplicates <- compounds %>%
        dplyr::count(!!rlang::sym(x)) %>%
        dplyr::filter(n > 1)

      if (y) {
        var_duplicates <- var_duplicates %>%
          dplyr::filter(!is.na(!!rlang::sym(x)))
      }

      if (nrow(var_duplicates) != 0) {
        stop(glue::glue("{nrow(var_duplicates)} duplicates of {x} in standards database: {paste(var_duplicates[[x]], collapse = ', ')}"))
      }
    }
  )

  return(invisible(0))
}
