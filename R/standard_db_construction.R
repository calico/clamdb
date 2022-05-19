#' Populate Standard DB Initialization
#'
#' Add a set of tables to a standards database which will be samples /
#'   compound independent.
#'
#' @inheritParams standard_db_populate_data
#' @param overwrite overwrite existing tables
#'
#' @return setup standards_db SQL database schema
#'
#' @export
standard_db_init <- function(mass_spec_standards_con, overwrite = FALSE) {
  stopifnot(
    class(overwrite) == "logical",
    length(overwrite) == 1,
    overwrite %in% c(TRUE, FALSE)
  )

  Rsql_class_maps <- sql_class_lookup_table(mass_spec_standards_con)

  # initialize each table with an appropriate SQL class
  # currently we are not generating several tables which are part of the
  # goal schema: methods, ionAttr, and coelutions

  add_sql_table(
    mass_spec_standards_con,
    Rsql_class_maps,
    "compounds",
    indicies = "compoundId",
    c(
      "compoundId" = "integer",
      "compoundName" = "character",
      "mass" = "numeric",
      "formula" = "character",
      "inchikey" = "character",
      "smiles" = "character",
      "systematicCompoundId" = "integer",
      "match_var" = "character"
    ),
    overwrite
  )

  add_sql_table(
    mass_spec_standards_con,
    Rsql_class_maps,
    "compoundAnnotSource",
    indicies = "compoundId",
    c(
      "compoundId" = "integer",
      "annotationSource" = "character",
      "compoundName" = "character",
      "mass" = "numeric",
      "formula" = "character",
      "inchikey" = "character",
      "smiles" = "character"
    ),
    overwrite
  )

  add_sql_table(
    mass_spec_standards_con,
    Rsql_class_maps,
    "ions",
    indicies = c("ionId", "compoundId"),
    c(
      "ionId" = "integer",
      "compoundId" = "integer",
      "mode" = "character",
      "adductName" = "character",
      "precursorMz" = "numeric",
      "precursorCharge" = "integer"
    ),
    overwrite
  )

  add_sql_table(
    mass_spec_standards_con,
    Rsql_class_maps,
    "samples",
    indicies = "sampleId",
    c(
      "sampleId" = "integer",
      "chemicalClass" = "character",
      "sampleType" = "character",
      "annotationSource" = "character",
      "mode" = "character",
      "analyticalMethod" = "character",
      "chromatographicMethod" = "character",
      "date" = "Date",
      "path" = "character"
    ),
    overwrite
  )

  add_sql_table(
    mass_spec_standards_con,
    Rsql_class_maps,
    "elutions",
    indicies = c("elutionId", "ionId", "sampleId"),
    c(
      "elutionId" = "integer",
      "ionId" = "integer",
      "sampleId" = "integer",
      "rtObserved" = "numeric",
      "rtAligned" = "numeric"
    ),
    overwrite
  )

  add_sql_table(
    mass_spec_standards_con,
    Rsql_class_maps,
    "fragmentations",
    indicies = c("fragId", "ionId", "sampleId"),
    c(
      "fragId" = "integer",
      "ionId" = "integer",
      "sampleId" = "integer",
      "msLevel" = "integer",
      "energy" = "character"
    ),
    overwrite
  )

  add_sql_table(
    mass_spec_standards_con,
    Rsql_class_maps,
    "fragmentationData",
    indicies = "fragId",
    c(
      "fragId" = "integer",
      "mz" = "numeric",
      "ic" = "numeric",
      "label" = "character"
    ),
    overwrite
  )

  add_sql_table(
    mass_spec_standards_con,
    Rsql_class_maps,
    "compoundAttributes",
    indicies = "compoundId",
    c(
      "compoundId" = "integer",
      "attributeType" = "character",
      "attributeValue" = "character"
    ),
    overwrite
  )

  return(invisible(0))
}

#' Populate Standard DB Support Tables
#'
#' @inheritParams standard_db_populate_data
#' @param methods_tables a list containing reduced methods and full methods
#' @param standards_table a table of internal and external standards
#' @param overwrite overwrite existing tables
#'
#' @examples
#' \dontrun{
#' methods_tables <- read_Calico_MS_methods()
#' }
#'
#' @export
standard_db_populate_support_tables <- function(mass_spec_standards_con,
                                                methods_tables = NULL,
                                                standards_table = NULL,
                                                overwrite = FALSE) {
  # chromatography methods
  if (any(class(methods_tables) != "NULL")) {
    DBI::dbWriteTable(
      mass_spec_standards_con,
      "methods",
      methods_tables$reduced_methods,
      overwrite = overwrite,
      row.names = FALSE
    )
    DBI::dbWriteTable(
      mass_spec_standards_con,
      "methodsFull",
      methods_tables$full_methods,
      overwrite = overwrite,
      row.names = FALSE
    )
  }

  # external and internal standards
  if (any(class(standards_table) != "NULL")) {
    DBI::dbWriteTable(
      mass_spec_standards_con,
      "qcStandards",
      standards_table,
      overwrite = overwrite,
      row.names = FALSE
    )
  }

  # add compound attributes
  lipid_attributes <- populate_lipid_attributes(mass_spec_standards_con)
  DBI::dbWriteTable(
    mass_spec_standards_con,
    "compoundAttributes",
    lipid_attributes,
    overwrite = overwrite,
    row.names = FALSE
  )

  # adduct types

  return(invisible(0))
}

#' Populate Standard DB Data
#'
#' Add new mass spectrometry standard data present in \code{ms_data_list} to
#'   the \code{mass_spec_standards_con} database respecting existing entries
#'   and looking up unique compounds in \code{ms_data_list}.
#'
#' @details
#' All supplied tables must include a .entry_id identifier which indicates
#'   the relationship between tables in the original dataset. This function
#'   will build mappings between .entry_id and the correponding sampleId,
#'   compoundId, ionId, elutionId and fragmentationId and append new entries
#'   to these tables as necessary.
#'
#' All samples and measurements will be added but only ions and compounds
#'   which are not already found in the database.
#'
#' @param ms_data_list a list containg mass spectrometry experimental data.
#'   Three entires are required: ions, samples and measurements:
#' \itemize{
#'   \item{ions}
#'   \item{samples}
#'   \item{elutions [optional]}
#'   \item{fragmentations [optional]}
#'   \item{fragmentationData [optional]}
#' }
#' @param mass_spec_standards_con A connection to an SQL database which stores
#'   experimental mass spectrometry data.
#' @param systematic_compounds_con A connection to an SQL databasing housing
#'   names, structure, aliases and systematic identifiers of known compounds.
#' @param annotationSource A unique character string to label the data source.
#' @param write_to_db if FALSE then return the entries which would be
#'   appended to the database, if TRUE then append these to
#'   \code{mass_spec_standards_con}.
#' @param strict return warnings in some cases rather than errors.
#' @param match_var_order identifier priority order for matching.
#' 
#' @return a list of tables (if \code{write_to_db} is FALSE) or if
#'   \code{write_to_db} is TRUE: append new table entries to
#'   \code{mass_spec_standards_con} as side-effects and return 1.
#'
#' @seealso \code{\link{standard_db_init}},
#'   \code{\link{standard_db_populate_support_tables}}
#'
#' @export
standard_db_populate_data <- function(ms_data_list,
                                      mass_spec_standards_con,
                                      systematic_compounds_con,
                                      annotationSource,
                                      match_var_order = c(
                                        "inchikey",
                                        "name",
                                        "alias",
                                        "inchikey_noproto",
                                        "inchikey_connectivity"
                                      ),
                                      write_to_db = FALSE,
                                      strict = TRUE) {

  # check formatting
  checkmate::assertClass(ms_data_list, "list")
  checkmate::assertClass(match_var_order, "character")
  checkmate::assertString(annotationSource)
  checkmate::assertLogical(write_to_db, len = 1)
  checkmate::assertLogical(strict, len = 1)

  elution_data_present <- ifelse(
    "elutions" %in% names(ms_data_list),
    TRUE,
    FALSE
  )
  fragmentation_data_present <- ifelse(
    all(c("fragmentations", "fragmentationData") %in% names(ms_data_list)),
    TRUE,
    FALSE
  )

  standard_db_test_ms_data_list(ms_data_list, strict)

  # summarize samples

  sample_summaries <- standard_db_summarize_samples(
    ms_data_list$samples,
    mass_spec_standards_con,
    annotationSource
  )

  # Processing ions and compounds
  # distinct combinations of compound_defining_fields will be unique entries
  #   in the compound database. identifying_compound_defining_fields are used
  #   to search against existing compounds. if no
  #   identifying_compound_defining_fields are present than the compound is
  #   added as a new unknown

  compound_summaries <- standard_db_summarize_compounds(
    ms_data_list$ions,
    mass_spec_standards_con,
    systematic_compounds_con,
    match_var_order,
    annotationSource
  )
  check_standards_database_compounds(compound_summaries$compound_attr)

  # ions are matched to their compound and then defined by [compoundId, adduct,
  #   mode, charge]

  ion_summaries <- standard_db_summarize_ions(
    ms_data_list$ions,
    compound_summaries,
    mass_spec_standards_con
  )

  # Summarizing measurements - elutions &/or fragmentation

  if (elution_data_present) {
    elution_summary_fields <- c("rtObserved", "rtAligned")
    included_elution_summary_fields <- intersect(
      elution_summary_fields,
      colnames(ms_data_list$elutions)
    )
    if (length(included_elution_summary_fields) == 0) {
      stop(
        "Either rtObserved or rtAligned must be included in the elutions table
        if this table is present in ms_data_list"
      )
    }

    max_elutionId <- tbl(
      mass_spec_standards_con,
      sql("SELECT MAX(elutionId) FROM elutions")
    ) %>%
      dplyr::collect() %>%
      unlist() %>%
      unname()
    if (is.na(max_elutionId)) {
      max_elutionId <- 0
    }

    appended_elutions <- ms_data_list$elutions %>%
      dplyr::left_join(sample_summaries$entry_sample_dict, by = ".entry_id") %>%
      dplyr::left_join(ion_summaries$entry_ion_dict, by = ".entry_id") %>%
      dplyr::mutate(elutionId = max_elutionId + 1:dplyr::n()) %>%
      dplyr::select(!!!quos(!!!rlang::syms(c("elutionId", "ionId", "sampleId", included_elution_summary_fields))))
  }

  if (fragmentation_data_present) {
    fragmentation_summary_fields <- c("msLevel", "energy")
    included_fragmentation_summary_fields <- intersect(fragmentation_summary_fields, colnames(ms_data_list$fragmentations))

    fragmentationData_summary_fields <- c("mz", "ic", "label")
    included_fragmentationData_summary_fields <- intersect(fragmentationData_summary_fields, colnames(ms_data_list$fragmentationData))

    max_fragId <- tbl(mass_spec_standards_con, sql("SELECT MAX(fragId) FROM fragmentations")) %>%
      dplyr::collect() %>%
      unlist() %>%
      unname()
    if (is.na(max_fragId)) {
      max_fragId <- 0
    }

    fragmentation_summary <- ms_data_list$fragmentations %>%
      dplyr::left_join(sample_summaries$entry_sample_dict, by = ".entry_id") %>%
      dplyr::left_join(ion_summaries$entry_ion_dict, by = ".entry_id") %>%
      dplyr::mutate(fragId = max_fragId + 1:dplyr::n())

    if (nrow(fragmentation_summary) != nrow(ms_data_list$fragmentations)) {
      stop("fragmentation dicts are misformatted - this is a bug, please submit a GitHub issue...")
    }
    if (any(is.na(fragmentation_summary$sampleId) | is.na(fragmentation_summary$ionId))) {
      stop("some sample_dict entries or ion_dict entries could not be matched to fragmentaiton data - this is a bug, please submit a GitHub issue...")
    }

    appended_fragmentations <- fragmentation_summary %>%
      dplyr::select(!!!quos(!!!rlang::syms(c("fragId", "ionId", "sampleId", included_fragmentation_summary_fields))))

    appended_fragmentationData <- ms_data_list$fragmentationData %>%
      dplyr::left_join(fragmentation_summary %>%
        dplyr::select(.entry_id, fragId), by = ".entry_id") %>%
      dplyr::select(!!!quos(!!!rlang::syms(c("fragId", included_fragmentationData_summary_fields))))

    na_fragIds <- sum(is.na(appended_fragmentationData$fragId))
    if (na_fragIds != 0) {
      stop(na_fragIds, " fragmentationData entries did not contain a valid fragId entry")
    }
  }

  # append summary tables to the mass_spec_standards database

  output_tables <- list(
    samples = sample_summaries$appended_samples,
    compounds = compound_summaries$compound_attr,
    compoundAnnotSource = compound_summaries$compound_annotation_sources,
    ions = ion_summaries$appended_ions
  )

  if (elution_data_present) {
    output_tables$elutions <- appended_elutions
  }
  if (fragmentation_data_present) {
    output_tables$fragmentations <- appended_fragmentations
    output_tables$fragmentationData <- appended_fragmentationData
  }

  if (write_to_db) {

    # adding compounds replaces the past entries for a compoundId if they are present
    # to leave the most current summary of metadata
    standard_db_remove_old_compounds(compound_summaries, mass_spec_standards_con)

    for (a_table in names(output_tables)) {
      an_output_table <- output_tables[[a_table]] %>%
        # convert character entries to unicode to avoid non-permissive characeters
        dplyr::mutate_if(is.character, iconv, from = "UTF-8", to = "UTF-8", sub = "xx")

      # append new entries to SQL database
      DBI::dbWriteTable(mass_spec_standards_con, a_table, an_output_table, append = TRUE, row.names = FALSE)
    }

    # when compounds are renamed, update to new compound name (and similarly for ions which depend on compounds)
    standard_db_convert_compounds(compound_summaries$compoundId_conversions, mass_spec_standards_con)
  } else {
    return(
      append(
        output_tables,
        list(
          standard_inconsistencies = compound_summaries$standard_inconsistencies,
          systematic_inconsistencies = compound_summaries$systematic_inconsistencies
        )
      )
    )
  }

  std_out_list <- list()
  if (nrow(compound_summaries$standard_inconsistencies) != 0) {
    std_out_list$standard_inconsistencies <- compound_summaries$standard_inconsistencies
  }
  if (nrow(compound_summaries$systematic_inconsistencies) != 0) {
    std_out_list$systematic_inconsistencies <- compound_summaries$systematic_inconsistencies
  }

  # check database for problems
  closeout_check <- try(
    check_standards_database(mass_spec_standards_con),
    silent = TRUE
  )

  if ("try-error" %in% class(closeout_check)) {
    warning(
      "After updating the database an error occurred in
      check_standards_database()
    "
    )

    warning(attr(closeout_check, "condition")$message)
  }

  if (length(std_out_list) == 0) {
    return(list())
  } else {
    return(std_out_list)
  }
}

standard_db_summarize_samples <- function(ms_data_samples, mass_spec_standards_con, annotationSource) {

  # Processing samples
  # test samples
  sample_defining_fields <- c("chemicalClass", "sampleType", "methodsId", "mode", "analyticalMethod", "chromatographicMethod", "date", "path")
  included_sample_defining_fields <- intersect(colnames(ms_data_samples), sample_defining_fields)

  # determine the largest sampleId in the dataset currently
  max_sampleId <- dplyr::tbl(mass_spec_standards_con, sql("SELECT MAX(sampleId) FROM samples")) %>%
    dplyr::collect() %>%
    unlist() %>%
    unname()
  if (is.na(max_sampleId)) {
    max_sampleId <- 0
  }

  # test for validity of annotation source
  past_annotationSources <- tbl(mass_spec_standards_con, sql("SELECT DISTINCT(annotationSource) FROM samples")) %>%
    dplyr::collect() %>%
    unlist() %>%
    unname()
  if (annotationSource %in% past_annotationSources) {
    stop(annotationSource, " is already defined in annotationSource; choose a new unique identifer for this annotation source")
  }

  # find unique samples and autoincrement sampleId
  appended_samples <- ms_data_samples %>%
    dplyr::distinct(!!!quos(!!!rlang::syms(included_sample_defining_fields))) %>%
    dplyr::mutate(
      sampleId = 1:dplyr::n() + max_sampleId,
      annotationSource = annotationSource
    )

  # match .entry_ids to sampleIds
  entry_sample_dict <- ms_data_samples %>%
    dplyr::left_join(appended_samples, by = included_sample_defining_fields) %>%
    dplyr::select(.entry_id, sampleId)

  missing_dict_entries <- length(union(entry_sample_dict$.entry_id, ms_data_samples$.entry_id)) - length(intersect(entry_sample_dict$.entry_id, ms_data_samples$.entry_id))
  if (missing_dict_entries != 0) {
    stop(missing_dict_entries, " .entry_ids could not be matched when matching .entry_ids to sampleIds - this is a bug, please submit a GitHub issue...")
  }

  list(
    entry_sample_dict = entry_sample_dict,
    appended_samples = appended_samples
  )
}

#' Summarize Compounds
#'
#' Reduce a set of compounds associated with data to unique entries by name,
#'   structures, and database aliases.
#'
#' @param ms_data_ions ions table from \code{ms_data_list}.
#' @param match_var_order identifier priority order for matching.
#' @inheritParams standard_db_populate_data
#'
#' @details
#' Compounds are split into two categories - those with identifying fields
#'   (name and/or structure) and those without any identifying fields.
#'   ions without identifying fields are each stored as a unique entry.
#'   ions with identifying fields are grouped by name + structure and
#'   compared to standards and systematic compounds.
#'
#' @return
#' \describe{
#'   \item{entry_compound_dict}{entry_id to compoundId}
#'   \item{compound_attr}{compoundId attributes - e.g., systematic id,
#'     structure, name, formula}
#'   \item{compound_annotation_sources}{compound_annotation_sources -
#'     annotations associated with a compoundId for a given annotationSource}
#'   \item{compoundId_conversions}{rewrite databased compoundIds (if they
#'     fused with another w/ further annotation)}
#'   \item{standard_inconsistencies}{multiple compoundIds fusing based on
#'     annotations}
#'   \item{systematic_inconsistencies}{multiple systematc compounds matching
#'     to a single compoundId}
#'   }
standard_db_summarize_compounds <- function(ms_data_ions,
                                            mass_spec_standards_con,
                                            systematic_compounds_con,
                                            match_var_order = c(
                                              "inchikey",
                                              "name",
                                              "alias",
                                              "inchikey_noproto",
                                              "inchikey_connectivity"
                                            ),
                                            annotationSource) {
  purrr::walk(
    unique(ms_data_ions$mode),
    checkmate::assertChoice,
    choices = c("positive", "negative")
  )

  # evaluate the state of the database for autoincrementing IDs
  max_compoundId <- tbl(
    mass_spec_standards_con,
    sql("SELECT MAX(compoundId) FROM compounds")
  ) %>%
    dplyr::collect() %>%
    unlist() %>%
    unname()
  if (is.na(max_compoundId)) {
    max_compoundId <- 0
  }

  mass_spec_standards_past_compounds <- dplyr::tbl(
    mass_spec_standards_con,
    dbplyr::sql("SELECT compoundName FROM compounds")
  ) %>% dplyr::collect()

  if (nrow(mass_spec_standards_past_compounds) == 0) {
    max_unknown_index <- 0
  } else {
    max_unknown_index <- stringr::str_match(
      mass_spec_standards_past_compounds$compoundName,
      "unknown ([0-9]+)"
    ) %>%
      {
        .[, 2]
      } %>%
      {
        .[!is.na(.)]
      } %>%
      as.integer() %>%
      {
        c(0, .)
      } %>%
      max()
  }

  # defining ions and compounds variables
  # all fields that will be used
  compound_fields <- c("compoundName", "inchikey", "smiles", "formula", "mass")
  included_compound_fields <- intersect(
    colnames(ms_data_ions),
    compound_fields
  )

  # fields which will be used for search
  identifying_compound_fields <- c("compoundName", "inchikey")
  # require that all identifying_compound_defining_fields are included
  missing_identifying_compound_fields <- setdiff(
    identifying_compound_fields,
    colnames(ms_data_ions)
  )

  if (length(missing_identifying_compound_fields) != 0) {
    stop(
      paste(missing_identifying_compound_fields, collapse = " & "),
      " are required fields in the ions table of ms_data_list:
         If you do not have this information, add columns of NAs"
    )
  }

  # define a list which can be passed into and returned by functions to
  #   increment entries
  db_compound_vars <- list(
    max_compoundId = max_compoundId,
    max_unknown_index = max_unknown_index,
    included_compound_fields = included_compound_fields,
    identifying_compound_fields = identifying_compound_fields
  )

  # find ions which do not contain any identifying features
  # these will just be stored as unknowns
  unidentifiable_ions <- dplyr::filter_at(
    ms_data_ions,
    quos(!!!rlang::syms(identifying_compound_fields)),
    dplyr::all_vars(is.na(.))
  )

  unidentifiable_ion_summary <- standard_db_summarize_unidentifiable_compounds(
    unidentifiable_ions,
    db_compound_vars
  )
  db_compound_vars <- unidentifiable_ion_summary$db_compound_vars

  # summarize unique compounds based on name and structure matches
  identifiable_ions <- ms_data_ions %>%
    dplyr::anti_join(unidentifiable_ions, by = ".entry_id")

  identifiable_ion_summary <- standard_db_summarize_identifiable_compounds(
    identifiable_ions,
    mass_spec_standards_con,
    systematic_compounds_con,
    match_var_order,
    db_compound_vars
  )

  entry_compound_dict <- dplyr::bind_rows(
    unidentifiable_ion_summary$entry_compound_dict,
    identifiable_ion_summary$entry_compound_dict
  )

  missing_dict_entries <- length(
    union(entry_compound_dict$.entry_id, ms_data_ions$.entry_id)
  ) - length(
    intersect(entry_compound_dict$.entry_id, ms_data_ions$.entry_id)
  )
  if (missing_dict_entries != 0) {
    stop(
      missing_dict_entries,
      " .entry_ids could not be matched when matching .entry_ids to compoundIds
        this is a bug, please submit a GitHub issue..."
    )
  }
  if (nrow(entry_compound_dict) != nrow(ms_data_ions)) {
    stop("some ions were matched to multiple compounds")
  }
  if (sum(is.na(entry_compound_dict$compoundId)) != 0) {
    stop(
      sum(is.na(new_entry_compound_dict$compoundId)),
      " compound .entry_ids could not be matched to compounds"
    )
  }

  compound_annotation_sources <- entry_compound_dict %>%
    dplyr::left_join(ms_data_ions, by = ".entry_id") %>%
    dplyr::select(!!!rlang::syms(c(
      "compoundId",
      db_compound_vars$included_compound_fields
    ))) %>%
    dplyr::distinct() %>%
    dplyr::mutate(annotationSource = annotationSource)

  list(
    entry_compound_dict = entry_compound_dict,
    compound_attr = dplyr::bind_rows(
      unidentifiable_ion_summary$compound_attr,
      identifiable_ion_summary$compound_attr
    ),
    compound_annotation_sources = compound_annotation_sources,
    compoundId_conversions =
      identifiable_ion_summary$compoundId_conversions,
    standard_inconsistencies =
      identifiable_ion_summary$standard_inconsistencies,
    systematic_inconsistencies =
      identifiable_ion_summary$systematic_inconsistencies
  )
}

#' Standard DB Summarize Unidentifiable Compounds
#'
#' @param unidentifiable_ions filtered ion table from \code{ms_data_list}
#' @param db_compound_vars list of parameters based on state of currently processed standards and database.
standard_db_summarize_unidentifiable_compounds <- function(unidentifiable_ions,
                                                           db_compound_vars) {
  stopifnot(
    class(db_compound_vars) == "list",
    all(
      c("max_unknown_index", "max_compoundId", "included_compound_fields") %in%
        names(db_compound_vars)
    )
  )

  # add unidentifiable ions as unknowns and append them to compound list
  if (nrow(unidentifiable_ions) == 0) {
    return(list(
      entry_compound_dict = NULL,
      compound_attr = NULL,
      db_compound_vars = db_compound_vars
    ))
  }

  # create a new unique compoundName for compounds without a name or structure
  max_unknown_index <- db_compound_vars$max_unknown_index
  max_compoundId <- db_compound_vars$max_compoundId

  unidentifiable_ions <- unidentifiable_ions %>%
    dplyr::mutate(
      compoundName = paste0("unknown ", max_unknown_index + 1:dplyr::n()),
      compoundId = max_compoundId + 1:dplyr::n()
    )

  db_compound_vars$max_unknown_index <- max_unknown_index + nrow(unidentifiable_ions)
  db_compound_vars$max_compoundId <- max_compoundId + nrow(unidentifiable_ions)

  list(
    entry_compound_dict = unidentifiable_ions %>%
      dplyr::select(.entry_id, compoundId),
    compound_attr = unidentifiable_ions %>%
      dplyr::select(!!!rlang::syms(
        c(
          "compoundId",
          db_compound_vars$included_compound_fields
        )
      )),
    db_compound_vars = db_compound_vars
  )
}


#' Query Systematic Compounds
#'
#' @param systematicCompoundIds NULL to query all compounds or a numeric
#'   vector of systematic compound IDs
#' @param systematic_compounds_con A connection to an SQL database which
#'   stores the structure and systematic IDs of well-studied compounds.
#' @param db_update_version multiple versions of the database exist that can be
#'   differentiated by update date.
#' \describe{
#'   \item{all}{Queries all versions of the database. This is useful when
#'     previously-created systematicCompoundIds need to be matched}
#'   \item{recent}{Queries the most recent version of the database. This is
#'     useful when matching to compounds by name or structure.}
#'   \item{YYYY-MM-DD}{Provide a specific date with the format YYYY-MM-DD
#'     to query database entries from this date.}
#'   }
#'
#' @return a list containing:
#' \describe{
#'   \item{distinct_compounds}{unique database entries}
#'   \item{database_ids}{systematic IDs of compounds in KEGG, HMDB, ...}
#'   \item{compound_aliases}{alternative names for a compound}
#'   \item{metabolic_pathways}{KEGG pathways associated with a metabolite}
#'   }
#'
#' @export
query_systematic_compounds <- function(systematic_compounds_con,
                                       systematicCompoundIds = NULL,
                                       db_update_version = "recent") {
  stopifnot(all(
    class(systematicCompoundIds) %in% c("NULL", "numeric", "integer")
  ))

  checkmate::assertString(db_update_version)

  is_date <- stringr::str_detect(
    db_update_version,
    "^[0-9]{4}-[0-9]{2}-[0-9]{2}$"
  )

  if (!(db_update_version %in% c("recent", "all") || is_date)) {
    stop("db_update_version must be either \"recent\", \"all\", or a date")
  }

  db_versions <- dplyr::tbl(
    systematic_compounds_con,
    dbplyr::sql("SELECT update_date, COUNT(update_date) as N FROM distinct_compounds GROUP BY update_date")
  ) %>%
    dplyr::collect()

  if (is_date) {
    if (db_update_version %in% db_versions$update_date) {
      db_versions <- db_versions %>%
        dplyr::filter(update_date %in% db_update_version)
    } else {
      stop(
        "db_update_version is formatted as a date, but does not match any
         of the available update dates: ",
        paste(db_versions$update_date, collapse = ", ")
      )
    }
  } else if (db_update_version == "recent") {
    db_versions <- db_versions %>%
      dplyr::arrange(desc(update_date)) %>%
      dplyr::slice(1)
  }

  db_table_names <- DBI::dbListTables(systematic_compounds_con)
  db_tables <- purrr::map(
    db_table_names,
    function(x, systematic_compounds_con) {
      dplyr::tbl(systematic_compounds_con, x) %>%
        dplyr::collect()
    }, systematic_compounds_con
  )
  names(db_tables) <- db_table_names

  db_tables$distinct_compounds <- db_tables$distinct_compounds %>%
    # remove fields which should be derived by clamdb::canonical_structures()
    dplyr::select(-inchikey_connectivity, -inchikey_isotauto) %>%
    # filter to appropriate version(s) of DB
    dplyr::filter(update_date %in% db_versions$update_date)

  # filter other tables to respect versioned distinct compounds
  db_tables <- purrr::map(db_tables, function(x) {
    x %>%
      dplyr::semi_join(
        db_tables$distinct_compounds,
        by = "systematicCompoundId"
      )
  })

  # verify that distinct compounds have unique defining fields

  if (db_update_version != "all") {
    stopifnot(
      nrow(db_tables$distinct_compounds) ==
        length(unique(db_tables$distinct_compounds$systematicCompoundId))
    )
    stopifnot(
      nrow(db_tables$distinct_compounds) ==
        length(unique(
          toupper(stringr::str_trim(db_tables$distinct_compounds$name))
        ))
    )
    stopifnot(
      nrow(db_tables$distinct_compounds) ==
        length(unique(db_tables$distinct_compounds$smiles))
    )
    stopifnot(
      nrow(db_tables$distinct_compounds) ==
        length(unique(db_tables$distinct_compounds$inchikey))
    )
    stopifnot(
      nrow(db_tables$distinct_compounds) ==
        length(unique(db_tables$distinct_compounds$inchi))
    )
  }

  if (is.null(systematicCompoundIds)) {
    return(db_tables)
  }

  reduced_db_tables <- purrr::map(
    db_tables,
    function(x, sys_compounds_ids) {
      x %>%
        dplyr::filter(systematicCompoundId %in% sys_compounds_ids)
    },
    sys_compounds_ids = systematicCompoundIds
  )

  return(reduced_db_tables)
}

#' Standard DB Summarize Identifiable Compounds
#'
#' Match ions to systematic database, then past standards to determine which
#'   unique compounds are present. Then determine which to add, or overwrite
#'   past info.
#'
#' @param identifiable_ions filtered ion table from \code{ms_data_list}.
#' @inheritParams standard_db_summarize_compounds
#' @inheritParams standard_db_summarize_unidentifiable_compounds
standard_db_summarize_identifiable_compounds <- function(identifiable_ions,
                                                         mass_spec_standards_con,
                                                         systematic_compounds_con,
                                                         match_var_order,
                                                         db_compound_vars) {
  mass_spec_standards_past_compounds <- dplyr::tbl(
    mass_spec_standards_con,
    "compounds"
  ) %>%
    dplyr::collect()

  mass_spec_standards_annot_sources <- dplyr::tbl(
    mass_spec_standards_con,
    "compoundAnnotSource"
  ) %>%
    dplyr::collect()

  # load annotations from the most recent version of the systematic compounds database
  systematic_compounds <- query_systematic_compounds(
    systematic_compounds_con = systematic_compounds_con
  )

  max_compoundId <- db_compound_vars$max_compoundId

  # matches to systematic compounds database
  #   ions which match a structure, name, or alias in the systematic
  #   compounds database should be merged

  shared_systematic_match <- identifiable_ions %>%
    dplyr::rename(name = compoundName) %>%
    match_structures(
      target_compound_structures = systematic_compounds$distinct_compounds,
      target_compound_aliases = systematic_compounds$compound_aliases,
      match_var_order = match_var_order,
      drop_ambiguous = TRUE
    ) %>%
    dplyr::select(.entry_id, systematicCompoundId, match_var) %>%
    dplyr::filter(!is.na(systematicCompoundId)) %>%
    dplyr::mutate(match_tag = paste0("s", systematicCompoundId))

  # matches to current database

  if (nrow(mass_spec_standards_past_compounds) != 0) {
    shared_standard_match <- identifiable_ions %>%
      dplyr::rename(name = compoundName) %>%
      match_structures(
        target_compound_structures = mass_spec_standards_past_compounds %>%
          dplyr::select(
            systematicCompoundId = compoundId,
            name = compoundName,
            inchikey
          ) %>%
          dplyr::mutate(name = iconv(name, "UTF-8", "UTF-8", sub = "xx")),
        match_var_order = setdiff(match_var_order, "alias"),
        drop_ambiguous = TRUE
      ) %>%
      # filter to matches
      dplyr::filter(!is.na(systematicCompoundId)) %>%
      dplyr::select(.entry_id, compoundId = systematicCompoundId, match_var) %>%
      dplyr::mutate(match_tag = paste0("d", compoundId))

    distinct_annot_sources <- mass_spec_standards_annot_sources %>%
      dplyr::distinct(compoundId, inchikey, compoundName) %>%
      # remove non-UTF characters if they exist
      dplyr::mutate(
        compoundName = iconv(compoundName, "UTF-8", "UTF-8", sub = "xx"),
        compoundName = stringr::str_replace(
          toupper(stringr::str_trim(compoundName)),
          "\\-",
          " "
        )
      )

    # aggregating over all past names and structures
    #   in addition to the current consensus meta-data
    shared_standard_annot_match_structure <- identifiable_ions %>%
      dplyr::select(.entry_id, inchikey) %>%
      dplyr::inner_join(
        distinct_annot_sources %>%
          dplyr::filter(!is.na(inchikey)) %>%
          dplyr::distinct(compoundId, inchikey) %>%
          dplyr::mutate(match_tag = paste0("d", compoundId)),
        by = "inchikey"
      )

    shared_standard_annot_match_name <- identifiable_ions %>%
      dplyr::select(.entry_id, compoundName) %>%
      dplyr::mutate(compoundName = stringr::str_replace(
        toupper(stringr::str_trim(compoundName)),
        "\\-",
        " "
      )) %>%
      dplyr::inner_join(
        distinct_annot_sources %>%
          dplyr::filter(!is.na(compoundName)) %>%
          dplyr::distinct(compoundId, compoundName) %>%
          dplyr::mutate(match_tag = paste0("d", compoundId)),
        by = "compoundName"
      )

    # enforce that compounds which match the same systematic compound
    #   are also aggregated. systematic compounds have more names/aliases
    #   so some old and new compounds may match a systematic compound without
    #   matching one another.

    shared_standard_systematic_id <- shared_systematic_match %>%
      dplyr::inner_join(
        mass_spec_standards_past_compounds %>%
          dplyr::select(compoundId, systematicCompoundId),
        by = "systematicCompoundId"
      ) %>%
      dplyr::distinct(.entry_id, compoundId, systematicCompoundId) %>%
      dplyr::mutate(match_tag = paste0("d", compoundId))

    standard_matches <- do.call(
      dplyr::bind_rows,
      list(
        shared_standard_match,
        shared_standard_annot_match_structure,
        shared_standard_annot_match_name,
        shared_standard_systematic_id
      )
    )
  } else {
    standard_matches <- tibble::tibble(
      .entry_id = seq_along(NULL),
      match_tag = NA_character_[-1]
    )
  }

  # internal consistency
  # enforce that entries with matching inchikey or compoundName
  #   are aggregated

  shared_structure <- identifiable_ions %>%
    dplyr::select(.entry_id, inchikey) %>%
    dplyr::inner_join(
      identifiable_ions %>%
        dplyr::filter(!is.na(inchikey)) %>%
        dplyr::distinct(inchikey) %>%
        dplyr::mutate(match_tag = safe_tag_paste(., "i")),
      by = "inchikey"
    )

  shared_name <- identifiable_ions %>%
    dplyr::select(.entry_id, compoundName) %>%
    dplyr::mutate(compoundName = stringr::str_replace(
      toupper(enc2utf8(compoundName)),
      "\\-",
      " "
    )) %>%
    dplyr::inner_join(
      identifiable_ions %>%
        dplyr::filter(!is.na(compoundName)) %>%
        dplyr::mutate(compoundName = stringr::str_replace(
          toupper(enc2utf8(compoundName)),
          "\\-",
          " "
        )) %>%
        dplyr::distinct(compoundName) %>%
        dplyr::mutate(match_tag = safe_tag_paste(., "n")),
      by = "compoundName"
    )

  # find all compound groups

  coannotation_graph <- list(
    shared_systematic_match,
    standard_matches,
    shared_structure,
    shared_name
  ) %>%
    purrr::map(select, .entry_id, match_tag) %>%
    dplyr::bind_rows() %>%
    igraph::graph_from_data_frame(directed = FALSE) %>%
    {
      tibble::tibble(
        node = igraph::V(.)$name,
        group = igraph::clusters(.)$membership
      )
    } %>%
    dplyr::mutate(node_type = dplyr::case_when(
      stringr::str_detect(node, "^[0-9]+$") ~ "entry",
      stringr::str_detect(node, "^s[0-9]+$") ~ "systematic_match",
      stringr::str_detect(node, "^d[0-9]+$") ~ "standard_match",
      TRUE ~ "linker"
    ))

  # bring in compound attributes from all grouping variables

  manual_annotations <- coannotation_graph %>%
    dplyr::filter(node_type == "entry") %>%
    dplyr::left_join(identifiable_ions %>%
      dplyr::select(!!!rlang::syms(c(
        ".entry_id",
        db_compound_vars$included_compound_fields
      ))) %>%
      dplyr::mutate(.entry_id = as.character(.entry_id)),
    by = c("node" = ".entry_id")
    )

  standard_annotations <- coannotation_graph %>%
    dplyr::filter(node_type == "standard_match") %>%
    dplyr::mutate(compoundId = as.integer(
      stringr::str_replace(node, "^d([0-9]+)$", "\\1")
    )) %>%
    dplyr::left_join(mass_spec_standards_past_compounds, by = "compoundId")

  systematic_annotations <- coannotation_graph %>%
    dplyr::filter(node_type == "systematic_match") %>%
    dplyr::left_join(
      shared_systematic_match %>%
        dplyr::distinct(match_tag, systematicCompoundId),
      by = c("node" = "match_tag")
    ) %>%
    dplyr::left_join(
      systematic_compounds$distinct_compounds,
      by = "systematicCompoundId"
    ) %>%
    dplyr::select(
      node,
      group,
      node_type,
      systematicCompoundId,
      inchikey,
      compoundName = name,
      smiles = smiles,
      mass = mw,
      formula = molform
    )

  coannotations <- do.call(dplyr::bind_rows, list(
    systematic_annotations,
    standard_annotations,
    manual_annotations
  )) %>%
    dplyr::mutate(node_type = ordered(
      node_type,
      levels = c("systematic_match", "standard_match", "entry")
    ))

  # correspondence between group and compoundId

  old_compounds <- coannotations %>%
    dplyr::filter(node_type %in% "standard_match") %>%
    dplyr::select(group, compoundId)

  # 1-1 unambiguous standard IDs
  standard_matches_unambiguous <- old_compounds %>%
    dplyr::group_by(group) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::ungroup()

  # 1-2+ ambiguous standard IDs
  # when a group is matched to multiple standards, convert  them into a single group

  standard_matches_ambiguous <- old_compounds %>%
    dplyr::group_by(group) %>%
    dplyr::filter(dplyr::n() >= 2) %>%
    dplyr::arrange(compoundId) %>%
    dplyr::mutate(consensus_compoundId = compoundId[1])

  standard_matches_ambiguous_consensus <- standard_matches_ambiguous %>%
    dplyr::distinct(group, consensus_compoundId) %>%
    dplyr::rename(compoundId = consensus_compoundId)

  standard_inconsistencies <- standard_matches_ambiguous %>%
    dplyr::arrange(group) %>%
    dplyr::left_join(mass_spec_standards_past_compounds, by = "compoundId")

  # compounds which will be converted from an old to a new compoundId

  compoundId_conversions <- standard_matches_ambiguous %>%
    dplyr::ungroup() %>%
    dplyr::filter(compoundId != consensus_compoundId) %>%
    dplyr::select(
      old_compoundId = compoundId,
      new_compoundId = consensus_compoundId
    )

  if (nrow(compoundId_conversions) != 0) {
    warning(
      "  ", nrow(compoundId_conversions),
      " existing compoundIds will be merged due to ambiguous structures/naming;
        see details in standard_inconsistencies"
    )
  }

  # 1-0 previously undefined standard IDs

  new_compounds <- coannotations %>%
    dplyr::anti_join(coannotations %>%
      dplyr::filter(node_type %in% "standard_match"),
    by = "group"
    ) %>%
    dplyr::distinct(group)

  if (nrow(new_compounds) != 0) {
    new_compounds <- new_compounds %>%
      # add a new compoundId
      dplyr::mutate(compoundId = max_compoundId + 1:dplyr::n())
  } else {
    new_compounds <- new_compounds %>%
      dplyr::mutate(compoundId = NA_integer_)
  }

  db_compound_vars$max_compoundId <- ifelse(
    nrow(new_compounds) != 0,
    max(new_compounds$compoundId),
    db_compound_vars$max_compoundId
  )

  # generate a mapping between each group and a single old/new compoundId
  group_compounds <- do.call(
    dplyr::bind_rows,
    list(
      standard_matches_unambiguous,
      standard_matches_ambiguous_consensus,
      new_compounds
    )
  )

  # compound attributes

  compound_attr <- coannotations %>%
    dplyr::group_by(group) %>%
    # prefer systematic names w/ early IDs first
    dplyr::arrange(node_type, systematicCompoundId) %>%
    dplyr::summarize_all(first_defined) %>%
    dplyr::select(-node, -node_type, -compoundId, -match_var) %>%
    dplyr::left_join(group_compounds, by = "group") %>%
    dplyr::select(-group)

  entry_compound_dict <- coannotations %>%
    dplyr::filter(node_type == "entry") %>%
    dplyr::select(.entry_id = node, group) %>%
    dplyr::mutate(.entry_id = as.numeric(.entry_id)) %>%
    dplyr::left_join(group_compounds, by = "group") %>%
    dplyr::select(.entry_id, compoundId)

  # quick checks on dimensionality

  stopifnot(nrow(entry_compound_dict) == nrow(identifiable_ions))
  stopifnot(all(!is.na(entry_compound_dict$compoundId)))
  stopifnot(
    length(entry_compound_dict$.entry_id) ==
      length(unique(entry_compound_dict$.entry_id))
  )

  # find instances when a single compound group is annotated with multiple
  #   systematic compound IDs suggesting abberant naming / structural
  #   definition

  ambiguous_systematic_ids <- coannotations %>%
    dplyr::filter(!is.na(systematicCompoundId)) %>%
    dplyr::distinct(group, systematicCompoundId) %>%
    dplyr::count(group) %>%
    dplyr::filter(n > 1)

  if (nrow(ambiguous_systematic_ids) != 0) {
    warning("  ", nrow(ambiguous_systematic_ids), " compounds match 2+ systematic compounds
  This is likely due to inaccurate specification of measurements' structures and/or name.
  See details in systematic_inconsistencies")
  }

  systematic_inconsistencies <- group_compounds %>%
    dplyr::inner_join(coannotations %>%
      dplyr::semi_join(ambiguous_systematic_ids, by = "group") %>%
      dplyr::select(-node, -compoundId) %>%
      dplyr::distinct(), by = "group") %>%
    dplyr::select(-group) %>%
    dplyr::arrange(compoundId, inchikey, compoundName, node_type)

  # return entry lookup table, attributes, and ID collisions errors

  list(
    entry_compound_dict = entry_compound_dict,
    compound_attr = compound_attr,
    compoundId_conversions = compoundId_conversions,
    standard_inconsistencies = standard_inconsistencies,
    systematic_inconsistencies = systematic_inconsistencies,
    db_compound_vars = db_compound_vars
  )
}

standard_db_summarize_ions <- function(ms_data_ions, compound_summaries, mass_spec_standards_con) {
  mass_spec_standards_past_ions <- dplyr::tbl(
    mass_spec_standards_con,
    dbplyr::sql("SELECT * FROM ions")
  ) %>%
    dplyr::collect()

  max_ionId <- tbl(
    mass_spec_standards_con,
    sql("SELECT MAX(ionId) FROM ions")
  ) %>%
    dplyr::collect() %>%
    unlist() %>%
    unname()
  if (is.na(max_ionId)) {
    max_ionId <- 0
  }

  ion_fields <- c("adductName", "mode", "precursorMz", "precursorCharge")
  included_ion_fields <- intersect(colnames(ms_data_ions), ion_fields)

  defining_ion_fields <- c("adductName", "mode", "precursorCharge")
  missing_defining_ion_fields <- setdiff(defining_ion_fields, colnames(ms_data_ions))
  if (length(missing_defining_ion_fields) != 0) {
    stop(
      paste(missing_defining_ion_fields, collapse = " & "),
      " are required variables in the ions table of ms_data_list:
        add these variables adding NAs where necessary"
    )
  }

  # match all ions to their compound and reduce to unique ions
  unique_ions <- ms_data_ions %>%
    # add compound to ions
    dplyr::left_join(compound_summaries$entry_compound_dict, by = ".entry_id") %>%
    # collapse over precursorMz for matched  (largely in case of rounding differences)
    dplyr::group_by(!!!quos(!!!rlang::syms(c("compoundId", defining_ion_fields)))) %>%
    dplyr::summarize_all(first_defined) %>%
    dplyr::ungroup() %>%
    dplyr::select(!!!quos(!!!rlang::syms(c("compoundId", included_ion_fields))))

  if (any(is.na(unique_ions$compoundId))) {
    stop("some ions could not be matched to their compound - this is a bug, please submit a GitHub issue...")
  }

  # search matched_ions against current ions table to determine whether any are already present
  # use this search to create an appended ions table and a dictionary between .entry_id and ionId

  past_ion_matches <- unique_ions %>%
    # join existing ions (which have an existing ionId)
    dplyr::left_join(mass_spec_standards_past_ions %>%
      dplyr::select(!!!quos(!!!rlang::syms(c("ionId", "compoundId", defining_ion_fields)))),
    by = c("compoundId", defining_ion_fields)
    )

  # combine existing ionId with new annotations
  identified_ion_matches <- past_ion_matches %>%
    # past ions
    dplyr::filter(!is.na(ionId))

  new_ion_matches <- past_ion_matches %>%
    # no current ionId
    dplyr::filter(is.na(ionId))

  if (nrow(new_ion_matches) > 0) {
    # give new ions a new id
    identified_ion_matches <- identified_ion_matches %>%
      dplyr::bind_rows(new_ion_matches %>%
        dplyr::mutate(ionId = max_ionId + 1:dplyr::n()))
  }

  entry_ion_dict <- ms_data_ions %>%
    dplyr::left_join(compound_summaries$entry_compound_dict, by = ".entry_id") %>%
    dplyr::left_join(identified_ion_matches %>%
      dplyr::select(!!!quos(!!!rlang::syms(c("ionId", "compoundId", defining_ion_fields)))),
    by = c("compoundId", defining_ion_fields)
    ) %>%
    dplyr::select(.entry_id, ionId)

  missing_dict_entries <- length(union(entry_ion_dict$.entry_id, ms_data_ions$.entry_id)) - length(intersect(entry_ion_dict$.entry_id, ms_data_ions$.entry_id))
  if (missing_dict_entries != 0) {
    stop(missing_dict_entries, " .entry_ids could not be matched when matching .entry_ids to ionIds - this is a bug, please submit a GitHub issue...")
  }
  if (nrow(entry_ion_dict) != nrow(ms_data_ions)) {
    stop("some ions were matched to multiple compounds")
  }
  if (sum(is.na(entry_ion_dict$ionId)) != 0) {
    stop(sum(is.na(entry_ion_dict$ionId)), " compound .entry_ids could not be matched to compounds")
  }

  appended_ions <- identified_ion_matches %>%
    dplyr::distinct(!!!quos(!!!rlang::syms(c("ionId", "compoundId", included_ion_fields)))) %>%
    # do not append already present ions
    dplyr::filter(ionId > max_ionId)

  list(
    entry_ion_dict = entry_ion_dict,
    appended_ions = appended_ions
  )
}

#' Convert Compound
#'
#' Convert an old compoundId to a new value throughout the database and update ions accordingly
#'
#' @param compoundId_conversions tibble containing \code{old_compoundIds} and \code{new_compoundId}
#' @inheritParams standard_db_populate_data
#'
#' @return update compounds and associated ions in \code{mass_spec_standards_con}
standard_db_convert_compounds <- function(compoundId_conversions, mass_spec_standards_con) {
  if (nrow(compoundId_conversions) == 0) {
    return(invisible(0))
  }

  # update compounds

  # update compounds - compoundAnnotSource
  select_old_annot_query <- sprintf(
    "SELECT * FROM compoundAnnotSource WHERE compoundId IN (%s)",
    paste(as.character(compoundId_conversions$old_compoundId), collapse = ", ")
  )
  old_annot_sources <- dplyr::tbl(mass_spec_standards_con, dbplyr::sql(select_old_annot_query)) %>% dplyr::collect()

  annot_sources_update <- old_annot_sources %>%
    dplyr::left_join(compoundId_conversions, by = c("compoundId" = "old_compoundId")) %>%
    dplyr::select(-compoundId) %>%
    dplyr::rename(compoundId = new_compoundId)

  # remove old compound entriess
  delete_old_annot_statement <- sprintf(
    "DELETE FROM compoundAnnotSource WHERE compoundId IN (%s)",
    paste(as.character(compoundId_conversions$old_compoundId), collapse = ", ")
  )
  DBI::dbSendStatement(mass_spec_standards_con, delete_old_annot_statement)

  # add new versions
  insert_new_annot_statement <- sprintf(
    "INSERT FROM compoundAnnotSource WHERE compoundId IN (%s)",
    paste(as.character(compoundId_conversions$old_compoundId), collapse = ", ")
  )
  DBI::dbWriteTable(mass_spec_standards_con, "compoundAnnotSource", annot_sources_update, append = TRUE, row.names = FALSE)

  # update ions

  # update match ions to old or new existing ionIds

  mass_spec_standards_past_ions <- dplyr::tbl(mass_spec_standards_con, dbplyr::sql("SELECT * FROM ions")) %>%
    dplyr::collect()
  max_ionId <- max(mass_spec_standards_past_ions$ionId)

  defining_ion_fields <- c("adductName", "mode", "precursorCharge")

  old_ions <- compoundId_conversions %>%
    dplyr::inner_join(mass_spec_standards_past_ions, by = c("old_compoundId" = "compoundId")) %>%
    dplyr::rename(old_ionId = ionId, old_precursorMz = precursorMz)

  new_ions <- compoundId_conversions %>%
    dplyr::distinct(new_compoundId) %>%
    dplyr::inner_join(mass_spec_standards_past_ions, by = c("new_compoundId" = "compoundId")) %>%
    dplyr::rename(new_ionId = ionId, new_precursorMz = precursorMz)

  old_new_ion_match <- old_ions %>%
    dplyr::left_join(new_ions, by = c("new_compoundId", defining_ion_fields))

  newly_defined_ions <- old_new_ion_match %>%
    dplyr::filter(is.na(new_ionId))

  # update ions (remove old and add new ions which weren't defined)
  if (nrow(newly_defined_ions) != 0) {
    # generate new ion entries for the compound if it is present for the old compound but not the new compound.
    newly_defined_ions <- newly_defined_ions %>%
      dplyr::select(-new_ionId) %>%
      dplyr::mutate(new_ionId = max_ionId + 1:dplyr::n())

    # add entry to ions table if new ions are defined
    appended_ions <- newly_defined_ions %>%
      dplyr::rename(compoundId = new_compoundId, ionId = new_ionId, precursorMz = old_precursorMz) %>%
      dplyr::select(-old_compoundId, -old_ionId, -new_precursorMz)

    # add new ions (which weren't defined for the new compound)
    DBI::dbWriteTable(mass_spec_standards_con, "ions", appended_ions, append = TRUE, row.names = FALSE)
  }

  # remove ions of the old compound
  delete_ions_statement <- sprintf(
    "DELETE FROM ions WHERE compoundId IN (%s)",
    paste(as.character(unique(old_ions$old_compoundId)), collapse = ", ")
  )
  DBI::dbSendStatement(mass_spec_standards_con, delete_ions_statement)

  ion_updates <- dplyr::bind_rows(
    newly_defined_ions,
    old_new_ion_match %>%
      dplyr::filter(!is.na(new_ionId))
  ) %>%
    dplyr::select(old_ionId, new_ionId)

  # overwrite ionId in fragmentation and elutions

  feature_tables <- c("fragmentations", "elutions")
  for (a_feature_table in feature_tables) {

    # select old ion's feature-level data
    past_feature_query <- sprintf(
      "SELECT * FROM %s WHERE ionId IN (%s)",
      a_feature_table,
      paste(as.character(unique(ion_updates$old_ionId)), collapse = ", ")
    )
    past_features <- dplyr::tbl(mass_spec_standards_con, dbplyr::sql(past_feature_query)) %>%
      dplyr::collect()

    if (nrow(past_features) == 0) {
      next
    }

    # remove old ion's feature-level data from db
    past_feature_remove_statement <- sprintf(
      "DELETE FROM %s WHERE ionId IN (%s)",
      a_feature_table,
      paste(as.character(unique(ion_updates$old_ionId)), collapse = ", ")
    )
    DBI::dbSendStatement(mass_spec_standards_con, past_feature_remove_statement)

    # add updasted ion's feature-level data to db

    updated_features <- past_features %>%
      dplyr::inner_join(ion_updates, by = c("ionId" = "old_ionId")) %>%
      dplyr::select(-ionId) %>%
      dplyr::rename(ionId = new_ionId)
    DBI::dbWriteTable(mass_spec_standards_con, a_feature_table, updated_features, append = TRUE, row.names = FALSE)
  }

  return(invisible(0))
}

#' Remove Old Compounds
#'
#' Compounds tables are updated as the standard database is populated. Whenever a compound is detected remove the
#' previous compound summary in case the entry is updated.
#'
#' @param compound_summaries output of \code{\link{standard_db_summarize_compounds}}
#' @inheritParams standard_db_populate_data
#'
#' @return delete compounds from mass_spec_standards_db
standard_db_remove_old_compounds <- function(compound_summaries, mass_spec_standards_con) {

  # compounds updated or removed when processing the current dataset

  removed_compoundIds <- union(compound_summaries$compound_attr$compoundId, compound_summaries$compoundId_conversions$old_compoundId)

  # sql doesn't handle mass deletions well so if many entries are being overwritten, break them into multiple queries
  delete_batch_size <- 1000

  tibble::tibble(batch = ceiling(seq(removed_compoundIds) / delete_batch_size), compoundId = removed_compoundIds) %>%
    plyr::dlply(.variables = "batch") %>%
    purrr::walk(function(x) {
      delete_statement <- sprintf(
        "DELETE from compounds WHERE compoundId IN (%s)",
        paste(as.character(x$compoundId), collapse = ", ")
      )

      DBI::dbSendStatement(mass_spec_standards_con, delete_statement)
    })
}


first_defined <- function(x) {
  first <- x[!is.na(x)]
  ifelse(length(x) > 0, x[1], NA)
}

safe_tag_paste <- function(., tag) {
  if (nrow(.) == 0) {
    NA_character_[-1]
  } else {
    paste0(tag, 1:nrow(.))
  }
}
