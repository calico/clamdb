#' Retrieve Standards Subset
#'
#' Retrieve a relevant subset of the \code{mass_spec_standards_con} database based on the class of chemicals and mode.
#'
#' If possible, a pre-computed subset (saved .RDS file) will be retrieved, instead of querying the database directly.
#' If no applicable .RDS file can be retrieved, the database will be queried directly.
#'
#' @inheritParams standard_db_init
#' @param matched_method string representation encompassing chemical_class, ms_mode, and chromatographic_method information
#' @param r_scripts_path system path containing R scripts, notably mzkit.R
#' @param chemical_class Type of chemical of interest: e.g., metabolites
#' @param ms_mode Mode mass spec is running in: Negative and Positive
#' @param chromatographic_method Chromatography method used
#'
#' @return a list containing relevant entries of the \code{mass_spec_standards} database.
#'
#' @export
get_standards_subset <- function(mass_spec_standards_con, matched_method, r_scripts_path, chemical_class, ms_mode, chromatographic_method) {

  # path to an existing/to-be-created standards .Rds
  saved_standards_path <- create_saved_standards_path(r_scripts_path, matched_method)

  # Use saved file instead of connection to DB for testing
  if (class(mass_spec_standards_con) == "character" && mass_spec_standards_con == "use-rds-file") {
    if (file.exists(saved_standards_path)) {
      standards_data <- readRDS(saved_standards_path)
      return(standards_data)
    } else {
      stop(paste0("Cached .rds database file is missing! looked for file: ", saved_standards_path))
    }
  }

  standards_data <- read_standards_subset(mass_spec_standards_con, matched_method, r_scripts_path)

  if (is.null(standards_data)) {
    standards_data <- query_standards_subset(
      mass_spec_standards_con,
      chemical_class,
      ms_mode,
      chromatographic_method
    )

    tryCatch(
      {
        saveRDS(standards_data, file = saved_standards_path)
        debugr::dwatch(msg = paste0("Saved standards_data as .rds file:\n", saved_standards_path))
      },
      error = function(err) {
        debugr::dwatch(msg = paste0("Unable to save standards_data as .rds file:\n", saved_standards_path))
      }
    )
  }

  standards_data
}

create_saved_standards_path <- function(r_scripts_path, matched_method) {
  do.call(
    file.path,
    as.list((c(
      stringr::str_split(
        r_scripts_path,
        pattern = .Platform$file.sep
      )[[1]] %>%
        {
          .[-length(.)]
        }, "assets", "standards_data",
      stringr::str_replace_all(matched_method, ".msp", ".rds")
    )))
  )
}

#' Read Saved Standards Subset
#'
#' Retrieve a relevant subset of the \code{mass_spec_standards_con} database based on the saved matched_method
#'
#' or NULL, if no corresponding file can be found.
#'
#' @inheritParams standard_db_init
#' @param matched_method string representation of configuration details, from standards DB
#' @param r_scripts_path system path containing R scripts, notably mzkit.R
#'
#' @export
read_standards_subset <- function(mass_spec_standards_con, matched_method, r_scripts_path) {
  saved_standards_data <- do.call(file.path, as.list((c(stringr::str_split(r_scripts_path, pattern = .Platform$file.sep)[[1]] %>%
    {
      .[-length(.)]
    }, "assets", "standards_data", stringr::str_replace_all(matched_method, ".msp", ".rds")))))

  if (file.exists(saved_standards_data)) {
    before <- Sys.time()
    standards_data <- readRDS(saved_standards_data)
    after <- Sys.time()

    debugr::dwatch(msg = paste("Retrieved saved standards_data object in", difftime(after, before, units = "secs"), "seconds"))
    # Retrieved saved standards_data object in 29.8922140598297 seconds

    current_db_create_time_string <- get_db_create_time(mass_spec_standards_con)

    if (is.null(current_db_create_time_string)) {
      debugr::dwatch(msg = paste("current_db_create_time_string found to be NULL.\nWill re-extract standards_data from standards database."))
    } else if (!is.null(standards_data$createTime) && standards_data$createTime == current_db_create_time_string) {
      debugr::dwatch(msg = paste("saved_standards_subset verified to correspond to most recent version."))
      return(standards_data)
    } else {
      debugr::dwatch(msg = paste("saved_standards_subset does not agree with most recent version of standards database.\nWill re-extract standards_data from standards database."))
    }
  } else {
    debugr::dwatch(msg = paste("Unable to find saved_standards_subset file\n", saved_standards_data, "\nWill re-extract standards_data from standards database."))
  }

  NULL
}

#' Retrieve update string
#'
#' Retrieve date and time of last known update to database.
#'
#' Used as a version, to determine if standards need to be re-extracted, or can safely be re-parsed from saved .rds file.
#'
#' @param mass_spec_standards_con connection to cloudSQL database
#'
#' @export
get_db_create_time <- function(mass_spec_standards_con) {
  db_info <- RMySQL::dbGetInfo(mass_spec_standards_con)

  db_create_time_string <- mass_spec_standards_con %>%
    dplyr::tbl(dbplyr::sql(sprintf("SELECT CREATE_TIME
                                  FROM information_schema.tables
                                  WHERE TABLE_SCHEMA ='%s' AND TABLE_NAME = 'fragmentationData'", db_info$dbname))) %>%
    dplyr::collect()

  db_create_time_string$CREATE_TIME
}

#' Query Standards Subset
#'
#' Extract a relevant subset of the \code{mass_spec_standards_con} database based on class of chemicals and mode.
#'
#' @inheritParams standard_db_init
#' @param chemical_class Type of chemical of interest: e.g., metabolites
#' @param ms_mode Mode mass spec is running in: Negative and Positive
#' @param chromatographic_method Chromatography method used
#' @param sample_filter_quosure An optional quosure that can be used to filter
#'   samples before extracting elutions and spectra
#'
#' @return a list containing relevant entries of the \code{mass_spec_standards} database.
#'
#' @export
query_standards_subset <- function(mass_spec_standards_con,
                                   chemical_class,
                                   ms_mode,
                                   chromatographic_method,
                                   sample_filter_quosure = NULL) {
  checkmate::assertCharacter(chemical_class, len = 1)
  checkmate::assertCharacter(ms_mode, len = 1)

  formatted_chromatographic_method <- parse_chromatographic_method(
    chromatographic_method
  )

  sample_summary <- dplyr::tbl(
    mass_spec_standards_con,
    dbplyr::sql("SELECT * FROM samples")
  ) %>%
    dplyr::collect()

  if (!is.null(sample_filter_quosure)) {
    # reduce samples based on a provided quosure
    checkmate::assertClass(sample_filter_quosure, "quosure")
    sample_summary <- sample_summary %>%
      dplyr::filter(!!sample_filter_quosure)
  }

  if (nrow(sample_summary %>% dplyr::filter(chemicalClass == chemical_class)) == 0) {
    stop("zero entries for chemical_class = ", chemical_class, " in standard database")
  }
  if (nrow(sample_summary %>% dplyr::filter(mode == ms_mode)) == 0) {
    stop("zero entries for mode = ", ms_mode, " in standard database")
  }

  chromatographic_method_number_matches <- sample_summary %>%
    dplyr::filter(!is.na(chromatographicMethod)) %>%
    dplyr::mutate(formatted_chromatographic_method = purrr::map(chromatographicMethod, parse_chromatographic_method)) %>%
    tidyr::unnest_wider(formatted_chromatographic_method) %>%
    dplyr::filter(chromatographic_method_number == formatted_chromatographic_method$chromatographic_method_number)

  if (nrow(chromatographic_method_number_matches) == 0) {
    warning("zero entries for chromatographic_method = ", chromatographic_method, " in standard database.\nMethods which just differ in version number (e.g., A vs. B) were also not present.")
  }

  # #################################################################### #
  before <- Sys.time()

  selected_samples <- mass_spec_standards_con %>%
    dplyr::tbl(dbplyr::sql(sprintf("SELECT *
                FROM samples
                WHERE chemicalClass = '%s'", chemical_class))) %>%
    dplyr::collect() %>%
    dplyr::mutate(date = lubridate::ymd(date))

  if (!is.null(sample_filter_quosure)) {
    # reduce samples based on a provided quosure
    checkmate::assertClass(sample_filter_quosure, "quosure")
    selected_samples <- selected_samples %>%
      dplyr::filter(!!sample_filter_quosure)
  }

  # separate sample chromatographic methods into a number and version

  n_samples_with_chromatography <- selected_samples %>%
    dplyr::filter(!is.na(chromatographicMethod)) %>%
    nrow()

  if (n_samples_with_chromatography == 0) {
    selected_samples <- selected_samples %>%
      dplyr::mutate(
        chromatographic_method_number = NA_character_,
        chromatographic_method_version = NA_character_
      )
  } else {
    selected_samples <- selected_samples %>%
      dplyr::left_join(
        selected_samples %>%
          dplyr::filter(!is.na(chromatographicMethod)) %>%
          dplyr::mutate(formatted_chromatographic_method = purrr::map(chromatographicMethod, parse_chromatographic_method)) %>%
          tidyr::unnest_wider(formatted_chromatographic_method) %>%
          dplyr::select(sampleId, chromatographic_method_number, chromatographic_method_version),
        by = "sampleId"
      )
  }

  after <- Sys.time()
  debugr::dwatch(msg = paste("Determined selected samples from standards DB in", difftime(after, before, units = "secs"), "seconds"))
  # Determined selected samples from standards DB in 0.0607609748840332 seconds (2019-05-07)
  # #################################################################### #

  # #################################################################### #
  before <- Sys.time()

  selected_fragmentations <- mass_spec_standards_con %>%
    dplyr::tbl(dbplyr::sql(sprintf("SELECT fragId, ionId, fragmentations.sampleId, msLevel, energy
                    FROM fragmentations
                    INNER JOIN samples ON fragmentations.sampleId = samples.sampleId
                    WHERE samples.mode = '%s' AND samples.chemicalClass = '%s'", ms_mode, chemical_class))) %>%
    dplyr::collect() %>%
    dplyr::semi_join(selected_samples, by = "sampleId")

  # drop stacked collision energies

  invalid_fragmentations <- selected_fragmentations %>%
    dplyr::filter(!stringr::str_detect(energy, "^[0-9]+\\.?[0-9]*$"))

  if (nrow(invalid_fragmentations) != 0) {
    invalid_energies <- unique(invalid_fragmentations$energy)
    warning(glue::glue(
      "{nrow(invalid_fragmentations)} fragmentations were defined by
      energies that cannot be used; these fragmentations will be discarded.
      energies: {paste(invalid_energies, collapse = ', ')}"
    ))

    selected_fragmentations <- selected_fragmentations %>%
      dplyr::anti_join(invalid_fragmentations, by = "fragId")
  }

  after <- Sys.time()
  debugr::dwatch(msg = paste("Determined selected fragmentations (metadata only, not MS/MS peaks) from standards DB in", difftime(after, before, units = "secs"), "seconds"))
  # Determined selected fragmentations from standards DB in 8.55233407020569 seconds (2019-05-07)
  # #################################################################### #

  # #################################################################### #
  before <- Sys.time()

  selected_fragmentatData <- mass_spec_standards_con %>%
    dplyr::tbl(dbplyr::sql(sprintf("SELECT fd.fragId, mz, ic, label
                     FROM fragmentationData fd
                     INNER JOIN fragmentations f ON fd.fragId = f.fragId
                     INNER JOIN samples s ON f.sampleId = s.sampleId
                     WHERE s.mode = '%s' AND s.chemicalClass = '%s'", ms_mode, chemical_class))) %>%
    dplyr::collect() %>%
    dplyr::semi_join(selected_fragmentations, by = "fragId")


  after <- Sys.time()
  debugr::dwatch(msg = paste("Determined selected fragmentation data (MS/MS peaks) from standards DB in", difftime(after, before, units = "secs"), "seconds"))
  # Determined selected fragmentation data from standards DB in 199.042061090469 seconds (2019-05-07)
  # Determined selected fragmentation data (MS/MS peaks) from standards DB in 2928.65496206284 seconds [48.8 min] (2019-08-27)
  # #################################################################### #

  # TODO: This isn't offering any huge speedups - leaving alone for time being.

  # before <- Sys.time()
  #
  # all_frag_data <- mass_spec_standards_con %>%
  #   dplyr::tbl(dplyr::sql(sprintf("SELECT fd.fragId as fragId, fd.mz as mz, fd.ic as ic, fd.label as label, #fragmentationData
  #                                 f.ionId as ionId, f.sampleId as sampleId, f.msLevel as msLevel, f.energy as energy #fragmentations
  #                                 FROM fragmentationData fd
  #                                 INNER JOIN fragmentations f ON fd.fragId = f.fragId
  #                                 INNER JOIN samples s ON f.sampleId = s.sampleId
  #                                 WHERE s.mode = '%s' AND s.chemicalClass = '%s'", ms_mode, chemical_class))) %>%
  #   dplyr::collect()
  #
  # after <- Sys.time()
  # debugr::dwatch(msg = paste("Determined selected fragmentation data (MS/MS peaks) from standards DB in",difftime(after,before,units="secs"),"seconds"))

  # #################################################################### #
  before <- Sys.time()

  selected_elutions <- mass_spec_standards_con %>%
    dplyr::tbl(dbplyr::sql(sprintf("SELECT elutionId, ionId, elutions.sampleId, rtObserved, rtAligned
                    FROM elutions
                    INNER JOIN samples ON elutions.sampleId = samples.sampleId
                    WHERE samples.chemicalClass = '%s'", chemical_class))) %>%
    dplyr::collect() %>%
    # filter elutions to only those generated using the correct chromaographic method value (but they may differ in the version)
    dplyr::semi_join(selected_samples %>%
      dplyr::filter(chromatographic_method_number == formatted_chromatographic_method$chromatographic_method_number),
    by = "sampleId"
    )

  after <- Sys.time()
  debugr::dwatch(msg = paste("Determined selected elutions from standards DB in", difftime(after, before, units = "secs"), "seconds"))
  # Determined selected elutions from standards DB in 0.0599739551544189 seconds (2019-05-07)
  # #################################################################### #

  charge_suffix <- ifelse(ms_mode == "positive", "+", "-")

  all_ions <- dplyr::tbl(mass_spec_standards_con, dbplyr::sql(sprintf("SELECT ionId, compoundId, mode, adductName, precursorMz, precursorCharge FROM ions where ions.mode = '%s'", ms_mode))) %>%
    dplyr::collect() %>%
    dplyr::mutate(adductName = ifelse(endsWith(adductName, charge_suffix) == TRUE, adductName, paste0(adductName, charge_suffix)))

  # reduce samples and ions based on measurements
  selected_samples <- selected_samples %>%
    dplyr::filter(sampleId %in% unique(c(selected_fragmentations$sampleId, selected_elutions$sampleId)))

  selected_ions <- all_ions %>%
    dplyr::filter(ionId %in% unique(c(selected_fragmentations$ionId, selected_elutions$ionId)))

  # #################################################################### #
  before <- Sys.time()

  all_compounds <- dplyr::tbl(mass_spec_standards_con, "compounds") %>%
    dplyr::collect()

  after <- Sys.time()
  debugr::dwatch(msg = paste("Retrieved all compounds from standards DB in", difftime(after, before, units = "secs"), "seconds"))
  # Retrieved all compounds from standards DB in 6.07924199104309 seconds (2019-05-07)
  # #################################################################### #

  selected_compounds <- all_compounds %>%
    dplyr::semi_join(selected_ions, by = "compoundId")

  standards_data <- list(
    samples = selected_samples,
    ions = selected_ions,
    compounds = selected_compounds,
    elutions = selected_elutions,
    fragmentations = selected_fragmentations,
    fragmentationData = selected_fragmentatData,
    createTime = get_db_create_time(mass_spec_standards_con)
  )

  standards_data
}

parse_chromatographic_method <- function(chromatographic_method) {
  # turn a chromatographic method string into a list containing the method
  # number and the version letter

  checkmate::assertCharacter(chromatographic_method, len = 1)
  if (!stringr::str_detect(chromatographic_method, "MS\\-Chrom\\-[0-9]{3}\\-[A-Z]")) {
    stop("invalid chromatographic_method: ", chromatographic_method, "\nIt must match the following regular expression: MS\\-Chrom\\-[0-9]{3}\\-[A-Z]")
  }

  chromatographic_method_number <- stringr::str_extract(chromatographic_method, "MS\\-Chrom\\-[0-9]{3}")
  chromatographic_method_version <- stringr::str_match(chromatographic_method, "MS\\-Chrom\\-[0-9]{3}\\-([A-Z])")[2]

  return(list(
    chromatographic_method_number = chromatographic_method_number,
    chromatographic_method_version = chromatographic_method_version
  ))
}

#' Extract All Standard Ions
#'
#' @param standards_data output of \code{\link{query_standards_subset}}.
#' @param query_ionIds a vector of ionIds (or all ions if NULL).
#' @param chromatographic_method NULL or a length 1 character vector naming a chromatography method to use for retention time querying. Retention times are drawn from any ions and mode of a compound where the chromatographic method matches.
#' @param collision_energies NULL or a length 1, comma-separated set of collision energies.
#' @param is_lipid flag to describe if this is lipid data in the database, or metabolite data. If it is, the extraction is much simpler and faster.
#'
#' @return a list containing:
#' \itemize{
#'   \item{query_ions: ions & their retention time (for the specified chromatographic_method).}
#'   \item{query_fragmentations: fragmentation measurements that are relevant for the specified ions.}
#'   \item{query_fragmentationData: measurement data [mz, ic] for the query_fragmentations.}
#'   }
#'
#' @examples
#' \dontrun{
#' ms_mode <- "Negative" # this is available in the polarity field of each scan
#' chemical_class <- "metabolites"
#'
#' standards_data <- query_standards_subset(mass_spec_standards_con, chemical_class, ms_mode)
#' }
#'
#' @export
extract_standard_ions <- function(standards_data,
                                  query_ionIds = NULL,
                                  chromatographic_method = NULL,
                                  collision_energies = NULL,
                                  is_lipid = FALSE) {
  missing_required_tables <- setdiff(c("samples", "ions", "compounds", "elutions", "fragmentations", "fragmentationData", "createTime"), names(standards_data))
  if (length(missing_required_tables) != 0) {
    stop(paste(missing_required_tables, collapse = " & "), " are required tables which were not found in standards_data")
  }

  stopifnot(class(query_ionIds) %in% c("NULL", "integer", "numeric"))
  stopifnot(class(chromatographic_method) %in% c("NULL", "character"))
  if (class(chromatographic_method) == "character") {
    stopifnot(length(chromatographic_method) == 1)
    formatted_chromatographic_method <- parse_chromatographic_method(chromatographic_method)
  } else {
    chromatographic_method <- NA_character_
    formatted_chromatographic_method <- list(
      chromatographic_method_number = NA_character_,
      chromatographic_method_version = NA_character_
    )
  }

  stopifnot(class(collision_energies) %in% c("NULL", "character"))
  if (class(collision_energies) == "character") {
    stopifnot(length(collision_energies) == 1)
  } else {
    collision_energies <- NA_character_
  }

  # select ions of interest as a subset of the ions standards table
  query_ions <- if (is.null(query_ionIds)) {
    standards_data$ions
  } else {
    standards_data$ions %>%
      dplyr::filter(ionId %in% query_ionIds)
  }

  before <- Sys.time()

  if (!all(c("chromatographic_method_number", "chromatographic_method_version") %in% colnames(standards_data$samples))) {
    stop("chromatographic_method_number and/or chromatographic_method_version
          are missing from standards_data. This is likely because standards_data
          was read from an existing .Rds. Delete this file and re-try.")
  }

  # for RT, we don't care what type of ion it is; RT for any ion of the same compound (with a matched chromatography method will be valid)
  standard_retention_times <- query_ions %>%
    dplyr::left_join(standards_data$ions %>%
      dplyr::select(ionId_samplecmpd = ionId, compoundId), by = "compoundId") %>%
    dplyr::inner_join(standards_data$elutions, by = c("ionId_samplecmpd" = "ionId")) %>%
    dplyr::left_join(standards_data$samples, by = "sampleId") %>%
    # require a chromatographic method match
    dplyr::filter(chromatographic_method_number == formatted_chromatographic_method$chromatographic_method_number) %>%
    # calculate the version difference between each version # and the version number of chromatographic_method
    dplyr::mutate(chromatographic_method_diff = purrr::map_int(chromatographic_method_version,
      chromatographic_method_version_diff,
      version2 = formatted_chromatographic_method$chromatographic_method_version
    )) %>%
    dplyr::group_by(ionId) %>%
    # we may want to sort by sampleType here too, if there is some logical ordering of sampleTypes
    dplyr::arrange(chromatographic_method_diff, desc(date)) %>%
    # only take the most relevant RT (closest method version, then most recent)
    dplyr::slice(1) %>%
    dplyr::select(ionId, rtObserved, rtAligned)

  after <- Sys.time()
  debugr::dwatch(msg = paste("[clamdb<standard_db_query.R>::extract_standard_ions()]: Determined standard_retention_times in ", difftime(after, before, units = "secs"), "seconds"))

  # find the relevant subsets of energy for constructing fragment prediction bases
  # the measured energies which are closest to the experimental ionization energy

  fragmentation_types <- standards_data$fragmentations %>%
    dplyr::left_join(standards_data$samples %>%
      dplyr::select(sampleId, sampleType),
    by = "sampleId"
    )

  # Speedup for cases with only one MS/MS spectrum per compound
  if (is_lipid) {
    before <- Sys.time()

    standard_fragmentations <- query_ions %>%
      dplyr::mutate(collision_energy = collision_energies) %>%
      dplyr::inner_join(standards_data$fragmentations, by = "ionId") %>%
      dplyr::left_join(standards_data$samples, by = c("sampleId", "mode")) %>%
      dplyr::select(ionId, collision_energy, fragId, sampleType) %>%
      dplyr::mutate(energy_class = sampleType) %>%
      dplyr::select(-sampleType) %>%
      dplyr::mutate(n_energy_matches = 0) %>%
      dplyr::mutate(is_valid = TRUE) %>%
      dplyr::mutate(energy = NA)

    after <- Sys.time()

    debugr::dwatch(msg = paste("[clamdb<standard_db_query.R>::extract_standard_ions()]: Determined standard_fragmentations for lipids in", difftime(after, before, units = "secs"), "seconds"))
    # [clamdb<standard_db_query.R>::extract_standard_ions()]: Determined standard_fragmentations for lipids in 2.38062882423401 seconds
  } else {
    before <- Sys.time()

    fragmentation_data <- query_ions %>%
      dplyr::mutate(collision_energy = collision_energies) %>%
      dplyr::inner_join(standards_data$fragmentations, by = "ionId") %>%
      dplyr::left_join(standards_data$samples, by = c("sampleId", "mode"))

    nested_fragmentation_data <- fragmentation_data %>%
      tidyr::nest(fragmentations = c(-ionId, -collision_energy))

    after <- Sys.time()

    debugr::dwatch(msg = paste("[clamdb<standard_db_query.R>::extract_standard_ions()]: Determined nested_fragmentation_data in", difftime(after, before, units = "secs"), "seconds"))
    # [clamdb<standard_db_query.R>::extract_standard_ions()]: Determined nested_fragmentation_data in  26.720715045929 seconds

    if (nrow(nested_fragmentation_data) == 0) {
      stop("no relevant standard fragmentation data was found - matching cannot proceed")
    }

    before <- Sys.time()

    standard_fragmentations <- nested_fragmentation_data %>%
      # expand experimental energies into a vector (from comma delimited values)
      dplyr::mutate(collision_energy_vector = purrr::map(collision_energy, function(i) {
        as.numeric(trimws(strsplit(i, split = ",")[[1]]))
      })) %>%
      # for each ion select (1) fragmentations of pure compounds which are closest to individual experimental energies
      #                     (2) fragmentation of pure compounds across all energies
      #                     (3) fragmentations generated from samples with the same mixed collision energies
      #                     (4) theoretical spectra which are not tied to specific energies
      dplyr::mutate(selected_fragmentations = purrr::map2(collision_energy_vector, fragmentations, select_best_fragmentations)) %>%
      dplyr::select(ionId, collision_energy, selected_fragmentations) %>%
      tidyr::unnest(selected_fragmentations)

    after <- Sys.time()
    debugr::dwatch(msg = paste("[clamdb<standard_db_query.R>::extract_standard_ions()]: Determined standard_fragmentations for non-lipids in", difftime(after, before, units = "secs"), "seconds"))
    # On M0052-M005A-all: [clamdb<standard_db_query.R>::extract_standard_ions()]: Determined standard_fragmentations in  2370.41975188255 seconds
  }

  before <- Sys.time()

  ## summaries of relevant ions
  query_ions <- query_ions %>%
    dplyr::left_join(standard_retention_times, by = "ionId") %>%
    dplyr::left_join(standards_data$compounds, by = "compoundId")
  query_fragmentations <- standard_fragmentations
  query_fragmentationData <- standards_data$fragmentationData %>%
    dplyr::semi_join(standard_fragmentations, by = "fragId")

  after <- Sys.time()

  debugr::dwatch(msg = paste("[clamdb<standard_db_query.R>::extract_standard_ions()]: Determined query_ions in", difftime(after, before, units = "secs"), "seconds"))
  # [clamdb<standard_db_query.R>::extract_standard_ions()]: Determined query_ions in  2.00654697418213 seconds

  # TODO: update to tidyr::nest_legacy() -> tidyr::nest() and tidyr::unnest_legacy() -> tidyr::unnest() once tidyR fixes performance bug
  before <- Sys.time()

  # nest standards to one line per ion
  extracted_standard_ions <- query_ions %>%
    tidyr::nest_legacy(-ionId, -compoundId, -compoundName, -smiles, .key = "stdCompounds") %>%
    dplyr::left_join(
      query_fragmentations %>%
        dplyr::inner_join(query_fragmentationData %>%
          tidyr::nest_legacy(-fragId, .key = "fragData"),
        by = "fragId"
        ) %>%
        tidyr::nest_legacy(-ionId, .key = "stdFragmentations"),
      by = "ionId"
    )

  after <- Sys.time()

  debugr::dwatch(msg = paste("[clamdb<standard_db_query.R>::extract_standard_ions()]: Determined extracted_standard_ions in", difftime(after, before, units = "secs"), "seconds"))
  # [clamdb<standard_db_query.R>::extract_standard_ions()]: Determined extracted_standard_ions in 51.368586063385 seconds (M005A)

  extracted_standard_ions
}

chromatographic_method_version_diff <- function(version1, version2) {

  # takes two characters and computes their distance in the alphabet

  stopifnot(nchar(version1) == 1, version1 %in% LETTERS)
  stopifnot(nchar(version2) == 1, version2 %in% LETTERS)

  abs(which(LETTERS == version1) - which(LETTERS == version2))
}

#' Select Best Fragmentatiions
#'
#' @param target_energies a length one string of NCE fragmentation energies - if using stacked collision energies then separate energies with ","
#' @param fragmentations a set of fragmentations of different provenances with characteristic fragIds
#'
#' @return
#' fragmentations to be used in diverse categories:
#' \itemize{
#'   \item{matched_single_energy: match the closest energies to pure standards.}
#'   \item{all_single_energy: include all pure energies.}
#'   \item{mixed_energy: match all of the energies of a mixed energy samples [not implemented].}
#'   \item{insilico_w_intensity: match a theoretical spectra w/ intensities.}
#'   \item{insilico_no_intensity: match a theoretical spectra w/ only masses.}
#'   }
select_best_fragmentations <- function(target_energies, fragmentations) {
  fragmentation_energies <- fragmentations %>%
    dplyr::mutate(n_energies = stringr::str_count(energy, ",") + 1)

  fragmentation_summaries <- list()

  if (!is.na(target_energies) && nrow(fragmentation_energies %>% dplyr::filter(n_energies == 1)) != 0) {

    # analyze single energy spectra to match the target energies
    chosen_energy_matches <- fragmentation_energies %>%
      dplyr::distinct(sampleType, energy, n_energies) %>%
      dplyr::filter(n_energies == 1) %>%
      dplyr::distinct(energy) %>%
      dplyr::mutate(energy = as.numeric(energy)) %>%
      {
        expand.grid(energy = .$energy, target_energy = target_energies, stringsAsFactors = FALSE)
      } %>%
      dplyr::mutate(energy_abs_diff = abs(energy - target_energy))

    # find up to n standard fragmentations measurements with energies closest to n target energies (used to generate an experimental spectrum)
    # this is implemented as an iterative approach which sequentially takes the best overall standard energy
    # and then removes this option for other target energies

    continue <- TRUE
    matched_energies <- tibble::tibble(energy = 1, target_energy = 1) %>% dplyr::slice(-1)
    while (continue) {
      added_energy <- chosen_energy_matches %>%
        dplyr::anti_join(matched_energies, by = "target_energy") %>%
        dplyr::anti_join(matched_energies, by = "energy") %>%
        dplyr::arrange(energy_abs_diff) %>%
        dplyr::slice(1)

      if (nrow(added_energy) == 0) {
        continue <- FALSE
      } else {
        matched_energies <- dplyr::bind_rows(matched_energies, added_energy)
      }
    }

    single_energy_summary <- tibble::tibble(
      energy_class = "matched_single_energy",
      n_energy_matches = nrow(matched_energies),
      total_energy_abs_diff = sum(matched_energies$energy_abs_diff)
    ) %>%
      dplyr::mutate(is_valid = ifelse(n_energy_matches == length(target_energies) & total_energy_abs_diff <= 100, TRUE, FALSE))

    matched_single_energy_fragmentations <- fragmentation_energies %>%
      dplyr::semi_join(matched_energies %>% dplyr::mutate(energy = as.character(energy)), by = "energy") %>%
      dplyr::group_by(energy) %>%
      dplyr::arrange(sampleType, desc(date)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(fragId, energy) %>%
      cbind(single_energy_summary) %>%
      tibble::as_tibble()

    fragmentation_summaries$matched_single_energy <- matched_single_energy_fragmentations
  } else {
    fragmentation_summaries$matched_single_energy <- tibble::tibble(fragId = -1, energy = NA_character_, energy_class = "matched_single_energy", n_energy_matches = 0, is_valid = FALSE)
  }

  if (nrow(fragmentation_energies %>% dplyr::filter(n_energies == 1)) != 0) {
    fragmentation_summaries$all_single_energy <- fragmentation_energies %>%
      dplyr::filter(n_energies == 1) %>%
      dplyr::group_by(energy) %>%
      dplyr::arrange(sampleType, desc(date)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(fragId, energy) %>%
      dplyr::mutate(
        energy_class = "all_single_energy",
        n_energy_matches = dplyr::n(),
        is_valid = TRUE
      )
  } else {
    fragmentation_summaries$all_single_energy <- tibble::tibble(fragId = -1, energy = NA_character_, energy_class = "all_single_energy", n_energy_matches = 0, is_valid = FALSE)
  }

  # add methods for mixed energy
  fragmentation_summaries$mixed_energy <- tibble::tibble(fragId = -1, energy = NA_character_, energy_class = "mixed_energy", n_energy_matches = 0, is_valid = FALSE)

  # add methods for no energy (i.e., theoretical spectra)

  if (nrow(fragmentation_energies %>% dplyr::filter(sampleType == "insilico_w_intensity")) != 0) {
    fragmentation_summaries$insilico_w_intensity <- fragmentation_energies %>%
      dplyr::arrange(desc(date)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(fragId, energy) %>%
      dplyr::mutate(
        energy_class = "insilico_w_intensity",
        n_energy_matches = 0,
        is_valid = TRUE
      )
  } else {
    fragmentation_summaries$insilico_w_intensity <- tibble::tibble(fragId = -1, energy = NA_character_, energy_class = "insilico_w_intensity", n_energy_matches = 0, is_valid = FALSE)
  }

  if (nrow(fragmentation_energies %>% dplyr::filter(sampleType == "insilico_no_intensity")) != 0) {
    fragmentation_summaries$insilico_no_intensity <- fragmentation_energies %>%
      dplyr::arrange(desc(date)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(fragId, energy) %>%
      dplyr::mutate(
        energy_class = "insilico_no_intensity",
        n_energy_matches = 0,
        is_valid = TRUE
      )
  } else {
    fragmentation_summaries$insilico_no_intensity <- tibble::tibble(fragId = -1, energy = NA_character_, energy_class = "insilico_no_intensity", n_energy_matches = 0, is_valid = FALSE)
  }

  fragmentation_summaries %>%
    dplyr::bind_rows() %>%
    dplyr::filter(is_valid)
}

#' Match Method
#'
#' @param config_list a list of named parameters which includes "chromatographic_method", "mode", and "collision_energies" to match to an existing method
#' @inheritParams standard_db_init
#'
#' @return a matched method string used by \code{\link{get_standards_subset}}
#'
#' @export
match_method <- function(config_list, mass_spec_standards_con) {
  defined_methods <- dplyr::tbl(mass_spec_standards_con, "methods") %>%
    dplyr::collect()

  if (!config_list$chromatographic_method %in% defined_methods$chromatographicMethod) {
    stop(
      config_list$chromatographic_method,
      " is not a defined value for chromatographic_method.
          Existing methods are:\n",
      paste(unique(defined_methods$chromatographicMethod), collapse = "\n")
    )
  }

  matched_analytical_methods <- defined_methods %>%
    dplyr::filter(
      chromatographicMethod == config_list$chromatographic_method,
      collisionEnergies == config_list$collision_energies,
      mode == config_list$mode
    )

  # Try to rescue multiple matches by checking for exact match to methodId
  if (nrow(matched_analytical_methods) > 1) {
    exact_methodId_match <- matched_analytical_methods %>%
      dplyr::filter(methodId == config_list$methodId)

    if (nrow(exact_methodId_match) == 1) {
      matched_analytical_methods <- exact_methodId_match
    }
  }

  if (nrow(matched_analytical_methods) == 1) {
    matched_method <- paste0(gsub(" ", "_", matched_analytical_methods$analyticalMethod[1]), "-unnormalized-matched_single_energy.msp")
  } else if (nrow(matched_analytical_methods) == 0) {
    warning("no analytical methods were matched;")
    matched_method <- ""
  } else {
    warning(nrow(matched_analytical_methods), " analytical methods were matched")
    matched_method <- ""
  }

  return(matched_method)
}

#' Find Standards
#'
#' Search for standards by name or formula
#'
#' @param compounds a character vector of compound names or molecular formulas.
#' @param match_type type of match to use.
#' \describe{
#'   \item{name}{exact case-insensitive substring matching}
#'   \item{formula}{exact formula match}
#' }
#' @inheritParams standard_db_init
#'
#' @examples
#' \dontrun{
#' compound_names <- c("N-acetyl-L-aspart", "Acetyl CoA")
#' find_standards(compound_names, match_type = "name", mass_spec_standards_con)
#'
#' compound_formulas <- c("C6H9NO5", "C23H38N7O17P3S")
#' find_standards(compound_formulas, match_type = "formula", mass_spec_standards_con)
#' }
#'
#' @export
find_standards <- function(compounds, match_type = "name", mass_spec_standards_con) {
  checkmate::assertCharacter(compounds)
  checkmate::assertChoice(match_type, c("name", "formula"))

  all_compounds <- dplyr::tbl(mass_spec_standards_con, "compounds") %>%
    dplyr::collect()

  if (match_type == "name") {
    compound_matches <- tibble::tibble(query = compounds) %>%
      dplyr::mutate(matches = purrr::map(query, find_standards_filter, all_compounds)) %>%
      tidyr::unnest(matches) %>%
      dplyr::select(query, compoundId, compoundName, formula)
  } else {
    compound_matches <- tibble::tibble(query = compounds) %>%
      dplyr::inner_join(all_compounds, by = c("query" = "formula")) %>%
      dplyr::select(query, compoundId, compoundName)
  }

  ion_matches <- dplyr::tbl(mass_spec_standards_con, "ions") %>%
    dplyr::collect() %>%
    dplyr::semi_join(compound_matches, by = "compoundId") %>%
    tidyr::nest(ions = -compoundId)

  compound_matches <- compound_matches %>%
    dplyr::left_join(ion_matches, by = "compoundId")

  return(compound_matches)
}

find_standards_filter <- function(query, all_compounds) {

  # treat spaces and dashes equally

  query <- stringr::str_replace_all(query, pattern = "[ -]", replacement = "[ -]")

  all_compounds %>%
    dplyr::mutate(
      filter_name = iconv(compoundName, "UTF-8", "UTF-8", sub = ""),
      filter_name = toupper(filter_name)
    ) %>%
    dplyr::filter(stringr::str_detect(filter_name, toupper(query))) %>%
    dplyr::select(-filter_name)
}
