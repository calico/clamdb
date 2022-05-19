#' Attributes Liset to MS Data List
#'
#' @param ion_attributes a tbl containing .entry_id, compoundName, inchikey [optional], smiles [optional], adductName, mass, precursorMz, mode, precursorCharge, date [optional], path [optional], energy, msLevel, chromatographyMethod [optional], rtObserved [optional]
#' @param fragmentationData a tbl containing .entry_id, mz, ic, and optionally label
#' @param chemicalClass e.g., lipids, fatty acids, polar
#' @param sampleType e.g., "pure standard", "public", "observed", "insilico_no_intensity", "insilico_w_intensity"
#'
#' @return an \code{ms_data_list} used by \code{\link{standard_db_populate_data}}
#'
#' @export
dbimport_attributes_list_to_ms_data_list <- function(ion_attributes, fragmentationData, chemicalClass, sampleType) {
  stopifnot(class(chemicalClass) == "character", length(chemicalClass) == 1)
  stopifnot(class(sampleType) == "character", length(sampleType) == 1)

  ion_fields <- c(".entry_id", "compoundName", "inchikey", "smiles", "adductName", "mass", "precursorMz", "mode", "precursorCharge")
  included_ion_fields <- intersect(ion_fields, colnames(ion_attributes))

  ions <- ion_attributes %>%
    dplyr::select(!!!rlang::syms(included_ion_fields))

  sample_fields <- c(".entry_id", "mode", "date", "path", "chromatographicMethod")
  included_sample_fields <- intersect(sample_fields, colnames(ion_attributes))

  samples <- ion_attributes %>%
    dplyr::select(!!!rlang::syms(included_sample_fields)) %>%
    dplyr::mutate(
      chemicalClass = chemicalClass,
      sampleType = sampleType
    )

  fragmentation_fields <- c(".entry_id", "energy", "msLevel")
  included_fragmentation_fields <- intersect(fragmentation_fields, colnames(ion_attributes))

  fragmentations <- ion_attributes %>%
    dplyr::select(!!!rlang::syms(included_fragmentation_fields))

  fragmentationData_fields <- c(".entry_id", "mz", "ic", "label")
  included_fragmentation_fields <- intersect(fragmentationData_fields, colnames(fragmentationData))

  fragmentationData <- fragmentationData %>%
    dplyr::select(!!!rlang::syms(included_fragmentation_fields))

  output <- list(
    ions = ions,
    samples = samples,
    fragmentations = fragmentations,
    fragmentationData = fragmentationData
  )

  if (all(c("rtObserved", "chromatographicMethod") %in% colnames(ion_attributes))) {
    elutions <- ion_attributes %>%
      filter(
        !is.na(rtObserved),
        !is.na(chromatographicMethod)
      ) %>%
      dplyr::select(.entry_id, rtObserved)

    if (nrow(elutions) != 0) {
      output$elutions <- elutions
    }
  }

  standard_db_test_ms_data_list_classes(output)

  return(output)
}

#' Nest standards by key
#'
#' Join library_key to ms_data_list$standards by path and separate ms_data_list by libraray_id
#'
#' @inheritParams standard_db_populate_data
#' @param library_key a tibble containing path and library_id
#' @param chunk_size NA or numeric chunk size of entries to simultaneously database
#'
#' @return a nested tibble containing library_id and library_data (a list of ms_data_lists)
#'
#' @export
dbimport_nest_standard_by_key <- function(ms_data_list, library_key, chunk_size = NA_real_) {
  stopifnot(class(chunk_size) == "numeric", length(chunk_size) == 1)
  if (!is.na(chunk_size)) {
    stopifnot(chunk_size > 0)
  }

  library_key_entry_dict <- ms_data_list$samples %>%
    dplyr::left_join(library_key, by = "path") %>%
    dplyr::select(.entry_id, library_id)

  if (!is.na(chunk_size)) {
    library_key_entry_dict <- library_key_entry_dict %>%
      # split each library into chunks by chunk_size
      dplyr::group_by(library_id) %>%
      dplyr::mutate(chunk = ceiling(seq(0.0001, ceiling(dplyr::n() / chunk_size), length.out = dplyr::n()))) %>%
      dplyr::ungroup() %>%
      tidyr::unite(col = "library_id", library_id, chunk, sep = "__")
  }

  # nest all tibbles by library_id
  nested_standards <- ms_data_list$samples %>%
    dplyr::left_join(library_key_entry_dict, by = c(".entry_id")) %>%
    tidyr::nest(samples = -library_id) %>%
    dplyr::left_join(
      ms_data_list$ions %>%
        dplyr::left_join(library_key_entry_dict, by = c(".entry_id")) %>%
        tidyr::nest(ions = -library_id),
      by = "library_id"
    ) %>%
    dplyr::left_join(
      ms_data_list$fragmentations %>%
        dplyr::left_join(library_key_entry_dict, by = c(".entry_id")) %>%
        tidyr::nest(fragmentations = -library_id),
      by = "library_id"
    ) %>%
    dplyr::left_join(
      ms_data_list$fragmentationData %>%
        dplyr::left_join(library_key_entry_dict, by = c(".entry_id")) %>%
        tidyr::nest(fragmentationData = -library_id),
      by = "library_id"
    )

  # nest elutions upfront since it is optional

  if ("elutions" %in% names(ms_data_list)) {
    elution_nest <- ms_data_list$elutions %>%
      dplyr::left_join(library_key_entry_dict, by = c(".entry_id")) %>%
      tidyr::nest(elutions = -library_id)

    nested_standards <- nested_standards %>%
      dplyr::left_join(elution_nest, by = "library_id")
  }

  nest_lists <- list()
  for (i in 1:nrow(nested_standards)) {
    nest_list <- list(
      samples = nested_standards$samples[[i]],
      ions = nested_standards$ions[[i]],
      fragmentations = nested_standards$fragmentations[[i]],
      fragmentationData = nested_standards$fragmentationData[[i]]
    )

    if ("elutions" %in% colnames(nested_standards) && !is.null(nested_standards$elutions[[i]])) {
      nest_list$elutions <- nested_standards$elutions[[i]]
    }

    nest_lists[[i]] <- nest_list
  }

  nested_standards$library_list <- nest_lists

  return(nested_standards %>% dplyr::select(library_id, library_list))
}

#' Adduct to Charge
#'
#' @param adductName vector of adducts / precursor types
#'
#' @return a vector of integer charge values
#'
#' @export
dbimport_adduct_to_charge <- function(adductName) {
  adductName %>%
    stringr::str_extract("[0-9]?[+-]$") %>%
    {
      dplyr::case_when(
        . %in% c("-", "+") ~ 1L,
        TRUE ~ as.integer(stringr::str_extract(., "^[0-9]"))
      )
    }
}

#' Parse Collision Energies
#'
#' Determine whether energy is eV or NCE based on annotations and then summarize energies in NCE
#'
#' @inheritParams dbimport_attributes_list_to_ms_data_list
#' @param energy_name Variable name (quoted) specifying collision energies
#' @param default_energy_type Default energy type (options are eV, NCE or NA [default])
#'
#' @return a character vector of energies (to allow for stacked collision energies) formatted in NCEs
#'
#' @export
dbimport_parse_collision_energies <- function(ion_attributes, energy_name, default_energy_type = NA_character_) {
  stopifnot("precursorMz" %in% colnames(ion_attributes))
  stopifnot("precursorCharge" %in% colnames(ion_attributes))
  stopifnot(energy_name %in% colnames(ion_attributes))

  stopifnot(class(default_energy_type) == "character", length(default_energy_type) == 1, default_energy_type %in% c("eV", "NCE", NA_character_))

  energy_name_sym <- dplyr::sym(energy_name)

  ion_energies <- ion_attributes %>%
    dplyr::select(!!energy_name_sym, precursorMz, precursorCharge) %>%
    # add .row_number so we can reconstruct the order later on
    dplyr::mutate(
      .row_number = 1:dplyr::n(),
      # determine energy convention used
      energy_type = dplyr::case_when(
        stringr::str_detect(!!energy_name_sym, "(eV)|(V)|(ev)|(EV)") ~ "eV",
        stringr::str_detect(!!energy_name_sym, "NCE") ~ "NCE",
        stringr::str_detect(!!energy_name_sym, "\\%") ~ "NCE",
        stringr::str_detect(!!energy_name_sym, "\\(nominal\\)") ~ "NCE",
        stringr::str_detect(!!energy_name_sym, "normalized") ~ "NCE",
        TRUE ~ default_energy_type
      ),
      # determien numerical collision energy
      energy_value = dplyr::case_when(
        stringr::str_detect(!!energy_name_sym, "^[0-9.]+") ~ stringr::str_extract(!!energy_name_sym, "^[0-9.]+"),
        stringr::str_detect(!!energy_name_sym, "([0-9.]+[\\-\\/])+[0-9.]+") ~ stringr::str_extract(!!energy_name_sym, "([0-9.]+[\\-\\/])+[0-9.]+")
      )
    ) %>%
    # convert eV to NCE
    dplyr::select(mz = precursorMz, charge = precursorCharge, energy = energy_value, energy_type, .row_number) %>%
    dplyr::mutate(
      mz = as.numeric(mz),
      energy = as.character(energy)
    )

  # combine energies with different conventions
  transformed_energies <- do.call(dplyr::bind_rows, list(
    ion_energies %>%
      dplyr::filter(energy_type == "eV") %>%
      dplyr::mutate(energy_conv = as.character(dbimport_nce_conversion(., to_NCE = TRUE))),
    ion_energies %>%
      dplyr::filter(energy_type == "NCE") %>%
      dplyr::mutate(energy_conv = energy),
    ion_energies %>%
      dplyr::filter(is.na(energy_type) | !(energy_type %in% c("eV", "NCE"))) %>%
      dplyr::mutate(energy_conv = NA_character_)
  )) %>%
    dplyr::arrange(.row_number)

  if (!all(transformed_energies$.row_number == 1:nrow(transformed_energies))) {
    stop("some rows are missing in the transformed energy output, this is a bug")
  }

  return(transformed_energies$energy_conv)
}

#' NCE conversion
#'
#' @description Mapping between mass spec energy conventions - between normalized collision energies and electron volts.
#'
#' @param energy_table a tibble containing energy (in units of NCE or eV), mz and charge
#' @param to_NCE TRUE for eV -> NCE, FALSE for NCE -> eV
#'
#' @return a numeric vector corresponding to the rows in energy_table with converted energies
#'
#' @export
dbimport_nce_conversion <- function(energy_table, to_NCE = TRUE) {
  stopifnot(all(c("energy", "mz", "charge") %in% colnames(energy_table)))

  stopifnot(class(energy_table$energy) %in% c("integer", "numeric", "character"))
  stopifnot(class(energy_table$mz) %in% c("integer", "numeric"))
  stopifnot(class(energy_table$charge) == "integer")

  energy_table <- energy_table %>%
    dplyr::mutate(energy = as.numeric(energy))

  if (any(is.na(energy_table$energy))) {
    warning("some energies were converted to NA by dbimport_nce_conversion, likely because they were stacked collision energies which this function does not currently support")
  }

  charge_factors <- tibble::tribble(
    ~charge, ~charge_factor,
    1L, 1,
    2L, 0.9,
    3L, 0.85,
    4L, 0.8,
    5L, 0.75
  )

  missing_charges <- setdiff(energy_table$charge, charge_factors$charge)

  if (length(missing_charges) != 0) {
    stop("charge of ", paste(missing_charges, collapse = ", "), " missing from charge lookup table")
  }

  energy_conversions <- energy_table %>%
    dplyr::left_join(charge_factors, by = "charge") %>%
    dplyr::mutate(conversion_factor = mz / (500 * charge_factor))

  # Absolute energy (eV) = (settling NCE) x (Isolation center) / (500 m/z) x (charge factor)

  if (to_NCE) {
    # eV -> NCE
    return(energy_conversions$energy / energy_conversions$conversion_factor)
  } else {
    # NCE -> eV
    return(energy_conversions$energy * energy_conversions$conversion_factor)
  }
}
