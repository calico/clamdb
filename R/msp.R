#' Download .msp Library
#'
#' Download a public msp fragmentation library.
#'
#' @param savepath Output path
#' @param library_url Url of the library to download
#' @param unzip_dl True if the file is a .zip file
#' @param overwrite Overwrite existing results
#' @param method_dl method to download libraries from the internet (default: "libcurl")
#'
#' @returns 0 invisibly
#'
#' @export
download_msp_library <- function(savepath,
                                 library_url,
                                 unzip_dl,
                                 method_dl = "libcurl",
                                 overwrite = FALSE) {
  if (file.exists(savepath) & overwrite == FALSE) {
    warning(
      "library found at ",
      savepath,
      " and overwrite is FALSE, so library will not be modified\n"
    )
    return(invisible(0))
  }

  if (unzip_dl) {
    zipped_path <- file.path("/tmp", "msp_zip.zip")

    download.file(library_url, zipped_path, method = method_dl)
    zipped_files <- unzip(zipped_path, list = TRUE)
    if (nrow(zipped_files) > 1) {
      stop(
        "multiple files are zipped into ",
        library_url,
        " the enclosed library files cannot be automatically downloaded"
      )
    }

    unzip(zipped_path, exdir = "/tmp")
    file.copy(from = file.path("/tmp", zipped_files$Name), to = savepath)

    file.remove(zipped_path)
    file.remove(file.path("/tmp", zipped_files$Name))
  } else {
    download.file(library_url, savepath, method = method_dl)
  }

  invisible(0)
}

#' Parse All Msps
#'
#' @param msp_libraries tibble containing
#' \describe{
#'   \item{library}{A directory}
#'   \item{subset}{An individual library}
#'   \item{filename}{Filename of msps to parse}
#'   \item{savepath}{Full path of \code{filename}}
#' }
#' @inheritParams run_msp_parse
#' @param big_co Max .msp file size (in Gb) to automatically process
#'
#' @return 0 invisibly
#'
#' @export
parse_all_msps <- function(msp_libraries,
                           conda_env,
                           conda_path,
                           path_clamdb_repository,
                           big_co = 1) {
  libraries_to_parse <- msp_libraries %>%
    dplyr::mutate(
      file_noext = stringr::str_replace(savepath, "\\.([mM][sSgG][pPfF])$", ""),
      completed = dplyr::case_when(
        !file.exists(paste0(file_noext, "_attr.feather")) ~ FALSE,
        !file.exists(paste0(file_noext, "_frag.feather")) ~ FALSE,
        TRUE ~ TRUE
      ),
      file_size = file.size(savepath) / (2^30)
    ) %>%
    dplyr::filter(!completed)

  # run some libraries

  big_libraries_to_parse <- libraries_to_parse %>%
    # >1Gb
    dplyr::filter(file_size >= big_co)
  if (nrow(big_libraries_to_parse) > 0) {
    message(glue::glue(
      "{nrow(big_libraries_to_parse)} libraries are larger than {big_co}Gb
        and should be processed by directly calling mspParse.py:
        \n{paste(big_libraries_to_parse$savepath, collapse = '\n')}"
    ))
  }

  libraries_to_parse %>%
    dplyr::filter(file_size < big_co) %>%
    {
      purrr::walk(.$savepath, function(x) {
        file_noext <- stringr::str_replace(x, "\\.msp$", "")

        if (!all(file.exists(
          paste(file_noext, c("_attr.feather", "_frag.feather"), sep = "")
        ))) {
          print(glue::glue(
            "converting {x} to feather summaries of spectra attributes and fragmentations"
          ))
          run_msp_obj <- try(
            clamdb::run_msp_parse(x, conda_env, conda_path, path_clamdb_repository),
            silent = TRUE
          )
          if (class(run_msp_obj) == "try-error") {
            warning(
              x,
              " could not be processed
              run mspParse.py directly to debug your file\n"
            )
          }
        }
      })
    }

  invisible(0)
}

#' Run MSP parse
#'
#' Parse msp file into feather tables using python.
#'
#' @param msp_path path to msp fragmentation file
#' @param path_clamdb_repository path to this repository.
#'    Necessary python scripts are found in subdirectory \code{inst/python},
#'    and will need to be installed and configured for msp parsing to work correctly.
#' 
#' @inheritParams canonicalize_structures
#'
#' @return attributes and fragmentations .feather tables in the
#'   \code{msp_path} directory
#'
#' @export
run_msp_parse <- function(msp_path,
                          conda_env,
                          conda_path,
                          path_clamdb_repository) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop('The "reticulate" package must be installed to use this function',
      call. = FALSE
    )
  }

  reticulate::use_condaenv(conda_env, conda = conda_path)
  required_modules <- c("pandas", "feather")
  test_reticulate_py_dependences(required_modules)

  mspParse_fxns <- reticulate::import_from_path(
    module = "mspParse",
    path = file.path(path_clamdb_repository, "inst", "python"),
    convert = TRUE
  )

  mspParse_fxns$format_msp_and_save(msp_path, "UTF-8")

  invisible(0)
}

#' Read .msp Feathers
#'
#' Read attributes and fragmentation tables produced by
#'   \code{\link{run_msp_parse}}.
#'
#' @param savepath Path to directory containing msp files
#' @param file_id String name of previously parsed file
#' @param feather_type which type of feather to read: both, attr, or frag
#'
#' @return A list containing compound_attributes and
#'   compound_fragmentationData tibbles.
#'
#' @export
read_msp_feathers <- function(savepath, file_id, feather_type = "all") {
  stopifnot(
    "character" %in% class(feather_type), length(feather_type) == 1,
    feather_type %in% c("all", "attr", "frag")
  )

  path_noext <- stringr::str_replace(savepath, ".[a-zA-Z]+$", "")

  if (
    file.exists(paste0(path_noext, "_attr.feather")) &&
      file.exists(paste0(path_noext, "_frag.feather"))
  ) {
    output <- list()

    if (feather_type %in% c("all", "attr")) {
      output[["compound_attributes"]] <- feather::read_feather(
        paste0(path_noext, "_attr.feather")
      ) %>%
        dplyr::mutate(file_id = file_id)
    }

    if (feather_type %in% c("all", "frag")) {
      output[["compound_fragmentationData"]] <- feather::read_feather(
        paste0(path_noext, "_frag.feather")
      ) %>%
        dplyr::mutate(file_id = file_id)
    }

    return(output)
  } else {
    warning(path_noext, " feather outputs not found\n")
    return(NULL)
  }
}


#' Format MSP Libraries
#'
#' Reformat the output \code{read_msp_feathers} so that attributes use a common
#' terminology and flag entries which are missing key attributes.
#'
#' @param msp_compound_attributes a list of fragmentation attributes and
#'   data produced by \code{\link{read_msp_feathers}}
#' @param comment_function function to use when parsing comments
#' @param attribute_sanitization_table lookup tibble between attribute and
#'   new_attribute
#' @param reduce if TRUE then sample 100 entries for each file type to test
#'   whether attributes can be reliably matched. If this done then attribute
#'   -> new_attribute mappings will be returned rather than the primary output.
#'
#' @return a tibble of fragmentation attributes for valid entries which have
#'   a name, mass, precursorMz and adductName
#'
#' @export
format_msp_libraries <- function(msp_compound_attributes,
                                 comment_function,
                                 attribute_sanitization_table,
                                 reduce = FALSE) {
  if (reduce) {

    # look at a small number of entries to verify that attributes can be read

    reduced_entries <- msp_compound_attributes %>%
      dplyr::distinct(file_id, entry) %>%
      dplyr::group_by(file_id) %>%
      dplyr::sample_n(min(100, n()))

    msp_compound_attributes <- msp_compound_attributes %>%
      dplyr::semi_join(reduced_entries, by = c("file_id", "entry"))
  }

  # extract information from comments

  comments <- msp_compound_attributes %>%
    dplyr::filter(attribute == "Comment") %>%
    # parse comments into a tibble
    dplyr::mutate(Comment = purrr::map(value, comment_function)) %>%
    dplyr::select(entry, file_id, Comment) %>%
    tidyr::unnest(Comment)

  # add comment attributes to existing attributes

  library_attributes <- do.call(
    dplyr::bind_rows,
    list(
      msp_compound_attributes %>%
        dplyr::mutate(value = trimws(value)),
      comments
    )
  ) %>%
    dplyr::mutate(lc_variable = tolower(attribute)) %>%
    # standardize attributes of interest
    dplyr::inner_join(
      attribute_sanitization_table,
      by = c("lc_variable" = "attribute")
    ) %>%
    dplyr::select(file_id, entry, attribute = new_attribute, value = value)

  # add exact mass converted from formulas

  library_formulas <- library_attributes %>%
    dplyr::filter(attribute %in% dplyr::filter(
      attribute_sanitization_table,
      new_attribute == "formula"
    )$attribute)

  formula_exact_masses <- clamr::formula_monoisotopic_mass(
    library_formulas$value
  )

  library_formulas <- library_formulas %>%
    dplyr::inner_join(formula_exact_masses, by = c("value" = "formula")) %>%
    dplyr::select(-value) %>%
    dplyr::rename(value = exact_mass) %>%
    dplyr::mutate(
      value = as.character(value),
      attribute = "mass"
    )

  library_attributes <- library_attributes %>%
    dplyr::bind_rows(library_formulas)

  # add precursorMz values from exact mass and adduct

  adduct_masses <- library_attributes %>%
    dplyr::filter(attribute == "adductName") %>%
    dplyr::distinct(value) %>%
    {
      clamr::adduct_mass(.$value)
    }

  # keep at most one entry for most fields, and retain all classes
  formatted_attributes <- library_attributes %>%
    dplyr::filter(attribute != "class") %>%
    dplyr::group_by(file_id, entry, attribute) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    tidyr::spread(attribute, value) %>%
    # add classes
    dplyr::left_join(library_attributes %>%
      dplyr::filter(attribute == "class") %>%
      tidyr::nest(classes = -c(entry, file_id)),
    by = c("file_id", "entry")
    ) %>%
    # add ion adduct mass and z
    dplyr::left_join(adduct_masses, by = c("adductName" = "adduct")) %>%
    # update formatting
    dplyr::mutate(
      .entry_id = 1:n(),
      dataset_id = "CalicoLipids",
      mass = as.numeric(mass),
      mode = dplyr::case_when(
        stringr::str_detect(pattern = "(pos)|P", file_id) ~ "positive",
        stringr::str_detect(pattern = "(neg)|N", file_id) ~ "negative"
      ),
      precursorCharge = abs(adduct_charge),
      precursorMz = dplyr::case_when(
        !is.na(precursorMz) ~ as.numeric(precursorMz),
        TRUE ~ (mass + adduct_mass) / precursorCharge
      ),
      energy = NA_real_,
      smiles = NA_character_,
      inchikey = NA_character_,
      compoundName = stringr::str_replace(
        name,
        " \\[[A-Za-z0-9+-]+\\][-+]$",
        ""
      )
    )

  # check for existence of fragmentationData

  valid_entries <- formatted_attributes %>%
    # filter entries with too much missing information for databasing
    dplyr::filter(
      !is.na(mass),
      !is.na(precursorMz),
      !is.na(precursorCharge),
      !is.na(name),
      !is.na(adductName),
      !is.na(mode)
    )

  if (nrow(valid_entries) != nrow(formatted_attributes)) {
    removed_entries <- formatted_attributes %>%
      dplyr::anti_join(valid_entries, by = ".entry_id")

    invalid_attributes <- formatted_attributes %>%
      dplyr::select(file_id, mass, precursorMz, name, adductName, mode) %>%
      tidyr::gather(attribute, value, -file_id) %>%
      dplyr::filter(is.na(value)) %>%
      dplyr::count(file_id, attribute) %>%
      dplyr::mutate(tag = glue::glue("no {attribute} ({n})")) %>%
      dplyr::group_by(file_id) %>%
      dplyr::summarize(problems = paste(tag, collapse = ", "))

    invalid_attributes %>%
      glue::glue_data("{file_id} has the following problems: {problems}") %>%
      paste(collapse = "\n  ") %>%
      warning()
  }

  if (reduce) {
    output <- do.call(
      dplyr::bind_rows,
      list(
        msp_compound_attributes %>%
          dplyr::mutate(value = trimws(value)),
        comments
      )
    ) %>%
      dplyr::mutate(lc_variable = tolower(attribute)) %>%
      dplyr::left_join(
        attribute_sanitization_table,
        by = c("lc_variable" = "attribute")
      ) %>% # for debugging
      dplyr::select(file_id, entry, attribute, new_attribute, value) %>%
      dplyr::semi_join(removed_entries, by = c("file_id", "entry"))
  } else {
    output <- valid_entries
  }

  return(output)
}
