#' Canonicalize Structures
#'
#' Call a conda installation of RDKit using reticulate to canonicalize and
#'   interconvert structures.
#'
#' @param input_vector a character vector of SMILES IDs.
#' @param input_type smiles or inchi
#' @param output_types e.g., smiles, inchi, inchi-key
#' @param conda_env a conda environment where RDKit is installed.
#' @param conda_path a valid path to the conda binary or auto to let reticulate
#'   find a default installation
#' @param remove_salts TRUE/FALSE whether salts should be removed from a
#'   formula during conversion
#'
#' @return a table mapping the smiles_vector (smiles variable) on
#'   canonical SMILES
#'
#' @examples
#' \dontrun{
#' input_vector <- c("O=C(O)C(N)C(O)C", "CC(O)C(N)C(=O)O", "badID")
#' input_type <- "smiles"
#' output_types <- "smiles"
#' conda_env <- "rdkit37"
#'
#' canonicalize_structures(input_vector, input_type, output_types, conda_env)
#' }
#'
#' @export
canonicalize_structures <- function(input_vector,
                                    input_type,
                                    output_types,
                                    conda_env,
                                    conda_path = "auto",
                                    remove_salts = TRUE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop('The "reticulate" package must be installed to use this function',
      call. = FALSE
    )
  }

  reticulate::use_condaenv(conda_env, conda = conda_path, required = TRUE)

  if (reticulate::py_module_available("mschem") == FALSE) {
    stop(
      "The \"mschem\" package was not found in your conda environment
        please install it from inst/python/mschem"
    )
  }

  mschem <- reticulate::import("mschem")

  canonical_IDs <- mschem$rdkitms$canonicalize_input(
    input_vector,
    input_type,
    output_types,
    remove_salts
  )

  tibble::tibble(type = names(canonical_IDs), list_data = canonical_IDs) %>%
    dplyr::mutate(df_data = purrr::map(list_data, generate_rdkit_table)) %>%
    dplyr::select(-list_data) %>%
    tidyr::unnest(df_data)
}

generate_rdkit_table <- function(an_output) {
  an_output %>%
    unlist() %>%
    {
      tibble::tibble(input = names(.), canonical_ID = unname(.))
    }
}

#' Match Structures
#'
#' @param query_compounds a tibble of compounds containing some combination of:
#'   name, inchi, inchikey, smiles.
#' @param target_compound_structures a tibble of compounds containing a
#'   systematicCompoundId and some combination of: name, inchi, inchikey,
#'   smiles.
#' @param target_compound_aliases an optional tibble containing additional
#'   aliases for a compound (1-many) with fields systematicCompoundId and
#'   alias.
#' @param match_var_order a character vector of identifiers that will be
#'   sequentially searched for matches to a query_compound until a match is
#'   found.
#'
#' Options are \itemize{
#'  \item{inchikey: exact match the full inchikey}
#'  \item{inchikey_noproto: inchikey without the protonation layer (this field
#'    will be generated from the inchikey)}
#'  \item{inchikey_connectivity: exact match to the structure layer of the
#'    inchikey (this field will be generated from the inchikey)}
#'  \item{inchi: exact match to inchi}
#'  \item{smiles: exact match to smiles}
#'  \item{name: exact match to a compound's name}
#'  \item{alias: exact match to a compound's aliases (from
#'    target_compound_aliases)}
#'  }
#' @param target_primary_key variable in \code{target_compound_structures} and
#'   \code{target_compound_aliases} if provided which defines a unique
#'   compound.
#' @param drop_ambiguous if drop_ambiguous is TRUE then an ambiguous match of
#'   one query to 2+ distinct compounds in systematic_compounds_con for a
#'   given variable in match_var_order will not be considered a match.
#' If FALSE then multiple database matches may be returned for a given query.
#'
#' @return query_compounds with a systematicCompoundId column indicating which
#'   entry of target_compounds matches.
#'
#' @export
match_structures <- function(query_compounds,
                             target_compound_structures,
                             target_compound_aliases = NULL,
                             match_var_order = c(
                               "inchikey",
                               "inchikey_noproto",
                               "inchikey_connectivity",
                               "name"
                             ),
                             target_primary_key = "systematicCompoundId",
                             drop_ambiguous = TRUE) {
  checkmate::assertDataFrame(query_compounds)
  query_compounds <- query_compounds %>% dplyr::ungroup()
  checkmate::assertDataFrame(target_compound_structures)
  target_compound_structures <- target_compound_structures %>% dplyr::ungroup()

  valid_match_vars <- c(
    "inchikey",
    "inchikey_noproto",
    "inchikey_connectivity",
    "smiles",
    "inchi",
    "name",
    "alias"
  )

  checkmate::assertCharacter(match_var_order, unique = TRUE)
  purrr::walk(match_var_order, checkmate::assertChoice, valid_match_vars)
  checkmate::assertLogical(drop_ambiguous, len = 1)

  # test for necessary structure fields
  structure_vars <- tibble::tribble(
    ~match_var, ~required_var,
    "inchikey", "inchikey",
    "inchikey_noproto", "inchikey",
    "inchikey_connectivity", "inchikey",
    "smiles", "smiles",
    "inchi", "inchi",
    "name", "name"
  )

  missing_query_fields <- structure_vars %>%
    dplyr::filter(match_var %in% match_var_order) %>%
    dplyr::filter(!(required_var %in% colnames(query_compounds)))
  if (nrow(missing_query_fields) != 0) {
    stop(
      paste(unique(missing_query_fields$required_var), collapse = " & "),
      " must be provided in query_compounds when match_var_order contains ",
      paste(unique(missing_query_fields$match_var), collapse = ", ")
    )
  }

  missing_target_fields <- structure_vars %>%
    dplyr::filter(match_var %in% match_var_order) %>%
    dplyr::filter(!(required_var %in% colnames(target_compound_structures)))
  if (nrow(missing_target_fields) != 0) {
    stop(
      paste(unique(missing_target_fields$required_var), collapse = " & "),
      " must be provided in query_compounds when match_var_order contains ",
      paste(unique(missing_target_fields$match_var), collapse = ", ")
    )
  }

  # determine whether target_primary_key is valid
  checkmate::assertChoice(
    target_primary_key,
    colnames(target_compound_structures)
  )

  distinct_primary_keys <- target_compound_structures %>%
    dplyr::distinct(!!rlang::sym(target_primary_key))

  if (nrow(distinct_primary_keys) != nrow(target_compound_structures)) {
    stop("distinct_primary_keys must be unique in target_compound_structures")
  }

  # test that distinct_primary_keys is in target_compound_structures but
  # not query_compounds
  if (any(
    c(distinct_primary_keys, "match_var", ".query_id") %in% colnames(query_compounds)
  )) {
    stop(glue::glue(
      "{target_primary_key}, match_var & .query_id are reserved variables
          in query_compounds, do not include them"
    ))
  }

  # test for requirements if matching aliases
  if ("alias" %in% match_var_order) {
    stopifnot("data.frame" %in% class(target_compound_aliases))

    missing_alias_fields <- setdiff(
      c(target_primary_key, "alias"),
      colnames(target_compound_aliases)
    )
    if (length(missing_alias_fields) != 0) {
      stop(
        paste(missing_alias_fields, collapse = " & "),
        " must be provided in target_compound_aliases"
      )
    }
  }

  if (any(
    c("inchikey_noproto", "inchikey_connectivity") %in%
      colnames(query_compounds)
  )) {
    stop(
      'Do not include "inchikey_noproto" or "inchikey_connectivity" as
        variables in "query_compounds", it will be generated from the
        supplied inchikey'
    )
  }
  if (any(
    c("inchikey_noproto", "inchikey_connectivity") %in%
      colnames(target_compound_structures)
  )) {
    stop(
      'Do not include "inchikey_noproto" or "inchikey_connectivity" as
        variables in "target_compound_structures", it will be generated
        from the supplied inchikey'
    )
  }

  # test match_var_order validity
  # extract inchikey_noproto and inchikey_connectivity from inchikey

  if (
    any(c("inchikey_noproto", "inchikey_connectivity") %in% match_var_order)
  ) {
    # create inchikey-structure
    # extract inchikey_structure for queries
    query_compounds <- query_compounds %>%
      dplyr::mutate(
        inchikey_connectivity = stringr::str_match(inchikey, "^([A-Z]+)-")[, 2],
        inchikey_noproto = stringr::str_match(inchikey, "^([A-Z]+-[A-Z]+)-")[, 2]
      )

    # extract inchikey_structure for db
    target_compound_structures <- target_compound_structures %>%
      dplyr::mutate(
        inchikey_connectivity = stringr::str_match(inchikey, "^([A-Z]+)-")[, 2],
        inchikey_noproto = stringr::str_match(inchikey, "^([A-Z]+-[A-Z]+)-")[, 2]
      )
  }

  # iterate through match variables and try to match currently unmatched compounds
  unmatched_compounds <- query_compounds %>%
    dplyr::mutate(.query_id = 1:dplyr::n())

  compound_matches_list <- list()
  for (a_match_var in match_var_order) {
    if (a_match_var %in% c(
      "inchikey",
      "inchikey_noproto",
      "inchikey_connectivity",
      "inchi",
      "smiles"
    )) {
      compound_matches <- unmatched_compounds %>%
        # removing NA entries
        dplyr::filter(!is.na(!!(rlang::sym(a_match_var)))) %>%
        dplyr::left_join(
          target_compound_structures %>%
            # select match variable and compound id
            dplyr::select(!!!rlang::syms(c(target_primary_key, a_match_var))) %>%
            # remove NA
            dplyr::filter(!is.na(!!(rlang::sym(a_match_var)))),
          by = a_match_var
        )
    } else if (a_match_var %in% c("name", "alias")) {
      # join upper cases for case insensitivity

      target_names <- if (a_match_var == "name") {
        target_compound_structures %>%
          dplyr::select(systematicCompoundId, name)
      } else if (a_match_var == "alias") {
        target_compound_aliases %>%
          dplyr::rename(name = alias)
      } else {
        stop(a_match_var, " not recognized for matching")
      }

      # remove NA names and convert to upper-case for targets
      target_names <- target_names %>%
        dplyr::filter(!is.na(name)) %>%
        dplyr::mutate(uc_name = toupper(stringr::str_trim(name))) %>%
        dplyr::select(-name)

      # remove NA names and convert to upper-case for queries
      compound_matches <- unmatched_compounds %>%
        dplyr::filter(!is.na(name)) %>%
        dplyr::mutate(uc_name = toupper(stringr::str_trim(name))) %>%
        dplyr::left_join(target_names, by = "uc_name") %>%
        dplyr::select(-uc_name)
    } else {
      stop(a_match_var, " not recognized for matching")
    }

    # discard ambiguous annotations if drop_ambiguous is TRUE

    if (drop_ambiguous) {
      compound_matches <- compound_matches %>%
        dplyr::group_by(.query_id) %>%
        dplyr::mutate(systematicCompoundId = ifelse(
          dplyr::n() != 1,
          NA_integer_,
          systematicCompoundId
        )) %>%
        dplyr::group_by(.query_id) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
    }

    matches <- compound_matches %>%
      dplyr::filter(!is.na(systematicCompoundId)) %>%
      dplyr::mutate(match_var = a_match_var)

    compound_matches_list[[a_match_var]] <- matches

    # compounds which still need to be matched
    unmatched_compounds <- unmatched_compounds %>%
      dplyr::anti_join(matches, by = ".query_id")
  }

  compound_matches_list[["unmatched"]] <- unmatched_compounds %>%
    dplyr::mutate(systematicCompoundId = NA_integer_, match_var = "no match")

  out <- compound_matches_list %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(.query_id) %>%
    dplyr::select(-.query_id)

  if (drop_ambiguous && nrow(out) != nrow(query_compounds)) {
    stop(
      "The output has a different number of entries than query_compounds.
       These entries should match if drop_ambiguous is TRUE"
    )
  }

  out
}
