context("Match entries using names and structure")

test_that("Match compounds by names and aliases", {

  query_names <- tibble::tibble(name = LETTERS[1:10])
  target_names <- tibble::tibble(name = c("A", "a", "B", "b", "C")) %>%
    dplyr::mutate(compoundId = 1:dplyr::n())
  target_aliases <- tibble::tribble(~ compoundId, ~ alias,
                                    1, "A",
                                    3, "E",
                                    6, "E",
                                    7, "F")

  expect_equal(match_structures(query_names, target_names, match_var_order = "name", drop_ambiguous = TRUE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 1)

  expect_equal(match_structures(query_names, target_names, match_var_order = "name", drop_ambiguous = FALSE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 5)

  expect_equal(match_structures(query_names, target_names, target_aliases, match_var_order = c("name", "alias"), drop_ambiguous = TRUE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 3)

  expect_equal(match_structures(query_names, target_names, target_aliases, match_var_order = c("name", "alias"), drop_ambiguous = FALSE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 8)
  
  query_compounds <- tibble::tibble(name = c("Testosterone", "Epitestosterone"))
  target_compounds <- tibble::tibble(systematicCompoundId = 1, name = "Testosterone")
  matches <- match_structures(
    query_compounds,
    target_compounds,
    match_var_order = "name"
    ) %>%
    dplyr::filter(!is.na(systematicCompoundId))
  
  expect_equal(
    matches,
    tibble::tibble(
      name = "Testosterone",
      systematicCompoundId = 1,
      match_var = "name"
      ))
})

test_that("Match compounds by structure", {

  query_names <- tibble::tibble(name = c("A", "B", "C", "D", "E"),
                                smiles = c("OCCCC", "COCCC", "CCOCC", "CCCOC", "CCCCO"),
                                inchikey = c("NBSCHQHZLSJFNQ-GASJEMHNSA-N", "HXXFSFRBOHSIMQ-VFUOTHLCSA-N", "BGWGXPAPYGQALX-VRPWFDPXSA-N", "WQZGKKKJIJFFOK-GASJEMHNSA-N", "WQZGKKKJIJFFOK-VANKVMQKSA-N"))
  target_names <- query_names %>%
    dplyr::mutate(compoundId = 1:dplyr::n())

  target_names$name[c(2:5)] <- target_names$smiles[c(3:5)] <- NA
  target_names$inchikey[3] <- "BGWGXPAPYGQALX-STANDIN-N"

  # sequential matching findings everyone
  expect_equal(match_structures(query_names, target_names, match_var_order = c("name", "smiles", "inchikey", "inchikey_connectivity"), drop_ambiguous = TRUE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 5)

  expect_equal(match_structures(query_names, target_names, match_var_order = c("name", "smiles", "inchikey", "inchikey_connectivity"), drop_ambiguous = TRUE)$match_var,
               c("name", "smiles", "inchikey", "inchikey", "inchikey_connectivity"))

  # match with inchikey finds 4
  expect_equal(match_structures(query_names, target_names, match_var_order = "inchikey", drop_ambiguous = TRUE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 4)

  # match with connecitvity uniquely finds three matches
  expect_equal(match_structures(query_names, target_names, match_var_order = "inchikey_connectivity", drop_ambiguous = TRUE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 3)

  # matching with connectivity without requiring uniqueness finds 7 pairs of matches
  expect_equal(match_structures(query_names, target_names, match_var_order = "inchikey_connectivity", drop_ambiguous = FALSE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 7)

  expect_equal(match_structures(query_names, target_names, match_var_order = "name", drop_ambiguous = FALSE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 1)

  expect_equal(match_structures(query_names, target_names, match_var_order = "smiles", drop_ambiguous = FALSE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 2)
  
  
  
  # matching with smiles + names + alias was not returning non-matches #514
  # can't reproduce
  
  target_aliases <- tibble::tribble(~ compoundId, ~ alias,
                                    1, "A",
                                    3, "E",
                                    4, "F")

  match_structures(query_names, target_names, target_aliases, match_var_order = c("smiles", "name", "alias"), drop_ambiguous = TRUE)
  
})

test_that("Match compounds - test handling of NAs", {

  query_names <- tibble::tibble(name = c("G6", "G1P", "R6P", "Glu", "Glu", rep(NA, 5)),
                                smiles = c("OCCCC", "COCCC", "CCOCC", "CCCOC", "CCCCO", rep(NA, 5)),
                                inchikey = c("NBSCHQHZLSJFNQ-GASJEMHNSA-N", "HXXFSFRBOHSIMQ-VFUOTHLCSA-N", "BGWGXPAPYGQALX-VRPWFDPXSA-N", "WQZGKKKJIJFFOK-GASJEMHNSA-N", "WQZGKKKJIJFFOK-VANKVMQKSA-N", rep(NA, 5)))
  target_names <- query_names %>%
    dplyr::mutate(compoundId = 1:dplyr::n())

  # mix of valid and NAs
  expect_equal(match_structures(query_names, target_names, match_var_order = "inchikey", drop_ambiguous = FALSE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 5)

  expect_equal(match_structures(query_names, target_names, match_var_order = "name", drop_ambiguous = FALSE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 7)

  expect_equal(match_structures(query_names, target_names, match_var_order = "smiles", drop_ambiguous = FALSE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 5)

  expect_equal(match_structures(query_names, target_names, match_var_order = "inchikey_connectivity", drop_ambiguous = FALSE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 7)

  expect_equal(match_structures(query_names, target_names, match_var_order = "inchikey_connectivity", drop_ambiguous = TRUE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 3)

  # no non-NA entries
  expect_equal(match_structures(query_names %>% dplyr::filter(is.na(name)), target_names %>% dplyr::filter(is.na(name)), match_var_order = c("name", "smiles", "inchikey", "inchikey_connectivity"), drop_ambiguous = FALSE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 0)

  expect_equal(match_structures(query_names %>% dplyr::filter(is.na(name)), target_names %>% dplyr::filter(is.na(name)), match_var_order = c("name", "smiles", "inchikey", "inchikey_connectivity"), drop_ambiguous = TRUE) %>%
                 dplyr::filter(!is.na(compoundId)) %>%
                 nrow(), 0)

})

