#' Test Collision Energy
#'
#' @param energy_vector a character vector of collision of the form "20,40,80"
test_collision_energy <- function(energy_vector) {
  stopifnot("character" %in% class(energy_vector))

  valid_energy_expr <- "^[0-9]+(,[0-9]+)*$"

  valid_energy_expr_set <- stringr::str_detect(
    energy_vector,
    valid_energy_expr
  )
  valid_energy_unique_values <- purrr::map_lgl(
    energy_vector[valid_energy_expr_set],
    function(x) {
      strsplit(x, split = ",")[[1]] %>%
        {
          length(unique(.)) == length(.)
        }
    }
  )

  if (any(!valid_energy_expr_set) | any(!valid_energy_unique_values)) {
    invalid_entries <- c(
      energy_vector[!valid_energy_expr_set],
      energy_vector[valid_energy_expr_set][!valid_energy_unique_values]
    )

    stop(
      length(invalid_entries),
      " invalid collision energy specifications:
        collision energies should be of the form \"20,30,40\"
        invalid entries are:\n",
      paste(invalid_entries, collapse = "  ")
    )
  }

  invisible(0)
}

test_reticulate_py_dependences <- function(modules) {
  if (length(modules) == 0) {
    return(invisible(0))
  }

  stopifnot(all(class(modules) %in% "character", length(modules) >= 1))

  available_modules <- purrr::map_lgl(modules, reticulate::py_module_available)
  if (all(available_modules)) {
    return(invisible(0))
  }
  stop(
    "required modules ",
    paste(modules[!available_modules], collapse = ", "),
    " not found"
  )
}
