#' Create Consensus Spectra
#'
#' Aggregate a set of MS2 spectra across a set of grouping variables into
#'   consensus spectra.
#'
#' @param mz_features A set of scans for an ion produced by
#'   \code{\link[clamr]{tidy_mzR_extract_mzs}}.
#' @param grouping_vars NULL if no grouping or a character vector of variables
#'   to summarize over when collapsing MS2s (all variables need to be in
#'   \code{ms_scan_header} of \code{mz_features}).
#' @param MSXtol mass tolerance corresponding to the type of scan being
#'   aggregated. Provided by \link[clamr]{build_clamr_config}
#' @param ic_min m/zs which are not present with an average IC of more
#'   \code{ic_min} for any value of \code{grouping_vars} are filtered.
#' @param ion_scan_number_min Filter mzs which are present in fewer than
#'   \code{ion_scan_number_min} scans (across all \code{grouping_vars}). This
#'   filter is applied after filtering for \code{ic_min}.
#' @param plot_type should a plot be generated; options are: none, print and return.
#'
#' @return a list containing:
#' \itemize{
#'   \item{grouping_var_summary - Information about the scans that were
#'     collapsed over grouping_vars}
#'   \item{consensus_scans - masses and ic of consensus MS2 scans}
#'   \item{summary_plot (optionally) - a summary plot will be added if
#'     \code{plot_type} = "return"}
#'  }
#'
#' @details
#' Currently this function only supports a single mass_accuracy_tolerance for
#'   the purposes collapsing spectra, if you are collapsing MS1 and MS2
#'   spectra, call this function twice with different tolerance parameters.
#'
#' @export
create_consensus_spectra <- function(mz_features,
                                     grouping_vars = NULL,
                                     MSXtol,
                                     ic_min = 1000,
                                     ion_scan_number_min = 1L,
                                     plot_type = "none") {
  checkmate::assertClass(mz_features, "list")
  stopifnot(all(c("fragmentation_events", "fragmentation_data") %in% names(mz_features)))

  stopifnot(class(grouping_vars) %in% c("NULL", "character"))
  missing_grouping_vars <- grouping_vars[!(grouping_vars %in% colnames(mz_features$fragmentation_events))]
  if (length(missing_grouping_vars) != 0) {
    stop("grouping_vars ", paste(missing_grouping_vars, collapse = " & "), " were not present in scan_header")
  }

  checkmate::assertNumber(ic_min, lower = 0)
  checkmate::assertInteger(ion_scan_number_min, len = 1, lower = 1)
  checkmate::assertChoice(plot_type, c("none", "print", "return"))

  if (nrow(mz_features$fragmentation_data) == 0) {
    output <- list(
      grouping_var_summary = tibble::tibble(),
      consensus_scans = tibble::tibble()
    )
    if (plot_type == "return") {
      output$summary_plot <- NULL
    }
    return(output)
  }

  grouped_fragments <- mz_features$fragmentation_data %>%
    # add grouping variables to scan data
    dplyr::left_join(
      mz_features$fragmentation_events %>%
        dplyr::select(!!!c("scan", rlang::syms(grouping_vars))),
      by = "scan"
    ) %>%
    dplyr::arrange(mz) %>%
    # bin mzs based on mass tolerance
    dplyr::mutate(mz_bin = clamr::find_mz_jumps(
      .$mz,
      MSXtol$tol,
      MSXtol$absolute_or_relative,
      collapse = FALSE
    ))

  # Calculate cosine similarities between scans which share grouping_vars as a measure of signal quality.
  scan_similarity <- grouped_fragments %>%
    tidyr::nest(scan_data = c(!!!rlang::syms(setdiff(colnames(grouped_fragments), grouping_vars)))) %>%
    dplyr::mutate(scan_cosine_similarity = purrr::map_dbl(scan_data, cosine_sim_scans)) %>%
    dplyr::select(-scan_data)

  scan_counts <- grouped_fragments %>%
    dplyr::distinct(!!!quos(!!!rlang::syms(c("scan", grouping_vars)))) %>%
    dplyr::count(!!!quos(!!!rlang::syms(grouping_vars)), name = "n_scans")

  if (class(grouping_vars) == "NULL") {
    grouping_var_summary <- scan_similarity %>% dplyr::mutate(n_scans = scan_counts$n_scans)
  } else {
    grouping_var_summary <- dplyr::left_join(scan_similarity, scan_counts, by = grouping_vars)
  }

  consensus_scans <- grouped_fragments %>%
    # define mz by the bin average
    dplyr::group_by(mz_bin) %>%
    dplyr::mutate(mz = mean(mz * ic) / mean(ic)) %>%
    # sum over multiple measurements of a scan in the same bin and over all scans that share grouping_vars
    dplyr::group_by(!!!quos(!!!rlang::syms(c("mz_bin", grouping_vars)))) %>%
    dplyr::summarize(mz = mz[1], ic_sum = sum(ic) / dplyr::n(), .groups = "drop")

  # add # of scans corresponding to grouping_vars
  if (class(grouping_vars) == "NULL") {
    consensus_scans <- consensus_scans %>% dplyr::mutate(n_scans = grouping_var_summary$n_scans)
  } else {
    consensus_scans <- dplyr::left_join(consensus_scans, grouping_var_summary %>% dplyr::select(-scan_cosine_similarity), by = grouping_vars)
  }

  consensus_scans <- consensus_scans %>%
    dplyr::mutate(ic_mean = ic_sum / n_scans) %>%
    # remove ions which are low signal for all values of grouping_vars
    dplyr::group_by(mz_bin) %>%
    dplyr::filter(any(ic_mean >= ic_min)) %>%
    # remove ions present in few samples
    dplyr::filter(sum(n_scans) >= ion_scan_number_min) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(!!!quos(!!!rlang::syms(c(grouping_vars, "mz")))) %>%
    dplyr::select(!!!quos(!!!rlang::syms(c(grouping_vars, "mz", "ic_mean")))) %>%
    dplyr::rename(ic = ic_mean)

  if (plot_type != "none") {
    summary_plot <- consensus_scans %>%
      dplyr::arrange(mz) %>%
      dplyr::mutate(mz = signif(mz, 6)) %>%
      dplyr::mutate(mz = factor(mz, levels = unique(.$mz))) %>%
      ggplot(aes(x = mz, y = ic)) +
      geom_bar(width = 0.8, stat = "identity") +
      scale_y_continuous("Ion Current (averaged over scans)", trans = "sqrt") +
      scale_x_discrete("M/Z") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))

    if (class(grouping_vars) != "NULL") {
      summary_plot <- summary_plot + facet_grid(as.formula(paste0(paste(grouping_vars, collapse = " + "), " ~ .")), scales = "free_y")
    }
  }

  output <- list(
    grouping_var_summary = grouping_var_summary,
    consensus_scans = consensus_scans
  )

  if (plot_type != "none") {
    if (plot_type == "print") {
      print(summary_plot)
    } else if (plot_type == "return") {
      output$summary_plot <- summary_plot
    } else {
      stop('Invalid "plot_type"')
    }
  }

  output
}

#' Cosine Similarity Scans
#'
#' Calculates the median cosine similarity across all pairs of different scans
#'
#' @param scan_data from \code{read_ms_file_mzR}
#'
#' @return a length one numeric
cosine_sim_scans <- function(scan_data) {
  if (length(unique(scan_data$scan)) == 1) {
    return(NA_real_)
  }

  reshape2::acast(scan_data, formula = scan ~ mz_bin, value.var = "ic", fun.aggregate = sum) %>%
    cosine_similarity() %>%
    {
      .[upper.tri(.)]
    } %>%
    median()
}

cosine_similarity <- function(x) {
  x %*% t(x) / (sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
}

#' Summarize Loop Injection
#'
#' This function wraps \code{clamr::tidy_mzR_extract_mzs} and
#'   \code{\link{create_consensus_spectra}} to generate a summary of MS2
#'   spectra at various collision energies.
#' 
#' @param ms_file_path full path to MS sample file (mzML or mzXML).
#' @param mz_compound standard information.
#' @param clamr_config configuration file for clamr search.
#' @param minTIC Minimum MS2 TIC to consider for filtering MS2 scans.
#' 
#' @param ... additional arguments to pass to
#'   \code{\link{create_consensus_spectra}}
#'
#' @export
summarize_loop_injection <- function(ms_file_path,
                                     mz_compound,
                                     minTIC = 0,
                                     clamr_config,
                                     ...) {
  dots <- list(...)

  checkmate::assertNumber(minTIC, lower = 0)
  clamr::require_tolerances(clamr_config, c(1L, 2L))

  mz_features <- clamr::tidy_mzR_from_msfile(ms_file_path) %>%
    clamr::tidy_mzR_extract_mzs(
      clamr_config,
      mz_compounds = mz_compound
    )

  # pull out arguments to summarize_loop_injection provided with ...
  summarize_loop_injection_extra_args <- dots[intersect(
    names(formals(summarize_loop_injection)),
    names(dots)
  )]

  # filter to scans with minimum minTIC
  # minTIC is used instead of precursorInsensity because TIC is determined for
  #   each scan while precursorIntensity is from the last MS1 event and
  #   depending on the experiment MS1 events may be rare.

  mz_features$fragmentation_events <- mz_features$fragmentation_events %>%
    filter(totIonCurrent >= minTIC)

  mz_features$fragmentation_data <- mz_features$fragmentation_data %>%
    dplyr::semi_join(mz_features$fragmentation_events, by = "scan")

  # group by collision energy and average spectra
  collapsed_spectra <- do.call(
    create_consensus_spectra,
    append(
      list(
        mz_features = mz_features,
        grouping_vars = "collisionEnergy",
        MSXtol = clamr_config$MS2tol
      ),
      summarize_loop_injection_extra_args
    )
  )

  collapsed_spectra$eic <- mz_features$eic
  collapsed_spectra
}
