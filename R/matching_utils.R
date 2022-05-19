#' Summarizes of peakgroup overall quality
#'
#' @param aligned_labelled_features tibble of groupId, peaks and scans created by \code{\link{align_features_to_standards}} or \code{link[clamr]{nest_peakgroup_features}}
#' @param missing_peakAreaTop missing value ion count cutoff
#'
#' @return summaries of an experimental peakgroup
compare_aligned_features_peakgroup_qc <- function(aligned_labelled_features, missing_peakAreaTop = 1000) {

  # reduce to scans and peaks of unique peak groups
  unique_peakgroups <- aligned_labelled_features %>%
    dplyr::group_by(groupId) %>%
    dplyr::slice(1) %>%
    dplyr::select(groupId, peaks, scans) %>%
    dplyr::ungroup()

  n_samples <- unique_peakgroups %>%
    tidyr::unnest(peaks) %>%
    dplyr::distinct(sampleId) %>%
    nrow()

  # summarize scans associates with each peakgroup
  associated_scan_counts <- unique_peakgroups %>%
    dplyr::filter(!purrr::map_lgl(scans, is.null)) %>%
    dplyr::mutate(n_scans = purrr::map_int(scans, nrow)) %>%
    dplyr::select(groupId, n_scans)

  associated_scan_summaries <- unique_peakgroups %>%
    dplyr::left_join(associated_scan_counts, by = "groupId") %>%
    dplyr::mutate(n_scans = ifelse(is.na(n_scans), 0L, n_scans)) %>%
    dplyr::mutate(scans_per_sample = n_scans / n_samples) %>%
    # floor max scans per sample to avoid some extreme outliers
    dplyr::mutate(scans_per_sample = pmin(scans_per_sample, 5)) %>%
    dplyr::select(groupId, scans_per_sample) %>%
    dplyr::mutate(scans_per_sample_quad = scans_per_sample^2)

  # summarize peaks associated with each peakgroup

  all_peaks <- unique_peakgroups %>%
    dplyr::select(-scans) %>%
    tidyr::unnest(peaks)

  missing_peaks <- tidyr::expand_grid(
    groupId = unique_peakgroups$groupId,
    sampleId = unique(all_peaks$sampleId)
  ) %>%
    dplyr::anti_join(all_peaks, by = c("groupId", "sampleId")) %>%
    dplyr::mutate(
      peakAreaTop = missing_peakAreaTop,
      quality = 0
    )

  # combine detected peaks with peaks that were missing for some samples
  asssociated_peak_summaries <- dplyr::bind_rows(
    all_peaks,
    missing_peaks
  ) %>%
    dplyr::mutate(is_good_peak = ifelse(quality > 0.5, TRUE, FALSE)) %>%
    dplyr::group_by(groupId) %>%
    dplyr::summarize(
      detected_peak_fraction = sum(quality != 0) / n_samples,
      detected_peak_fraction_quad = detected_peak_fraction^2,
      good_peak_fraction = sum(is_good_peak) / n_samples,
      good_peak_fraction_quad = good_peak_fraction^2,
      max_peak_quality = max(quality),
      mean_peak_quality = mean(quality),
      mean_peak_quality_quad = mean_peak_quality^2,
      sd_peak_quality = sd(quality),
      mean_peakAreaTop = mean(peakAreaTop),
      median_peakAreaTop = median(peakAreaTop),
      cv_peakAreaTop = sd(peakAreaTop) / mean_peakAreaTop,
      log2_mean_peakAreaTop = log2(mean_peakAreaTop),
      log2_mean_peakAreaTop_quad = log2_mean_peakAreaTop^2,
      log2_median_peakAreaTop = log2(median_peakAreaTop),
      log2_median_peakAreaTop_quad = log2_median_peakAreaTop^2
    )

  # return all features for each groupId
  asssociated_peak_summaries %>%
    dplyr::left_join(associated_scan_summaries, by = "groupId")
}

#' Compare peakgroup m/z and rt to standard
#' 
#' @param peaks tibble of peaks information
#' @param stdCompounds standard compounds for spectral matching.
#' 
compare_aligned_features_fullscans <- function(peaks, stdCompounds) {
  stopifnot(nrow(stdCompounds) == 1)

  peaks %>%
    # comparison mz / rt from peaks
    dplyr::summarize(
      peakMz_mean = sum(peakMz * peakAreaTop) / sum(peakAreaTop),
      peakMz_sd = sd(peakMz),
      peakRt_mean = sum(peakRt * peakAreaTop) / sum(peakAreaTop),
      peakRt_sd = sd(peakRt),
      # standard rt
      standard_rt_missing = ifelse(is.na(stdCompounds$rtObserved), TRUE, FALSE),
      # mass contrast
      peakMz_diff = peakMz_mean - stdCompounds$precursorMz,
      peakMz_abs_diff = abs(peakMz_diff),
      peakMz_abs_diff_sqrt = sqrt(peakMz_abs_diff),
      # rt contrast
      peakRt_diff_observed = peakRt_mean - stdCompounds$rtObserved,
      peakRt_abs_diff_observed = abs(peakRt_diff_observed),
      peakRt_abs_diff_observed_sqrt = sqrt(peakRt_abs_diff_observed),
      peakRt_diff_aligned = peakRt_mean - stdCompounds$rtAligned,
      peakRt_abs_diff_aligned = abs(peakRt_diff_aligned),
      peakRt_abs_diff_aligned_sqrt = sqrt(peakRt_abs_diff_aligned)
    )
}

#' Compare experimental precursor and standard fragmentation
#' 
#' @param peaks peaks data for alignment.
#' @param scans MS2 scans associated with peak groups.
#' @param stdFragmentations fragmentation patterns associated with standard compounds.
#' @inheritParams compare_aligned_features
#' 
compare_aligned_features_fragmentations <- function(peaks, scans, stdFragmentations, MS2tol, n_top_spectra_summed = 3L, quality_weights = c("purity" = 2, "quality" = 1), frag_similarity_methods = "cosine",
                                                    remove_precursor = FALSE, minimum_ic_fraction = 0, print_plot = print_plot) {
  stopifnot("integer" %in% class(n_top_spectra_summed), length(n_top_spectra_summed) == 1)
  stopifnot(length(quality_weights) == 2, all(names(quality_weights) == c("purity", "quality")), all(quality_weights > 0))
  stopifnot(all(c("tol", "absolute_or_relative") %in% names(MS2tol)))

  # find high-quality experimental fragmentation events and summarize quality

  # peaks <<- peaks
  # scans <<- scans
  # stdFragmentations <<- stdFragmentations

  if (is.null(scans)) {
    fragmentation_summary <- tibble::tibble(no_experimental_fragmentation = TRUE)
  } else {
    best_peakgroup_fragmentations <- scans %>%
      dplyr::left_join(peaks %>%
        dplyr::select(peakId, quality),
      by = "peakId"
      ) %>%
      # for a given peakgroup, prioritize the best spectra to sum based on MS2 purity (precursorPurity) & peak quality (quality)
      dplyr::mutate(MS2_quality = (precursorPurity * quality_weights["purity"] + quality * quality_weights["quality"]) / sum(quality_weights)) %>%
      dplyr::arrange(desc(MS2_quality)) %>%
      # take up to the top n_top_spectra_summed best MS2s
      dplyr::slice(1:n_top_spectra_summed) %>%
      dplyr::select(precursorMz, precursorPurity, quality, MS2_quality, data)

    fragmentation_summary <- best_peakgroup_fragmentations %>%
      dplyr::summarize(
        mean_frag_peakPurity = mean(precursorPurity),
        mean_frag_peakQuality = mean(quality),
        mean_frag_quality = mean(MS2_quality),
        n_summed_spectra = dplyr::n()
      ) %>%
      dplyr::mutate(no_experimental_fragmentation = FALSE)
  }

  if (is.null(scans) || is.null(stdFragmentations)) {
    fragmentation_summary <- fragmentation_summary %>% dplyr::mutate(no_standard_fragmentation = TRUE)

    if (is.null(scans)) {
      warning("no experimental fragmentation data; featurize_fragmentations will return a minimal feature set\n")
    } else {
      warning("no standard fragmentation data; featurize_fragmentations will return a reduced feature set\n")
    }

    return(fragmentation_summary)
  }

  # extract fragmentationData from high-quality fragmentation events

  best_peakgroup_fragData <- best_peakgroup_fragmentations %>%
    # unpack fragmentation data (and bin fragments for a given peak group based on MS2 tolerance)
    dplyr::mutate(peak_num = 1:dplyr::n()) %>%
    tidyr::nest_legacy(.key = "frag") %>%
    dplyr::mutate(fragData = purrr::map(frag, clamr::extract_and_group_fragments, MS2tol = MS2tol)) %>%
    dplyr::select(-frag) %>%
    tidyr::unnest(fragData)

  precursorMz <- mean(best_peakgroup_fragmentations$precursorMz)

  fragmentationData_summary <- best_peakgroup_fragData %>%
    # sum the top n_top_spectra_summed spectra
    dplyr::group_by(mz_set) %>%
    dplyr::summarize(
      mz = sum(ms2_mz * ms2_ic) / sum(ms2_ic),
      ic = mean(c(ms2_ic, rep(0, nrow(best_peakgroup_fragmentations) - dplyr::n()))),
      mean_ms2_ic_frac = mean(c(ms2_ic_frac, rep(0, nrow(best_peakgroup_fragmentations) - dplyr::n()))),
      cv_ms2_ic_frac = sd(c(ms2_ic_frac, rep(0, nrow(best_peakgroup_fragmentations) - dplyr::n()))) / mean_ms2_ic_frac,
      fractionSamples = dplyr::n() / nrow(best_peakgroup_fragmentations)
    ) %>%
    dplyr::mutate(is_precursor = ifelse(abs(mz - precursorMz) < 0.01, TRUE, FALSE))

  # join experimental and standard spectra

  std_fragData <- stdFragmentations %>%
    tidyr::unnest(fragData) %>%
    dplyr::mutate(is_precursor = ifelse(abs(mz - precursorMz) < 0.01, TRUE, FALSE)) %>%
    tidyr::nest(stdFragData = -energy_class)

  spectra_scores <- std_fragData %>%
    dplyr::mutate(scores = purrr::map(stdFragData,
      clamr::score_fragmentation_similarity,
      an_expFrag = fragmentationData_summary,
      MS2tol = MS2tol,
      frag_similarity_methods = frag_similarity_methods,
      remove_precursor = remove_precursor,
      minimum_ic_fraction = minimum_ic_fraction,
      print_plot = print_plot
    )) %>%
    dplyr::filter(!purrr::map_lgl(scores, is.null))

  # check to make sure that spectra could be faithfully compared (and were not filtered due to minimum_ic_fraction or some other condition)

  if (nrow(spectra_scores) == 0) {
    fragmentation_summary <- fragmentation_summary %>% dplyr::mutate(no_standard_fragmentation = TRUE)
    warning("spectra fragments were filtered; featurize_fragmentations will return a reduced feature set\n")
    return(fragmentation_summary)
  }

  spectra_scores <- spectra_scores %>%
    dplyr::select(-stdFragData) %>%
    tidyr::unnest(scores) %>%
    tidyr::unite(score_type, energy_class, frag_similarity_method) %>%
    tidyr::spread(score_type, method_score) %>%
    dplyr::mutate(no_standard_fragmentation = FALSE)

  dplyr::bind_cols(fragmentation_summary, spectra_scores)
}

spectra_normalize <- function(fragData) {
  fragData %>%
    dplyr::mutate(ic = ic / sum(ic)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ic = round(ic, 6)) %>%
    dplyr::filter(ic != 0)
}


#' Determine best matches for library fragment peaks to experimental fragment peaks.
#'
#' @param library_peaks list of m/zs corresponding to fragment spectra, ordered in increasing order
#' @param exp_peaks list of m/zs correspond to experimental spectra, ordered in increasing m/z
#' @param maxMzDiffInAmu maximum m/z difference between two fragments allowable for matching
#' @param min_intensity_frac_exp smallest allowable intensity fraction of max intensity for an experimental fragment.
#' If the fragment is below the threshold, the peak will be excluded.
#'
#' @return rank list of integers asignments of parent peaks to distinct ms2 peak groups
#' rank[library_position] = experimental_position
#' -1 if no match
#'
#' Note that `findFragPairsGreedyMz` is an Rcpp function, needs to be sourced via
#'
#' Rcpp::sourceCpp("lib/cpp/clamdb/Rcpp_matching_utils.cpp")
#'
#' @export
spectral_matches <- function(library_peaks, exp_peaks, maxMzDiffInAmu = 0.01, min_intensity_frac_exp = 0.0001) {
  # TODO: implement intensity filtering, and possibly sorting (to remove sorting requirement)
  return(findFragPairsGreedyMz(library_peaks$mz, exp_peaks$mz, maxMzDiffInAmu, min_intensity_frac_exp))
}
