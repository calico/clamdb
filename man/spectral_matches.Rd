% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matching_utils.R
\name{spectral_matches}
\alias{spectral_matches}
\title{Determine best matches for library fragment peaks to experimental fragment peaks.}
\usage{
spectral_matches(
  library_peaks,
  exp_peaks,
  maxMzDiffInAmu = 0.01,
  min_intensity_frac_exp = 1e-04
)
}
\arguments{
\item{library_peaks}{list of m/zs corresponding to fragment spectra, ordered in increasing order}

\item{exp_peaks}{list of m/zs correspond to experimental spectra, ordered in increasing m/z}

\item{maxMzDiffInAmu}{maximum m/z difference between two fragments allowable for matching}

\item{min_intensity_frac_exp}{smallest allowable intensity fraction of max intensity for an experimental fragment.
If the fragment is below the threshold, the peak will be excluded.}
}
\value{
rank list of integers asignments of parent peaks to distinct ms2 peak groups
rank[library_position] = experimental_position
-1 if no match

Note that `findFragPairsGreedyMz` is an Rcpp function, needs to be sourced via

Rcpp::sourceCpp("lib/cpp/clamdb/Rcpp_matching_utils.cpp")
}
\description{
Determine best matches for library fragment peaks to experimental fragment peaks.
}
