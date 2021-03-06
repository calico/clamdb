% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matching_feature_engineering.R
\name{meta_featurize_scores}
\alias{meta_featurize_scores}
\title{Meta Featurize Scores}
\usage{
meta_featurize_scores(standard_match_scores)
}
\arguments{
\item{standard_match_scores}{comparison_features with a column for score based on score_standard_matches match_scores.}
}
\value{
summary features for each possible match between ionIds and groupIds
}
\description{
Take a set of possible matches between groups and ions and generate meta features which compare
summarize the # of matches scored for each ion, score compared to the top ion match, etc.
}
