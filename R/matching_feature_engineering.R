#' Meta Featurize Scores
#'
#' Take a set of possible matches between groups and ions and generate meta features which compare
#' summarize the # of matches scored for each ion, score compared to the top ion match, etc.
#'
#' @param standard_match_scores comparison_features with a column for score based on score_standard_matches match_scores.
#'
#' @return summary features for each possible match between ionIds and groupIds
meta_featurize_scores <- function(standard_match_scores) {

  # comparison to top
  groupId_summaries <- standard_match_scores %>%
    dplyr::group_by(groupId) %>%
    dplyr::summarize(n_ion_candidates = dplyr::n())

  ionId_summaries <- standard_match_scores %>%
    dplyr::group_by(ionId) %>%
    dplyr::summarize(n_ion_calls = dplyr::n())

  match_meta_summaries <- standard_match_scores %>%
    dplyr::select(groupId, ionId, score) %>%
    dplyr::arrange(desc(score)) %>%
    # compare ranks_by_ion
    dplyr::group_by(groupId) %>%
    # NAs will later be replaced 0
    # single_groupId_candidate included to identify entries which will be replaced with 0
    dplyr::mutate(
      single_groupId_candidate = ifelse(dplyr::n() == 1, 1, 0),
      top_groupId_score_diff = score - score[1],
      second_groupId_score_diff = score - score[2]
    ) %>%
    dplyr::group_by(ionId) %>%
    dplyr::mutate(
      single_ionId_candidate = ifelse(dplyr::n() == 1, 1, 0),
      top_ionId_score_diff = score - score[1],
      second_ionId_score_diff = score - score[2]
    ) %>%
    dplyr::ungroup()

  match_meta_summaries %>%
    dplyr::left_join(groupId_summaries, by = "groupId") %>%
    dplyr::left_join(ionId_summaries, by = "ionId")
}
