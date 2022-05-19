#' Score Standard Matches
#'
#' @param comparison_features Output of \link{compare_aligned_features}.
#' @param grouping_vars A character vector of variables which define a unique peak group entry.
#' @param score_method A method for scoring the features.
#' @param ... additional arguments to pass to scoring methods
#'
#' @return comparison_features with \code{score} and \code{is_match} added as variables.
#'
#' @export
score_standard_matches <- function(comparison_features, grouping_vars = "groupId", score_method, ...) {
  dots <- list(...)

  stopifnot(all(grouping_vars %in% colnames(comparison_features)))

  score_method_tbl <- tibble::tribble(
    ~short_name, ~method,
    "baseline", "score_standard_matches_baseline",
    "baseline_lipid", "score_standard_matches_baseline_lipids",
    "lasso", "score_standard_matches_lasso",
    "lasso_svm", "score_standard_matches_lasso_svm",
    "forest", "score_standard_matches_random_forest",
    "custom_function", "score_standard_matches_custom_function"
  )

  checkmate::assertString(score_method)
  if (!(score_method %in% score_method_tbl$short_name)) {
    stop(score_method, " is not a valid scoring method; valid methods are ", paste(score_method_tbl$short_name, collapse = ", "))
  }

  score_method_call <- score_method_tbl$method[score_method_tbl$short_name == score_method]

  # define method args and call score_method
  method_args <- dots[intersect(names(dots), names(formals(score_method_call)))]
  match_scores <- do.call(score_method_call, append(list(comparison_features = comparison_features), method_args))

  top1_matches <- match_scores %>%
    dplyr::group_by(!!!rlang::syms(grouping_vars)) %>%
    dplyr::arrange(desc(score)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  top1_discoveries <- top1_matches %>%
    dplyr::filter(is_match)

  list(
    match_scores = match_scores,
    top1_matches = top1_matches,
    top1_discoveries = top1_discoveries
  )
}

#' Score Standard Matches Baseline
#'
#' @inheritParams score_standard_matches
#' @param MS1sd SD of MS1 tolerance
#' @param RTsd SD of retention time tolerance
score_standard_matches_baseline <- function(comparison_features, MS1sd = 1e-6, RTsd = 0.5) {
  stopifnot(all(c("peakMz_diff", "peakMz_mean", "peakRt_diff_observed", "matched_single_energy_robust_cosine") %in% colnames(comparison_features)))

  baseline_cutoff <- dnorm(1, mean = 0, sd = 1, log = TRUE) + dnorm(0, 1, RTsd, log = TRUE) + log(0.8) * 40

  match_scores <- comparison_features %>%
    dplyr::mutate(
      mz_score = dnorm(peakMz_diff / peakMz_mean / 1e-6, mean = 0, sd = 1, log = TRUE),
      rt_score = dnorm(ifelse(is.na(peakRt_diff_observed), 4, peakRt_diff_observed), mean = 0, sd = RTsd, log = TRUE),
      frag_score = log(ifelse(is.na(matched_single_energy_robust_cosine), 0.5, matched_single_energy_robust_cosine)) * 20,
      score = mz_score + rt_score + frag_score,
      is_match = ifelse(score > baseline_cutoff, TRUE, FALSE)
    )

  match_scores
}

#' Score Standard Matches Lipid Baseline (no RT)
#'
#' @inheritParams score_standard_matches_baseline
score_standard_matches_baseline_lipids <- function(comparison_features, MS1sd = 1e-6) {
  stopifnot(all(c("peakMz_diff", "peakMz_mean", "matched_single_energy_fisher.") %in% colnames(comparison_features)))

  baseline_cutoff <- dnorm(1, mean = 0, sd = 1, log = TRUE) + 5 # a reasonable fisher exact -log10 p-value

  match_scores <- comparison_features %>%
    dplyr::mutate(
      mz_score = dnorm(peakMz_diff / peakMz_mean / 1e-6, mean = 0, sd = 1, log = TRUE),
      frag_score = matched_single_energy_fisher,
      score = mz_score + frag_score,
      is_match = ifelse(score > baseline_cutoff, TRUE, FALSE)
    )

  match_scores
}

#' Score Standard Matches Using LASSO Logistic Regression
#'
#' @inheritParams score_standard_matches
#' @param lasso_pars_df regression coefficients for lasso logistic regression
#' @param lasso_co cutoff for match FDR control
#' @param return_marginals TRUE/FALSE whether to add the marginal fit of each covariate
score_standard_matches_lasso <- function(comparison_features,
                                         lasso_pars_df,
                                         lasso_co = -Inf,
                                         return_marginals = FALSE) {
  stopifnot("data.frame" %in% class(comparison_features))
  stopifnot("data.frame" %in% class(lasso_pars_df))
  stopifnot(all(class(lasso_co) == "numeric"), length(lasso_co) == 1)
  stopifnot(all(class(return_marginals) == "logical"), length(return_marginals) == 1, return_marginals %in% c(TRUE, FALSE))

  # check for all variables with non-zero regression weight
  nonzero_predictors <- lasso_pars_df %>%
    dplyr::filter(coef != 0) %>%
    {
      setdiff(.$variable, "intercept")
    }

  missing_predictors <- nonzero_predictors[!(nonzero_predictors %in% colnames(comparison_features))]

  if (length(missing_predictors) != 0) {
    stop(length(missing_predictors), " needed regression features not found when predicting metabolite matches:\n", paste(missing_predictors, collapse = ",  "))
  }

  tidy_predictors <- comparison_features %>%
    dplyr::select(!!!rlang::syms(nonzero_predictors)) %>%
    dplyr::mutate(.row_number = 1:dplyr::n()) %>%
    tidyr::gather("variable", "value", -.row_number) %>%
    dplyr::filter(variable %in% nonzero_predictors) %>%
    dplyr::mutate(value = ifelse(is.na(value), 0, value))

  marginal_prediction <- tidy_predictors %>%
    dplyr::left_join(lasso_pars_df, by = "variable") %>%
    dplyr::mutate(marginal_prediction = coef * value)

  lasso_prediction <- marginal_prediction %>%
    dplyr::group_by(.row_number) %>%
    dplyr::summarize(
      score = sum(marginal_prediction) + lasso_pars_df$coef[lasso_pars_df$variable == "intercept"],
      is_match = ifelse(score > lasso_co, TRUE, FALSE)
    )

  if (return_marginals) {
    marginal_variables <- marginal_prediction %>%
      dplyr::mutate(variable = paste0(".marginal_", variable)) %>%
      dplyr::select(.row_number, variable, marginal_prediction) %>%
      tidyr::spread(variable, marginal_prediction)

    lasso_prediction <- lasso_prediction %>%
      dplyr::left_join(marginal_variables, by = ".row_number")
  }

  comparison_features %>%
    dplyr::mutate(.row_number = 1:dplyr::n()) %>%
    dplyr::left_join(lasso_prediction, by = ".row_number") %>%
    dplyr::select(-.row_number)
}

#' Score Standard Matches with LASSO and SVM
#'
#' @inheritParams score_standard_matches_lasso
#' @inheritParams score_meta_svm
#' @param svm_co cutoff for match FDR control
score_standard_matches_lasso_svm <- function(comparison_features,
                                             lasso_pars_df,
                                             svm_model,
                                             svm_co = 0.5) {
  stopifnot("data.frame" %in% class(comparison_features))
  stopifnot("data.frame" %in% class(lasso_pars_df))
  stopifnot("svm" %in% class(svm_model))
  stopifnot("numeric" %in% class(svm_co), length(svm_co) == 1)

  lasso_match_scores <- score_standard_matches_lasso(
    comparison_features,
    lasso_pars_df
  )

  meta_features <- meta_featurize_scores(lasso_match_scores) %>%
    dplyr::rename(lasso_score = score)

  meta_features$score <- score_meta_svm(meta_features, svm_model)
  meta_features$is_match <- meta_features$score >= svm_co

  comparison_features %>%
    dplyr::left_join(meta_features, by = c("groupId", "ionId"))
}

#' Score Meta Features using SVM
#'
#' @param meta_features tibble generated from \code{\link{meta_featurize_scores}}
#' @param svm_model parameters of \code{e1071} SVM
score_meta_svm <- function(meta_features,
                           svm_model) {
  stopifnot("data.frame" %in% class(meta_features))
  stopifnot("svm" %in% class(svm_model))

  # check for all variables with non-zero regression weight
  svm_vars <- colnames(svm_model$SV)

  missing_predictors <- svm_vars[!(svm_vars %in% colnames(meta_features))]

  if (length(missing_predictors) != 0) {
    stop(length(missing_predictors), " needed meta features not found when using an SVM to refine matches:\n", paste(missing_predictors, collapse = ",  "))
  }

  x <- meta_features[, svm_vars] %>%
    as.matrix()
  x[is.na(x)] <- 0

  if (!all(colnames(x) == svm_vars)) {
    stop("mismatched matrix and svm variables")
  }

  fitted_svm_match_pr <- predict(svm_model, x, probability = TRUE)
  attr(fitted_svm_match_pr, "probabilities")[, "+"]
}

#' Score Standard Matches with Random Forest
#'
#' @inheritParams score_standard_matches
#' @param forest_model random forest object from \code{randomForest}
#' @param forest_co cutoff for match FDR control
score_standard_matches_random_forest <- function(comparison_features,
                                                 forest_model,
                                                 forest_co = 0.5) {
  stopifnot("data.frame" %in% class(comparison_features))
  stopifnot("randomForest" %in% class(forest_model))
  stopifnot("numeric" %in% class(forest_co), length(forest_co) == 1)
  
  # variables containing in random forest when it was created
  forest_model_predictors <- attr(forest_model$terms, "term.labels")

  missing_forest_model_predictors <- setdiff(forest_model_predictors, colnames(comparison_features))
  if (length(missing_forest_model_predictors) != 0) {
    stop(length(missing_forest_model_predictors), " variables are missing from comparison_features which were present in the random forest classifier
          so this model cannot score metabolite matches. The missing variables are:
          ", paste(missing_forest_model_predictors, collapse = ", "), "\n")
  }

  reduced_comparison_features <- comparison_features %>%
    dplyr::select(!!!forest_model_predictors) %>%
    dplyr::mutate_all(function(x) {
      ifelse(is.na(x), 0, x)
    })

  comparison_features %>%
    dplyr::mutate(
      score = predict(forest_model, newdata = reduced_comparison_features, type = "prob") %>%
        {
          .[, "+"]
        } %>%
        unname(),
      is_match = ifelse(score >= forest_co, TRUE, FALSE)
    )
}

#' Score Standard Matches with an arbitrary function and cutoff
#'
#' @inheritParams score_standard_matches
#' @param custom_function a user-specified function
#' @param cutoff cutoff for match FDR control
score_standard_matches_custom_function <- function(comparison_features,
                                                   custom_function,
                                                   cutoff) {
  stopifnot(class(custom_function) == "function")
  stopifnot(class(cutoff) == "numeric", length(cutoff) == 1)

  custom_function_args <- names(formals(custom_function))
  missing_custom_function_args <- setdiff(custom_function_args, colnames(comparison_features))
  if (length(missing_custom_function_args) != 0) {
    stop("custom_function contains variables which are missing from \"comparision_features\":", paste(missing_custom_function_args, collapse = ", "))
  }

  # turn all args into a list
  custom_function_scores <- comparison_features %>%
    dplyr::select(!!!custom_function_args) %>%
    as.list() %>%
    # apply args with custom_function
    purrr::pmap_dbl(custom_function)

  comparison_features %>%
    dplyr::mutate(
      score = custom_function_scores,
      is_match = ifelse(score >= cutoff, TRUE, FALSE)
    )
}
