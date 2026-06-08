#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-05-28
#'
#' Desc:
#' -------
#' This file defines the 1:1 nearest-neighbor propensity-score matching
#' estimator for the ATT.


# Packages ------------------------------------------------------------------------------------
library(data.table)
source("estimators/utils.R")


# Estimator -----------------------------------------------------------------------------------

estimate_att_ps_matching <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    estimand = "ATT",
    trim = c(0.02, 0.98),
    M = 1,
    replace = TRUE,
    return_details = FALSE
) {
    # 1. calculate propensity score
    ps <- fit_propensity_score(
        data = data,
        treat_var = treat_var,
        covariates = covariates,
        trim = trim
    )

    matching_fit <- Matching::Match(
        Y = data[[outcome_var]],
        Tr = data[[treat_var]],
        X = ps,
        M = M,
        estimand = estimand,
        replace = replace
    )

    # 3. output
    if (!return_details) {
        return(matching_fit$est[[1]])
    }

    result <- list(
        estimate = matching_fit$est[[1]],
        fit = matching_fit
    )

    return(result)
}
