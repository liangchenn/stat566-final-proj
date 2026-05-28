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

estimate_att_ps_matching <- function(data,
                                     outcome_var,
                                     treat_var,
                                     covariates,
                                     estimand = "ATT",
                                     covariate_set = NA_character_,
                                     trim = c(0.02, 0.98),
                                     M = 1,
                                     replace = TRUE) {
    # estimand <- match_estimand(estimand)
    # if (estimand != "ATT") {
    #     stop("estimate_att_ps_matching() only supports estimand = 'ATT'.")
    # }
    # validate_estimator_input(data, outcome_var, treat_var, covariates)
    # 
    # if (!requireNamespace("Matching", quietly = TRUE)) {
    #     stop("Package `Matching` is required for propensity-score matching.")
    # }

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

    return(matching_fit$est[[1]])
}
