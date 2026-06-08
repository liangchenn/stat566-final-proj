#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-06-07
#'
#' Desc:
#' -------
#' This file defines the 1:1 nearest-neighbor Mahalanobis matching
#' estimator for the ATE, following the structure of estimators/ate/matching.R.

# Packages ------------------------------------------------------------------------------------
library(data.table)
source("estimators/utils.R")


# Estimator -----------------------------------------------------------------------------------

estimate_ate_mahalanobis_matching <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    estimand = "ATE",
    M = 1,
    replace = TRUE,
    return_deatails = FALSE
) {
    # 1. use all covariates given to compute distance
    X <- stats::model.matrix(
        as.formula(paste("~", make_rhs(covariates))),
        data = data
    )

    # X[, .("(Intercept)") := NULL] # only for data.table
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]

    # 2. Matching with X
    matching_fit <- Matching::Match(
        Y = data[[outcome_var]],
        Tr = data[[treat_var]],
        X = X,
        M = M,
        estimand = estimand,
        replace = replace,
        Weight = 2
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
