#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-05-28
#'
#' Desc:
#' -------
#' This file defines the Hajek-stabilized IPW estimator for the ATE.
#' The implementation follows examples/hw3-ans.Rmd and uses the shared
#' propensity-score helper in estimators/utils.R.


# Packages ------------------------------------------------------------------------------------
library(data.table)
source("estimators/utils.R")


# Estimator -----------------------------------------------------------------------------------

estimate_ate_hajek_ipw <- function(data,
                                   outcome_var,
                                   treat_var,
                                   covariates,
                                   estimand = "ATE",
                                   covariate_set = NA_character_,
                                   trim = c(0.02, 0.98)
) {
    
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]
    ps <- fit_propensity_score(
        data = data,
        treat_var = treat_var,
        covariates = covariates,
        trim = trim
    )

    treated_weight <- A / ps
    control_weight <- (1 - A) / (1 - ps)

    mu1 <- sum(treated_weight * Y) / sum(treated_weight)
    mu0 <- sum(control_weight * Y) / sum(control_weight)
    tau_hat <- mu1 - mu0

    return(tau_hat)
}
