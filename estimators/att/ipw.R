#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-05-28
#'
#' Desc:
#' -------
#' This file defines the Hajek-stabilized IPW estimator for the ATT.
#' The propensity score is fit with logistic regression using the supplied
#' covariates, then control outcomes are reweighted to the treated population.

# Packages ------------------------------------------------------------------------------------
library(data.table)
source("estimators/utils.R")


# Estimator -----------------------------------------------------------------------------------

estimate_att_hajek_ipw <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    estimand = "ATT",
    trim = c(0.02, 0.98),
    return_details = FALSE
) {
    # extract variables
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]

    # propensity score
    ps <- fit_propensity_score(
        data = data,
        treat_var = treat_var,
        covariates = covariates,
        trim = trim
    )

    control_weight <- (1 - A) * ps / (1 - ps)

    mu1 <- mean(Y[A == 1])
    mu0_treated <- sum(control_weight * Y) / sum(control_weight)
    tau_hat <- mu1 - mu0_treated

    if (!return_details) {
        return(tau_hat)
    }

    result <- list(
        estimate = tau_hat,
        ps = ps,
        control_weight = control_weight
    )

    return(result)
}
