#' README:
#' -------
#' - author: Jian Kang, Liang-Cheng Chen
#' - date: 2026-06-08
#'
#' Desc:
#' -------
#' This file defines the augmented inverse-probability-weighted (AIPW)
#' estimator for the ATT.
#'
#' estimate_att_aipw() is the base estimator for one fit/evaluate pass.
#' Bootstrapping can call it directly. Cross-fitting can call it fold by fold
#' with held-out prediction data and combine the returned scores.

# Packages ------------------------------------------------------------------------------------
library(data.table)
source("estimators/utils.R")


# Base estimator ------------------------------------------------------------------------------

estimate_att_aipw <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    estimand = "ATT",
    prediction_data = NULL,
    trim = c(0.02, 0.98),
    return_details = FALSE
) {
    required_vars <- unique(c(outcome_var, treat_var, covariates))

    train_data <- as.data.frame(data)
    train_data <- train_data[, required_vars, drop = FALSE]

    if (is.null(prediction_data)) {
        prediction_data <- train_data
    } else {
        prediction_data <- as.data.frame(prediction_data)
        prediction_data <- prediction_data[, required_vars, drop = FALSE]
    }

    control_train_data <- train_data[
        train_data[[treat_var]] == 0,
        ,
        drop = FALSE
    ]

    outcome_fit <- glm(
        make_formula(outcome_var, covariates),
        data = control_train_data,
        family = binomial("logit")
    )
    propensity_fit <- glm(
        make_formula(treat_var, covariates),
        data = train_data,
        family = binomial("logit")
    )

    mu0_hat <- stats::predict(
        outcome_fit,
        newdata = prediction_data,
        type = "response"
    )
    ps_hat <- stats::predict(
        propensity_fit,
        newdata = prediction_data,
        type = "response"
    )

    A <- prediction_data[[treat_var]]
    Y <- prediction_data[[outcome_var]]
    ps_hat <- pmin(pmax(ps_hat, trim[1]), trim[2])

    scores <- A *
        (Y - mu0_hat) -
        (1 - A) * ps_hat / (1 - ps_hat) * (Y - mu0_hat)
    tau_hat <- mean(scores) / mean(A)

    if (!return_details) {
        return(tau_hat)
    }

    result <- list(
        estimate = tau_hat,
        scores = scores,
        mu0_hat = mu0_hat,
        ps_hat = ps_hat
    )

    return(result)
}
