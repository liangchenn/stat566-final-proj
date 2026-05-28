#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-05-28
#'
#' Desc:
#' -------
#' This file defines the shared estimator contract for the project.
#' Every causal estimator should expose the same input pattern and return the
#' same output schema so analysis files can compare estimators directly.
#'
#' Required estimator interface:
#' -----------------------------
#' estimate_xxx <- function(
#'     data,
#'     outcome_var,
#'     treat_var,
#'     covariates,
#'     estimand = c("ATE", "ATT"),
#'     ...
#' )
#'
#' Required estimator output:
#' --------------------------
#' A one-row data.table with columns:
#' - estimator: method name, e.g. "AIPW", "IPW", "Outcome regression"
#' - estimand: target estimand, e.g. "ATE" or "ATT"
#' - estimate: point estimate on the outcome scale
#' - std_error: standard error if available, otherwise NA_real_
#' - ci_lower: lower confidence interval if available, otherwise NA_real_
#' - ci_upper: upper confidence interval if available, otherwise NA_real_
#' - n: sample size used
#' - outcome_var: outcome column used
#' - treat_var: treatment column used
#' - covariate_set: optional adjustment-set name
#' - covariates: comma-separated covariate names
#' - details: optional short note
#'
#' Practical convention:
#' ---------------------
#' Estimator-specific functions should compute a point estimate and call
#' new_estimator_result(). Bootstrap uncertainty should be added by
#' bootstrap_estimator(), which keeps uncertainty calculations consistent.


# Packages ------------------------------------------------------------------------------------
library(data.table)
source("estimators/utils.R")


# Validation ----------------------------------------------------------------------------------

validate_estimator_input <- function(data, outcome_var, treat_var, covariates = character()) {
    required_vars <- unique(c(outcome_var, treat_var, covariates))
    missing_vars <- setdiff(required_vars, names(data))

    if (length(missing_vars) > 0) {
        stop(sprintf("Missing variables: %s", paste(missing_vars, collapse = ", ")))
    }

    if (!all(data[[treat_var]] %in% c(0, 1))) {
        stop(sprintf("Treatment variable `%s` must be coded as 0/1.", treat_var))
    }

    if (!is.numeric(data[[outcome_var]])) {
        stop(sprintf("Outcome variable `%s` must be numeric.", outcome_var))
    }

    invisible(TRUE)
}

match_estimand <- function(estimand) {
    match.arg(estimand, choices = c("ATE", "ATT"))
}


# Result objects -------------------------------------------------------------------------------

new_estimator_result <- function(estimator,
                                 estimand,
                                 estimate,
                                 data,
                                 outcome_var,
                                 treat_var,
                                 covariates,
                                 std_error = NA_real_,
                                 ci_lower = NA_real_,
                                 ci_upper = NA_real_,
                                 covariate_set = NA_character_,
                                 details = NA_character_) {
    estimand <- match_estimand(estimand)

    data.table(
        estimator = estimator,
        estimand = estimand,
        estimate = as.numeric(estimate),
        std_error = as.numeric(std_error),
        ci_lower = as.numeric(ci_lower),
        ci_upper = as.numeric(ci_upper),
        n = nrow(data),
        outcome_var = outcome_var,
        treat_var = treat_var,
        covariate_set = covariate_set,
        covariates = paste(covariates, collapse = ", "),
        details = details
    )
}

is_estimator_result <- function(x) {
    required_cols <- c(
        "estimator", "estimand", "estimate", "std_error", "ci_lower",
        "ci_upper", "n", "outcome_var", "treat_var", "covariate_set",
        "covariates", "details"
    )

    is.data.frame(x) && all(required_cols %in% names(x)) && nrow(x) == 1
}

assert_estimator_result <- function(x) {
    if (!is_estimator_result(x)) {
        stop("Estimator must return a one-row result from new_estimator_result().")
    }
    invisible(TRUE)
}


# Bootstrap wrapper ---------------------------------------------------------------------------

bootstrap_estimator <- function(estimator_fn,
                                data,
                                outcome_var,
                                treat_var,
                                covariates,
                                estimand = c("ATE", "ATT"),
                                n_bootstrap = 100,
                                seed = 2026,
                                conf_level = 0.95,
                                ...) {
    estimand <- match_estimand(estimand)
    validate_estimator_input(data, outcome_var, treat_var, covariates)

    point_result <- estimator_fn(
        data = data,
        outcome_var = outcome_var,
        treat_var = treat_var,
        covariates = covariates,
        estimand = estimand,
        ...
    )
    assert_estimator_result(point_result)

    set.seed(seed)
    boot_estimates <- replicate(n_bootstrap, {
        idx <- sample.int(nrow(data), size = nrow(data), replace = TRUE)
        boot_data <- data[idx, , drop = FALSE]

        boot_result <- estimator_fn(
            data = boot_data,
            outcome_var = outcome_var,
            treat_var = treat_var,
            covariates = covariates,
            estimand = estimand,
            ...
        )
        assert_estimator_result(boot_result)
        boot_result$estimate
    })

    alpha <- 1 - conf_level
    point_result[, std_error := stats::sd(boot_estimates, na.rm = TRUE)]
    point_result[, ci_lower := stats::quantile(boot_estimates, alpha / 2, na.rm = TRUE)]
    point_result[, ci_upper := stats::quantile(boot_estimates, 1 - alpha / 2, na.rm = TRUE)]
    point_result[, details := paste0(
        ifelse(is.na(details), "", details),
        ifelse(is.na(details), "", "; "),
        "bootstrap n=", n_bootstrap,
        ", conf_level=", conf_level
    )]

    attr(point_result, "bootstrap_estimates") <- boot_estimates
    point_result
}


# Optional report helpers ---------------------------------------------------------------------

combine_estimator_results <- function(...) {
    results <- list(...)
    invisible(lapply(results, assert_estimator_result))
    rbindlist(results, fill = TRUE)
}

select_analysis_columns <- function(result) {
    result[, .(
        estimator,
        estimand,
        estimate,
        std_error,
        ci_lower,
        ci_upper,
        n,
        covariate_set,
        details
    )]
}
