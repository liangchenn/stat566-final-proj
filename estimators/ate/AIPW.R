#' README:
#' -------
#' - author: Jian Kang, Liang-Cheng Chen
#' - date: 2026-06-03
#'
#' Desc:
#' -------
#' This file defines the cross-fitted augmented inverse-probability-weighted
#' (AIPW) estimator for the ATE. The implementation follows notebooks/AIPW.Rmd:
#' logistic outcome and propensity-score models are fit on K - 1 folds and
#' evaluated on the held-out fold.


# Packages ------------------------------------------------------------------------------------
library(data.table)
source("estimators/utils.R")


# Estimator -----------------------------------------------------------------------------------

estimate_ate_aipw <- function(data,
                              outcome_var,
                              treat_var,
                              covariates,
                              estimand = "ATE",
                              covariate_set = NA_character_,
                              trim = c(0.02, 0.98),
                              K = 5,
                              seed = 2026,
                              conf_level = 0.95,
                              return_details = FALSE
) {
    required_vars <- unique(c(outcome_var, treat_var, covariates))
    missing_vars <- setdiff(required_vars, names(data))

    if (length(missing_vars) > 0) {
        stop(sprintf("Missing variables: %s", paste(missing_vars, collapse = ", ")))
    }
    if (estimand != "ATE") {
        stop("estimate_ate_aipw() only supports estimand = 'ATE'.")
    }
    if (length(trim) != 2 || trim[1] <= 0 || trim[2] >= 1 || trim[1] >= trim[2]) {
        stop("`trim` must contain two ordered values strictly between 0 and 1.")
    }
    if (length(K) != 1 || K < 2 || K != as.integer(K)) {
        stop("`K` must be an integer greater than or equal to 2.")
    }
    if (length(conf_level) != 1 || conf_level <= 0 || conf_level >= 1) {
        stop("`conf_level` must be strictly between 0 and 1.")
    }

    analysis_data <- as.data.frame(data)
    analysis_data <- analysis_data[, required_vars, drop = FALSE]
    analysis_data <- analysis_data[stats::complete.cases(analysis_data), , drop = FALSE]
    n <- nrow(analysis_data)

    if (n < K) {
        stop("`K` cannot exceed the number of complete observations.")
    }

    A <- analysis_data[[treat_var]]
    Y <- analysis_data[[outcome_var]]

    if (!all(A %in% c(0, 1))) {
        stop(sprintf("Treatment variable `%s` must be coded as 0/1.", treat_var))
    }
    if (!is.numeric(Y)) {
        stop(sprintf("Outcome variable `%s` must be numeric.", outcome_var))
    }
    if (!all(Y %in% c(0, 1))) {
        stop(sprintf(
            "Outcome variable `%s` must be binary (0/1) for logistic AIPW.",
            outcome_var
        ))
    }

    outcome_formula <- make_interacted_outcome_formula(
        outcome_var = outcome_var,
        treat_var = treat_var,
        covariates = covariates
    )
    ps_formula <- make_formula(treat_var, covariates)

    set.seed(seed)
    folds <- sample(rep(seq_len(K), length.out = n))
    scores <- numeric(n)

    for (k in seq_len(K)) {
        train <- analysis_data[folds != k, , drop = FALSE]
        test <- analysis_data[folds == k, , drop = FALSE]

        if (length(unique(train[[treat_var]])) < 2) {
            stop(sprintf("Training fold %d does not contain both treatment groups.", k))
        }

        outcome_fit <- stats::glm(
            outcome_formula,
            data = train,
            family = stats::binomial("logit")
        )
        propensity_fit <- stats::glm(
            ps_formula,
            data = train,
            family = stats::binomial("logit")
        )

        test_treated <- test
        test_control <- test
        test_treated[[treat_var]] <- 1
        test_control[[treat_var]] <- 0

        mu1 <- stats::predict(outcome_fit, newdata = test_treated, type = "response")
        mu0 <- stats::predict(outcome_fit, newdata = test_control, type = "response")
        ps <- stats::predict(propensity_fit, newdata = test, type = "response")
        ps <- pmin(pmax(ps, trim[1]), trim[2])

        A_test <- test[[treat_var]]
        Y_test <- test[[outcome_var]]
        scores[folds == k] <- (mu1 - mu0) +
            A_test * (Y_test - mu1) / ps -
            (1 - A_test) * (Y_test - mu0) / (1 - ps)
    }

    estimate <- mean(scores)

    if (!return_details) {
        return(estimate)
    }

    alpha <- 1 - conf_level
    critical_value <- stats::qnorm(1 - alpha / 2)
    std_error <- stats::sd(scores) / sqrt(n)

    return(list(
        estimate = estimate,
        std_error = std_error,
        ci_lower = estimate - critical_value * std_error,
        ci_upper = estimate + critical_value * std_error,
        n = n,
        K = K,
        trim = trim,
        scores = scores,
        folds = folds,
        outcome_var = outcome_var,
        treat_var = treat_var,
        covariate_set = covariate_set,
        covariates = covariates
    ))
}


# Explicit alias for callers that want the method name to emphasize cross-fitting.
estimate_ate_crossfit_aipw <- estimate_ate_aipw
