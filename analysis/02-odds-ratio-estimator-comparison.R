#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-05-28
#'
#' Desc:
#' -------
#' This file compares naive, regression-adjusted, IPW, and matching odds-ratio
#' estimates for RHC effects on 30-day and 180-day mortality.
#'
#' Input
#' -----
#' - data/processed/rhc-processed-data.csv
#' - data/processed/rhc-var-sets.rds
#'
#' Output
#' ------
#' - results/tables/odds-ratio-estimator-comparison.csv
#' - results/figures/odds-ratio-estimator-comparison.png

.FILE_NAME <- "02-odds-ratio-estimator-comparison.R"
cat(sprintf("...Running %s ...", .FILE_NAME))


# Packages ------------------------------------------------------------------------------------
library(data.table)
library(readr)
library(ggplot2)

source("estimators/utils.R")


# Main ----------------------------------------------------------------------------------------

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"
output_table_path <- "results/tables/odds-ratio-estimator-comparison.csv"
output_plot_path <- "results/figures/odds-ratio-estimator-comparison.png"

rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)

outcome_vars <- c("death30d", "death180d")
treat_var <- var_sets$treatment
covariates <- var_sets$adjustment_sets$basic


# Odds-ratio helpers --------------------------------------------------------------------------

make_or <- function(p1, p0, eps = 1e-6) {
    p1 <- pmin(pmax(p1, eps), 1 - eps)
    p0 <- pmin(pmax(p0, eps), 1 - eps)
    (p1 / (1 - p1)) / (p0 / (1 - p0))
}

estimate_naive_or <- function(data, outcome_var, treat_var) {
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]
    make_or(mean(Y[A == 1]), mean(Y[A == 0]))
}

estimate_regression_ate_or <- function(data, outcome_var, treat_var, covariates) {
    fit <- glm(
        make_interacted_outcome_formula(outcome_var, treat_var, covariates),
        data = data,
        family = binomial("logit")
    )

    data1 <- copy(data)
    data0 <- copy(data)
    data1[, (treat_var) := 1]
    data0[, (treat_var) := 0]

    p1 <- mean(predict(fit, newdata = data1, type = "response"))
    p0 <- mean(predict(fit, newdata = data0, type = "response"))
    make_or(p1, p0)
}

estimate_regression_att_or <- function(data, outcome_var, treat_var, covariates) {
    fit <- glm(
        make_interacted_outcome_formula(outcome_var, treat_var, covariates),
        data = data,
        family = binomial("logit")
    )

    treated_data <- data[data[[treat_var]] == 1, ]
    treated_data1 <- copy(treated_data)
    treated_data0 <- copy(treated_data)
    treated_data1[, (treat_var) := 1]
    treated_data0[, (treat_var) := 0]

    p1 <- mean(predict(fit, newdata = treated_data1, type = "response"))
    p0 <- mean(predict(fit, newdata = treated_data0, type = "response"))
    make_or(p1, p0)
}

estimate_ate_ipw_or <- function(data, outcome_var, treat_var, covariates, trim = c(0.02, 0.98)) {
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]
    ps <- fit_propensity_score(data, treat_var, covariates, trim = trim)

    w1 <- A / ps
    w0 <- (1 - A) / (1 - ps)
    p1 <- sum(w1 * Y) / sum(w1)
    p0 <- sum(w0 * Y) / sum(w0)
    make_or(p1, p0)
}

estimate_att_ipw_or <- function(data, outcome_var, treat_var, covariates, trim = c(0.02, 0.98)) {
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]
    ps <- fit_propensity_score(data, treat_var, covariates, trim = trim)

    w0 <- (1 - A) * ps / (1 - ps)
    p1 <- mean(Y[A == 1])
    p0 <- sum(w0 * Y) / sum(w0)
    make_or(p1, p0)
}

estimate_ps_matching_or <- function(data,
                                    outcome_var,
                                    treat_var,
                                    covariates,
                                    estimand = c("ATE", "ATT"),
                                    trim = c(0.02, 0.98),
                                    M = 1,
                                    replace = TRUE) {
    estimand <- match.arg(estimand)
    ps <- fit_propensity_score(data, treat_var, covariates, trim = trim)

    matching_fit <- Matching::Match(
        Y = data[[outcome_var]],
        Tr = data[[treat_var]],
        X = ps,
        M = M,
        estimand = estimand,
        replace = replace
    )

    p1 <- weighted.mean(
        data[[outcome_var]][matching_fit$index.treated],
        matching_fit$weights
    )
    p0 <- weighted.mean(
        data[[outcome_var]][matching_fit$index.control],
        matching_fit$weights
    )
    make_or(p1, p0)
}


# Estimation ----------------------------------------------------------------------------------

comparison <- rbindlist(lapply(outcome_vars, function(outcome_var) {
    data.table(
        outcome = outcome_var,
        estimand = c("Raw OR", "ATE", "ATE", "ATE", "ATT", "ATT", "ATT"),
        estimator = c(
            "Naive odds ratio",
            "Regression adjustment",
            "Hajek IPW",
            "PS matching",
            "Regression adjustment",
            "Hajek IPW",
            "PS matching"
        ),
        odds_ratio = c(
            estimate_naive_or(rhc, outcome_var, treat_var),
            estimate_regression_ate_or(rhc, outcome_var, treat_var, covariates),
            estimate_ate_ipw_or(rhc, outcome_var, treat_var, covariates),
            estimate_ps_matching_or(rhc, outcome_var, treat_var, covariates, estimand = "ATE"),
            estimate_regression_att_or(rhc, outcome_var, treat_var, covariates),
            estimate_att_ipw_or(rhc, outcome_var, treat_var, covariates),
            estimate_ps_matching_or(rhc, outcome_var, treat_var, covariates, estimand = "ATT")
        )
    )
}))

comparison[, log_odds_ratio := log(odds_ratio)]
fwrite(comparison, output_table_path)


# Plot ----------------------------------------------------------------------------------------

comparison[, label := paste(estimand, estimator, sep = ": ")]
comparison[, label := factor(label, levels = rev(unique(label)))]

or_plot <- ggplot(comparison, aes(x = label, y = odds_ratio, color = estimand)) +
    geom_point(size = 2.6) +
    geom_hline(yintercept = 1, color = "grey45", linetype = "dashed") +
    scale_y_log10() +
    facet_wrap(~ outcome) +
    coord_flip() +
    labs(
        x = NULL,
        y = "Estimated odds ratio, log scale",
        color = "Estimand",
        title = "RHC odds-ratio estimates for mortality outcomes"
    ) +
    theme_minimal()

ggsave(
    filename = output_plot_path,
    plot = or_plot,
    width = 8,
    height = 5,
    dpi = 300
)

rm(list = ls()); invisible(gc())
cat("...done\n")
