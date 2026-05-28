#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-05-28
#'
#' Desc:
#' -------
#' This file compares naive, regression-adjusted, IPW, and matching estimates
#' for the RHC treatment effect on 30-day mortality.
#'
#' Input
#' -----
#' - data/processed/rhc-processed-data.csv
#' - data/processed/rhc-var-sets.rds
#'
#' Output
#' ------
#' - results/tables/estimator-comparison-death30d.csv
#' - results/figures/estimator-comparison-death30d.png

.FILE_NAME <- "00-estimator-comparison.R"
cat(sprintf("...Running %s ...", .FILE_NAME))


# Packages ------------------------------------------------------------------------------------
library(data.table)
library(readr)
library(ggplot2)

source("estimators/ate/ipw.R")
source("estimators/ate/matching.R")
source("estimators/att/ipw.R")
source("estimators/att/matching.R")
source("estimators/utils.R")


# Main ----------------------------------------------------------------------------------------

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"
output_table_path <- "results/tables/estimator-comparison-death30d.csv"
output_plot_path <- "results/figures/estimator-comparison-death30d.png"

rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)

outcome_var <- "death30d"
treat_var <- var_sets$treatment
covariates <- var_sets$adjustment_sets$basic


# Estimator helpers ---------------------------------------------------------------------------

estimate_naive <- function(data, outcome_var, treat_var) {
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]
    mean(Y[A == 1]) - mean(Y[A == 0])
}

estimate_regression_ate <- function(data, outcome_var, treat_var, covariates) {
    fit <- lm(
        make_interacted_outcome_formula(outcome_var, treat_var, covariates),
        data = data
    )
    data1 <- copy(data)
    data0 <- copy(data)
    data1[, (treat_var) := 1]
    data0[, (treat_var) := 0]
    mu1 <- predict(fit, newdata = data1)
    mu0 <- predict(fit, newdata = data0)
    mean(mu1 - mu0)
}

estimate_regression_att <- function(data, outcome_var, treat_var, covariates) {
    fit <- lm(
        make_interacted_outcome_formula(outcome_var, treat_var, covariates),
        data = data
    )
    treated_data <- data[data[[treat_var]] == 1, ]
    treated_data1 <- copy(treated_data)
    treated_data0 <- copy(treated_data)
    treated_data1[, (treat_var) := 1]
    treated_data0[, (treat_var) := 0]
    mu1 <- predict(fit, newdata = treated_data1)
    mu0 <- predict(fit, newdata = treated_data0)
    mean(mu1 - mu0)
}


# Estimation ----------------------------------------------------------------------------------

comparison <- data.table(
    estimand = c("Raw gap", "ATE", "ATE", "ATE", "ATT", "ATT", "ATT"),
    estimator = c(
        "Naive difference",
        "Regression adjustment",
        "Hajek IPW",
        "PS matching",
        "Regression adjustment",
        "Hajek IPW",
        "PS matching"
    ),
    estimate = c(
        estimate_naive(rhc, outcome_var, treat_var),
        estimate_regression_ate(rhc, outcome_var, treat_var, covariates),
        estimate_ate_hajek_ipw(rhc, outcome_var, treat_var, covariates),
        estimate_ate_ps_matching(rhc, outcome_var, treat_var, covariates),
        estimate_regression_att(rhc, outcome_var, treat_var, covariates),
        estimate_att_hajek_ipw(rhc, outcome_var, treat_var, covariates),
        estimate_att_ps_matching(rhc, outcome_var, treat_var, covariates)
    )
)

fwrite(comparison, output_table_path)


# Plot ----------------------------------------------------------------------------------------

comparison[, label := paste(estimand, estimator, sep = ": ")]
comparison[, label := factor(label, levels = rev(label))]

comparison_plot <- ggplot(comparison, aes(x = label, y = estimate, color = estimand)) +
    geom_point(size = 2.8) +
    geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
    coord_flip() +
    labs(
        x = NULL,
        y = "Estimated risk difference",
        color = "Estimand",
        title = "RHC effect estimates for 30-day mortality"
    ) +
    theme_minimal()

ggsave(
    filename = output_plot_path,
    plot = comparison_plot,
    width = 7,
    height = 4.5,
    dpi = 300
)

rm(list = ls()); invisible(gc())
cat("...done\n")
