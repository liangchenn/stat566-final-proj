#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-05-28
#'
#' Desc:
#' -------
#' This file checks sensitivity of Hajek IPW estimates to propensity-score
#' trimming thresholds.
#'
#' Input
#' -----
#' - data/processed/rhc-processed-data.csv
#' - data/processed/rhc-var-sets.rds
#'
#' Output
#' ------
#' - results/tables/trimming-sensitivity-death30d.csv
#' - results/figures/trimming-sensitivity-death30d.png

.FILE_NAME <- "03-trimming-sensitivity.R"
cat(sprintf("...Running %s ...", .FILE_NAME))


# Packages ------------------------------------------------------------------------------------
library(data.table)
library(readr)
library(ggplot2)

source("estimators/ate/ipw.R")
source("estimators/att/ipw.R")
source("estimators/utils.R")


# Main ----------------------------------------------------------------------------------------

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"
output_table_path <- "results/tables/trimming-sensitivity-death30d.csv"
output_plot_path <- "results/figures/trimming-sensitivity-death30d.png"

rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)

outcome_var <- "death30d"
treat_var <- var_sets$treatment
covariates <- var_sets$adjustment_sets$basic

trim_grid <- list(
    c(0.001, 0.999),
    c(0.010, 0.990),
    c(0.020, 0.980),
    c(0.050, 0.950),
    c(0.100, 0.900)
)

effective_sample_size <- function(w) {
    sum(w)^2 / sum(w^2)
}

results <- rbindlist(lapply(trim_grid, function(trim) {
    ps <- fit_propensity_score(
        data = rhc,
        treat_var = treat_var,
        covariates = covariates,
        trim = trim
    )
    A <- rhc[[treat_var]]
    ate_weight <- A / ps + (1 - A) / (1 - ps)
    att_weight <- A + (1 - A) * ps / (1 - ps)

    rbind(
        data.table(
            estimand = "ATE",
            estimator = "Hajek IPW",
            trim_lower = trim[1],
            trim_upper = trim[2],
            estimate = estimate_ate_hajek_ipw(
                rhc, outcome_var, treat_var, covariates, trim = trim
            ),
            ess = effective_sample_size(ate_weight)
        ),
        data.table(
            estimand = "ATT",
            estimator = "Hajek IPW",
            trim_lower = trim[1],
            trim_upper = trim[2],
            estimate = estimate_att_hajek_ipw(
                rhc, outcome_var, treat_var, covariates, trim = trim
            ),
            ess = effective_sample_size(att_weight)
        )
    )
}))

results[, trim := sprintf("[%.3f, %.3f]", trim_lower, trim_upper)]
fwrite(results, output_table_path)

trim_plot <- ggplot(results, aes(x = trim, y = estimate, color = estimand, group = estimand)) +
    geom_point(size = 2.4) +
    geom_line() +
    geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
    labs(
        x = "Propensity-score trim",
        y = "Estimated risk difference",
        color = "Estimand",
        title = "IPW sensitivity to propensity-score trimming"
    ) +
    theme_minimal()

ggsave(
    filename = output_plot_path,
    plot = trim_plot,
    width = 7,
    height = 4.5,
    dpi = 300
)

rm(list = ls()); invisible(gc())
cat("...done\n")
