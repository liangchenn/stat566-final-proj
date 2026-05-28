#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-05-28
#'
#' Desc:
#' -------
#' This file checks whether ATE and ATT estimates are sensitive to the chosen
#' adjustment set.
#'
#' Input
#' -----
#' - data/processed/rhc-processed-data.csv
#' - data/processed/rhc-var-sets.rds
#'
#' Output
#' ------
#' - results/tables/adjustment-set-sensitivity-death30d.csv
#' - results/figures/adjustment-set-sensitivity-death30d.png

.FILE_NAME <- "04-adjustment-set-sensitivity.R"
cat(sprintf("...Running %s ...", .FILE_NAME))


# Packages ------------------------------------------------------------------------------------
library(data.table)
library(readr)
library(ggplot2)

source("estimators/ate/ipw.R")
source("estimators/ate/matching.R")
source("estimators/att/ipw.R")
source("estimators/att/matching.R")


# Main ----------------------------------------------------------------------------------------

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"
output_table_path <- "results/tables/adjustment-set-sensitivity-death30d.csv"
output_plot_path <- "results/figures/adjustment-set-sensitivity-death30d.png"

rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)

outcome_var <- "death30d"
treat_var <- var_sets$treatment

adjustment_sets <- var_sets$adjustment_sets[c("minimal", "basic", "detailed")]

results <- rbindlist(lapply(names(adjustment_sets), function(set_name) {
    covariates <- adjustment_sets[[set_name]]

    data.table(
        adjustment_set = set_name,
        estimand = c("ATE", "ATE", "ATT", "ATT"),
        estimator = c("Hajek IPW", "PS matching", "Hajek IPW", "PS matching"),
        estimate = c(
            estimate_ate_hajek_ipw(rhc, outcome_var, treat_var, covariates),
            estimate_ate_ps_matching(rhc, outcome_var, treat_var, covariates),
            estimate_att_hajek_ipw(rhc, outcome_var, treat_var, covariates),
            estimate_att_ps_matching(rhc, outcome_var, treat_var, covariates)
        ),
        n_covariates = length(covariates)
    )
}))

fwrite(results, output_table_path)

set_plot <- ggplot(
    results,
    aes(x = adjustment_set, y = estimate, color = estimator, group = estimator)
) +
    geom_point(size = 2.4) +
    geom_line() +
    geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
    facet_wrap(~ estimand) +
    labs(
        x = "Adjustment set",
        y = "Estimated risk difference",
        color = "Estimator",
        title = "Sensitivity to adjustment-set choice"
    ) +
    theme_minimal()

ggsave(
    filename = output_plot_path,
    plot = set_plot,
    width = 7,
    height = 4.5,
    dpi = 300
)

rm(list = ls()); invisible(gc())
cat("...done\n")
