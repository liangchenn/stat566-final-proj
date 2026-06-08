#' README:
#' -------
#' - author: Andy Shin, Liang-Cheng Chen
#' - date: 2026-05-28
#'
#' Desc:
#' -------
#' This file checks covariate balance before and after propensity-score
#' weighting for the basic adjustment set.
#'
#' Input
#' -----
#' - data/processed/rhc-processed-data.csv
#' - data/processed/rhc-var-sets.rds
#'
#' Output
#' ------
#' - results/tables/balance-diagnostics-basic.csv
#' - results/figures/balance-diagnostics-basic.png
#' - results/figures/propensity-score-overlap-basic.png

.FILE_NAME <- "01-balance-diagnostics.R"
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
output_table_path <- "results/tables/balance-diagnostics-basic.csv"
output_love_plot_path <- "results/figures/balance-diagnostics-basic.png"
output_overlap_plot_path <- "results/figures/propensity-score-overlap-basic.png"

rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)

treat_var <- var_sets$treatment
covariates <- var_sets$adjustment_sets$basic

ps <- fit_propensity_score(
    data = rhc,
    treat_var = treat_var,
    covariates = covariates
)

A <- rhc[[treat_var]]
rhc[, ps := ps]
rhc[, ate_weight := A / ps + (1 - A) / (1 - ps)]
rhc[, att_weight := A + (1 - A) * ps / (1 - ps)]


# Balance helpers -----------------------------------------------------------------------------

weighted_mean <- function(x, w) {
    sum(w * x, na.rm = TRUE) / sum(w[!is.na(x)])
}

weighted_var <- function(x, w) {
    m <- weighted_mean(x, w)
    sum(w * (x - m)^2, na.rm = TRUE) / sum(w[!is.na(x)])
}

weighted_smd <- function(data, var, treat_var, weight_var = NULL) {
    A <- data[[treat_var]]
    x <- data[[var]]

    if (is.null(weight_var)) {
        m1 <- mean(x[A == 1], na.rm = TRUE)
        m0 <- mean(x[A == 0], na.rm = TRUE)
        v1 <- stats::var(x[A == 1], na.rm = TRUE)
        v0 <- stats::var(x[A == 0], na.rm = TRUE)
    } else {
        w <- data[[weight_var]]
        m1 <- weighted_mean(x[A == 1], w[A == 1])
        m0 <- weighted_mean(x[A == 0], w[A == 0])
        v1 <- weighted_var(x[A == 1], w[A == 1])
        v0 <- weighted_var(x[A == 0], w[A == 0])
    }

    pooled_sd <- sqrt((v1 + v0) / 2)
    if (is.na(pooled_sd) || pooled_sd == 0) {
        return(0)
    }
    abs((m1 - m0) / pooled_sd)
}


# Diagnostics ---------------------------------------------------------------------------------

balance <- rbindlist(lapply(covariates, function(var) {
    data.table(
        covariate = var,
        raw = weighted_smd(rhc, var, treat_var),
        ate_ipw = weighted_smd(rhc, var, treat_var, "ate_weight"),
        att_ipw = weighted_smd(rhc, var, treat_var, "att_weight")
    )
}))

fwrite(balance[order(-raw)], output_table_path)

balance_long <- melt(
    balance,
    id.vars = "covariate",
    variable.name = "sample",
    value.name = "abs_smd"
)
balance_long[,
    covariate := factor(covariate, levels = balance[order(raw), covariate])
]

love_plot <- ggplot(
    balance_long,
    aes(x = covariate, y = abs_smd, color = sample)
) +
    geom_point(size = 2.2) +
    geom_hline(yintercept = 0.1, color = "grey45", linetype = "dashed") +
    coord_flip() +
    labs(
        x = NULL,
        y = "Absolute SMD",
        color = "Sample",
        title = "Covariate balance before and after IPW"
    ) +
    theme_minimal()

ggsave(
    filename = output_love_plot_path,
    plot = love_plot,
    width = 7,
    height = 5,
    dpi = 300
)

overlap_plot <- ggplot(rhc, aes(x = ps, fill = factor(get(treat_var)))) +
    geom_density(alpha = 0.45) +
    labs(
        x = "Estimated propensity score",
        y = "Density",
        fill = "RHC",
        title = "Propensity score overlap by treatment group"
    ) +
    theme_minimal()

ggsave(
    filename = output_overlap_plot_path,
    plot = overlap_plot,
    width = 7,
    height = 4.5,
    dpi = 300
)

rm(list = ls())
invisible(gc())
cat("...done\n")
