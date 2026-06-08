#' README:
#' -------
#' - author: Jian Kang, Liang-Cheng Chen
#' - date: 2026-06-01
#'
#' Desc:
#' -------
#' This file adds bootstrap uncertainty for the primary Hajek IPW risk-difference
#' estimate and computes E-values for the corresponding mortality risk ratio.
#'
#' Input
#' -----
#' - data/processed/rhc-processed-data.csv
#' - data/processed/rhc-var-sets.rds
#'
#' Output
#' ------
#' - results/tables/primary-ipw-estimate-death30d.csv
#' - results/tables/evalue-sensitivity-death30d.csv
#' - results/figures/evalue-sensitivity-death30d.png

.FILE_NAME <- "06-sensitivity-analysis.R"
cat(sprintf("...Running %s ...", .FILE_NAME))


# Packages ------------------------------------------------------------------------------------
library(data.table)
library(readr)
library(ggplot2)

source("estimators/utils.R")


# Main ----------------------------------------------------------------------------------------

set.seed(2026)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"
primary_table_path <- "results/tables/primary-ipw-estimate-death30d.csv"
evalue_table_path <- "results/tables/evalue-sensitivity-death30d.csv"
evalue_plot_path <- "results/figures/evalue-sensitivity-death30d.png"

rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)

outcome_var <- "death30d"
treat_var <- var_sets$treatment
covariates <- var_sets$adjustment_sets$basic
n_bootstrap <- 200


# Estimation helpers --------------------------------------------------------------------------

effective_sample_size <- function(w) {
    sum(w)^2 / sum(w^2)
}

make_evalue <- function(rr) {
    rr <- as.numeric(rr)
    rr_star <- ifelse(rr < 1, 1 / rr, rr)
    ifelse(rr_star <= 1, 1, rr_star + sqrt(rr_star * (rr_star - 1)))
}

estimate_ipw_contrast <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    estimand = c("ATE", "ATT"),
    trim = c(0.02, 0.98)
) {
    estimand <- match.arg(estimand)
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]
    ps <- fit_propensity_score(data, treat_var, covariates, trim = trim)

    if (estimand == "ATE") {
        w1 <- A / ps
        w0 <- (1 - A) / (1 - ps)
        ess <- effective_sample_size(w1 + w0)
    } else {
        w1 <- A
        w0 <- (1 - A) * ps / (1 - ps)
        ess <- effective_sample_size(w1 + w0)
    }

    p1 <- sum(w1 * Y) / sum(w1)
    p0 <- sum(w0 * Y) / sum(w0)

    data.table(
        estimand = estimand,
        risk_rhc = p1,
        risk_no_rhc = p0,
        risk_difference = p1 - p0,
        risk_ratio = p1 / p0,
        odds_ratio = (p1 / (1 - p1)) / (p0 / (1 - p0)),
        ess = ess,
        n = nrow(data)
    )
}

bootstrap_ipw <- function(data, estimand) {
    point <- estimate_ipw_contrast(
        data,
        outcome_var,
        treat_var,
        covariates,
        estimand
    )

    boot <- rbindlist(lapply(seq_len(n_bootstrap), function(i) {
        idx <- sample.int(nrow(data), nrow(data), replace = TRUE)
        estimate_ipw_contrast(
            data[idx],
            outcome_var,
            treat_var,
            covariates,
            estimand
        )
    }))

    alpha <- 0.05
    point[, `:=`(
        rd_se = stats::sd(boot$risk_difference, na.rm = TRUE),
        rd_ci_lower = stats::quantile(
            boot$risk_difference,
            alpha / 2,
            na.rm = TRUE
        ),
        rd_ci_upper = stats::quantile(
            boot$risk_difference,
            1 - alpha / 2,
            na.rm = TRUE
        ),
        rr_se = stats::sd(boot$risk_ratio, na.rm = TRUE),
        rr_ci_lower = stats::quantile(boot$risk_ratio, alpha / 2, na.rm = TRUE),
        rr_ci_upper = stats::quantile(
            boot$risk_ratio,
            1 - alpha / 2,
            na.rm = TRUE
        ),
        n_bootstrap = n_bootstrap
    )]

    point[]
}


# Primary estimates and E-values ---------------------------------------------------------------

primary <- rbind(
    bootstrap_ipw(rhc, "ATE"),
    bootstrap_ipw(rhc, "ATT")
)
primary[, estimator := "Hajek IPW"]
setcolorder(
    primary,
    c(
        "estimand",
        "estimator",
        "risk_rhc",
        "risk_no_rhc",
        "risk_difference",
        "rd_se",
        "rd_ci_lower",
        "rd_ci_upper",
        "risk_ratio",
        "rr_se",
        "rr_ci_lower",
        "rr_ci_upper",
        "odds_ratio",
        "ess",
        "n",
        "n_bootstrap"
    )
)

evalues <- copy(primary)
evalues[, `:=`(
    evalue_point = make_evalue(risk_ratio),
    evalue_ci_limit = fifelse(
        rr_ci_lower <= 1 & rr_ci_upper >= 1,
        1,
        pmin(make_evalue(rr_ci_lower), make_evalue(rr_ci_upper))
    ),
    interpretation = paste0(
        "An unmeasured confounder would need risk-ratio associations of at least ",
        sprintf("%.2f", make_evalue(risk_ratio)),
        " with both RHC receipt and 30-day mortality, conditional on measured covariates, ",
        "to explain away the point estimate."
    )
)]
evalues <- evalues[, .(
    estimand,
    estimator,
    risk_ratio,
    rr_ci_lower,
    rr_ci_upper,
    evalue_point,
    evalue_ci_limit,
    interpretation
)]

fwrite(primary, primary_table_path)
fwrite(evalues, evalue_table_path)


# Plot ----------------------------------------------------------------------------------------

evalue_plot_data <- melt(
    evalues,
    id.vars = c("estimand", "estimator"),
    measure.vars = c("evalue_point", "evalue_ci_limit"),
    variable.name = "quantity",
    value.name = "evalue"
)
evalue_plot_data[,
    quantity := fifelse(
        quantity == "evalue_point",
        "Point estimate",
        "CI limit"
    )
]

evalue_plot <- ggplot(
    evalue_plot_data,
    aes(x = estimand, y = evalue, fill = quantity)
) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 1, color = "grey45", linetype = "dashed") +
    labs(
        x = "Estimand",
        y = "E-value",
        fill = NULL,
        title = "E-value sensitivity analysis for hidden confounding"
    ) +
    theme_minimal()

ggsave(
    filename = evalue_plot_path,
    plot = evalue_plot,
    width = 7,
    height = 4.5,
    dpi = 300
)

rm(list = ls())
invisible(gc())
cat("...done\n")
