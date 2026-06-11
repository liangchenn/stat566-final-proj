#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-06-05
#'
#' Desc:
#' -------
#' This file compares naive, regression-adjusted, IPW, matching, and AIPW
#' odds-ratio estimates for the RHC effect on 30-day mortality.
#'
#' Input
#' -----
#' - data/processed/rhc-processed-data.csv
#' - data/processed/rhc-var-sets.rds
#'
#' Output
#' ------
#' - results/tables/odds-ratio-estimator-comparison-death30d.csv
#' - results/figures/odds-ratio-estimator-comparison-death30d.png

.FILE_NAME <- "02-odds-ratio-estimator-comparison.R"
cat(sprintf("...Running %s ...", .FILE_NAME))


# Packages ------------------------------------------------------------------------------------
library(data.table)
library(readr)
library(ggplot2)

source("estimators/utils.R")


# Setups --------------------------------------------------------------------------------------

input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"

output_table_path <- "results/tables/odds-ratio-estimator-comparison-death30d.csv"
output_plot_path <- "results/figures/odds-ratio-estimator-comparison-death30d.png"
output_tex_path <- "results/tables/odds-ratio-comparison-death30d.tex"

bootstrap_B <- 200
bootstrap_seed <- 2026
conf_level <- 0.95


# Main ----------------------------------------------------------------------------------------

rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)

outcome_var <- "death30d"
treat_var <- var_sets$treatment
covariates <- var_sets$adjustment_sets$basic


# Odds-ratio helpers --------------------------------------------------------------------------

make_odd_ratio <- function(p1, p0) {
    # eps = 1e-6
    # p1 <- pmin(pmax(p1, eps), 1 - eps)
    # p0 <- pmin(pmax(p0, eps), 1 - eps)
    odd_ratio <- (p1 / (1 - p1)) / (p0 / (1 - p0))
}

estimate_naive_or <- function(data, outcome_var, treat_var) {
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]
    make_odd_ratio(mean(Y[A == 1]), mean(Y[A == 0]))
}

estimate_naive_log_or_se <- function(
    data,
    outcome_var,
    treat_var,
    correction = 0.5
) {
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]

    n11 <- sum(A == 1 & Y == 1)
    n10 <- sum(A == 1 & Y == 0)
    n01 <- sum(A == 0 & Y == 1)
    n00 <- sum(A == 0 & Y == 0)

    sqrt(
        1 /
            (n11 + correction) +
            1 / (n10 + correction) +
            1 / (n01 + correction) +
            1 / (n00 + correction)
    )
}

estimate_regression_ate_or <- function(
    data,
    outcome_var,
    treat_var,
    covariates
) {
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
    make_odd_ratio(p1, p0)
}

estimate_regression_att_or <- function(
    data,
    outcome_var,
    treat_var,
    covariates
) {
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
    make_odd_ratio(p1, p0)
}

estimate_ate_ipw_or <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    trim = c(0.02, 0.98)
) {
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]
    ps <- fit_propensity_score(data, treat_var, covariates, trim = trim)

    w1 <- A / ps
    w0 <- (1 - A) / (1 - ps)
    p1 <- sum(w1 * Y) / sum(w1)
    p0 <- sum(w0 * Y) / sum(w0)
    make_odd_ratio(p1, p0)
}

estimate_att_ipw_or <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    trim = c(0.02, 0.98)
) {
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]
    ps <- fit_propensity_score(data, treat_var, covariates, trim = trim)

    w0 <- (1 - A) * ps / (1 - ps)
    p1 <- mean(Y[A == 1])
    p0 <- sum(w0 * Y) / sum(w0)
    make_odd_ratio(p1, p0)
}

estimate_ate_aipw_or <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    trim = c(0.02, 0.98)
) {
    outcome_fit <- glm(
        make_outcome_formula(outcome_var, treat_var, covariates),
        data = data,
        family = binomial("logit")
    )
    ps <- fit_propensity_score(data, treat_var, covariates, trim = trim)

    data1 <- copy(data)
    data0 <- copy(data)
    data1[, (treat_var) := 1]
    data0[, (treat_var) := 0]

    mu1 <- predict(outcome_fit, newdata = data1, type = "response")
    mu0 <- predict(outcome_fit, newdata = data0, type = "response")
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]

    p1 <- mean(mu1 + A * (Y - mu1) / ps)
    p0 <- mean(mu0 + (1 - A) * (Y - mu0) / (1 - ps))
    make_odd_ratio(p1, p0)
}

estimate_att_aipw_or <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    trim = c(0.02, 0.98)
) {
    control_data <- data[data[[treat_var]] == 0, ]
    outcome_fit <- glm(
        make_formula(outcome_var, covariates),
        data = control_data,
        family = binomial("logit")
    )
    ps <- fit_propensity_score(data, treat_var, covariates, trim = trim)

    mu0 <- predict(outcome_fit, newdata = data, type = "response")
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]

    p1 <- mean(Y[A == 1])
    p0 <- mean(A * mu0 + (1 - A) * ps / (1 - ps) * (Y - mu0)) / mean(A)
    make_odd_ratio(p1, p0)
}

estimate_ps_matching_or <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    estimand = c("ATE", "ATT"),
    trim = c(0.02, 0.98),
    M = 1,
    replace = TRUE
) {
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
    make_odd_ratio(p1, p0)
}

estimate_ps_matching_log_or_se <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    estimand = c("ATE", "ATT"),
    trim = c(0.02, 0.98),
    M = 1,
    replace = TRUE,
    correction = 0.5
) {
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

    y1 <- data[[outcome_var]][matching_fit$index.treated]
    y0 <- data[[outcome_var]][matching_fit$index.control]
    weights <- matching_fit$weights

    n11 <- sum(weights * y1)
    n10 <- sum(weights * (1 - y1))
    n01 <- sum(weights * y0)
    n00 <- sum(weights * (1 - y0))

    sqrt(
        1 /
            (n11 + correction) +
            1 / (n10 + correction) +
            1 / (n01 + correction) +
            1 / (n00 + correction)
    )
}

estimate_all_methods <- function(data, outcome_var, treat_var, covariates) {
    data.table(
        Estimand = c("Naive", rep("ATE", 4), rep("ATT", 4)),
        Method = c(
            "Naive comparison",
            rep(
                c(
                    "Regression adjustment",
                    "Hajek IPW",
                    "PS matching",
                    "AIPW"
                ),
                2
            )
        ),
        Estimate = c(
            estimate_naive_or(data, outcome_var, treat_var),
            estimate_regression_ate_or(
                data,
                outcome_var,
                treat_var,
                covariates
            ),
            estimate_ate_ipw_or(data, outcome_var, treat_var, covariates),
            estimate_ps_matching_or(
                data,
                outcome_var,
                treat_var,
                covariates,
                estimand = "ATE"
            ),
            estimate_ate_aipw_or(data, outcome_var, treat_var, covariates),
            estimate_regression_att_or(
                data,
                outcome_var,
                treat_var,
                covariates
            ),
            estimate_att_ipw_or(data, outcome_var, treat_var, covariates),
            estimate_ps_matching_or(
                data,
                outcome_var,
                treat_var,
                covariates,
                estimand = "ATT"
            ),
            estimate_att_aipw_or(data, outcome_var, treat_var, covariates)
        )
    )
}

estimate_bootstrap_methods <- function(
    data,
    outcome_var,
    treat_var,
    covariates
) {
    data.table(
        Estimand = c(rep("ATE", 3), rep("ATT", 3)),
        Method = rep(
            c(
                "Regression adjustment",
                "Hajek IPW",
                "AIPW"
            ),
            2
        ),
        log_or = log(c(
            estimate_regression_ate_or(
                data,
                outcome_var,
                treat_var,
                covariates
            ),
            estimate_ate_ipw_or(data, outcome_var, treat_var, covariates),
            estimate_ate_aipw_or(data, outcome_var, treat_var, covariates),
            estimate_regression_att_or(
                data,
                outcome_var,
                treat_var,
                covariates
            ),
            estimate_att_ipw_or(data, outcome_var, treat_var, covariates),
            estimate_att_aipw_or(data, outcome_var, treat_var, covariates)
        ))
    )
}


# Estimation ----------------------------------------------------------------------------------

comparison <- estimate_all_methods(
    rhc,
    outcome_var = outcome_var,
    treat_var = treat_var,
    covariates = covariates
)
comparison[, log_or := log(Estimate)]

set.seed(bootstrap_seed)
bootstrap_log_or <- replicate(
    bootstrap_B,
    {
        idx <- sample(nrow(rhc), replace = TRUE)
        bs_df <- rhc[idx, ]
        estimate_bootstrap_methods(
            bs_df,
            outcome_var = outcome_var,
            treat_var = treat_var,
            covariates = covariates
        )$log_or
    }
)

bootstrap_methods <- estimate_bootstrap_methods(
    rhc,
    outcome_var = outcome_var,
    treat_var = treat_var,
    covariates = covariates
)[, .(Estimand, Method)]
bootstrap_methods[, log_or_se := apply(bootstrap_log_or, 1, sd)]

comparison <- merge(
    comparison,
    bootstrap_methods,
    by = c("Estimand", "Method"),
    all.x = TRUE
)
comparison[
    Method == "Naive comparison",
    log_or_se := estimate_naive_log_or_se(
        rhc,
        outcome_var,
        treat_var
    )
]
comparison[
    Estimand == "ATE" & Method == "PS matching",
    log_or_se := estimate_ps_matching_log_or_se(
        rhc,
        outcome_var,
        treat_var,
        covariates,
        estimand = "ATE"
    )
]
comparison[
    Estimand == "ATT" & Method == "PS matching",
    log_or_se := estimate_ps_matching_log_or_se(
        rhc,
        outcome_var,
        treat_var,
        covariates,
        estimand = "ATT"
    )
]

z_value <- qnorm(1 - (1 - conf_level) / 2)
comparison[, SE := Estimate * log_or_se]
comparison[, ci_low := exp(log_or - z_value * log_or_se)]
comparison[, ci_high := exp(log_or + z_value * log_or_se)]
comparison[, `95%CI` := sprintf("[%.4f, %.4f]", ci_low, ci_high)]

comparison[, Estimand := factor(Estimand, levels = c("Naive", "ATE", "ATT"))]
comparison[,
    Method := factor(
        Method,
        levels = c(
            "Naive comparison",
            "Regression adjustment",
            "PS matching",
            "Hajek IPW",
            "AIPW"
        )
    )
]
setorder(comparison, Estimand, Method)

comparison_table <- comparison[, .(
    Estimand,
    Method,
    Estimate,
    SE,
    `95%CI`
)]

fwrite(comparison_table, output_table_path)

# Post-Processing for tex table ---------------------------------------------------------------

TAB_LABEL <- "tab:or_estimates"

x_tab <- xtable(
    comparison_table,
    caption = "Odds Ratio Results by different causal estimators",
    label = TAB_LABEL,
    align = "rllccc",
    digits = 4
)
# add notes
note_text <- "\\addlinespace \n \\multicolumn{5}{l}{\\small Note: The SE for IPW, AIPW were obtained with bootstrapping with 200 times.} \\\\ \n"
# output
print(
    x_tab,
    include.rownames = FALSE,
    booktabs = TRUE,
    caption.placement = "top",
    add.to.row = list(pos = list(nrow(x_tab)), command = note_text), # <- 核心在這
    file = output_tex_path
)


# Plot ----------------------------------------------------------------------------------------

comparison[, label := paste(Estimand, Method, sep = ": ")]
comparison[, label := factor(label, levels = rev(label))]

or_plot <- ggplot(comparison, aes(x = label, y = Estimate, color = Estimand)) +
    geom_errorbar(
        aes(ymin = ci_low, ymax = ci_high),
        width = 0.18,
        linewidth = 0.7
    ) +
    geom_point(size = 2.6) +
    geom_hline(yintercept = 1, color = "grey45", linetype = "dashed") +
    scale_y_log10() +
    coord_flip() +
    labs(
        x = NULL,
        y = "Estimated odds ratio, log scale",
        color = "Estimand",
        title = "RHC odds-ratio estimates for 30-day mortality",
        subtitle = sprintf(
            "Naive and matching use analytic log-OR SEs; other bars use %s bootstrap resamples",
            bootstrap_B
        )
    ) +
    theme_minimal() +
    theme(
        panel.grid.major.y = element_blank(),
        legend.position = "bottom"
    )

ggsave(
    filename = output_plot_path,
    plot = or_plot,
    width = 8,
    height = 5,
    dpi = 300
)

rm(list = ls())
invisible(gc())
cat("...done\n")
