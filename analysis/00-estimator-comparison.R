#' README:
#' -------
#' - author: Andy Shin, Jian Kang, Liang-Cheng Chen
#' - date: 2026-05-28
#'
#' Desc:
#' -------
#' This file compares naive, regression-adjusted, IPW, matching, and AIPW
#' estimates for the RHC treatment effect on 30-day mortality.
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
source("estimators/ate/AIPW.R")
source("estimators/att/ipw.R")
source("estimators/att/matching.R")
source("estimators/att/AIPW.R")
source("estimators/utils.R")


# Main ----------------------------------------------------------------------------------------

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"
output_table_path <- "results/tables/estimator-comparison-death30d.csv"
output_plot_path <- "results/figures/estimator-comparison-death30d.png"

bootstrap_B <- 200
bootstrap_seed <- 2026
conf_level <- 0.95

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

estimate_naive_se <- function(data, outcome_var, treat_var) {
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]
    Y1 <- Y[A == 1]
    Y0 <- Y[A == 0]
    sqrt(var(Y1) / length(Y1) + var(Y0) / length(Y0))
}

estimate_regression_ate <- function(data, outcome_var, treat_var, covariates) {
    fit <- lm(
        make_outcome_formula(outcome_var, treat_var, covariates),
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
        make_outcome_formula(outcome_var, treat_var, covariates),
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

extract_matching_se <- function(matching_fit) {
    if (!is.null(matching_fit$se)) {
        return(as.numeric(matching_fit$se[[1]]))
    } else {
        return(NA)
    }
}

estimate_matching_row <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    estimand
) {
    if (estimand == "ATE") {
        matching_result <- estimate_ate_ps_matching(
            data,
            outcome_var,
            treat_var,
            covariates,
            return_details = TRUE
        )
    } else if (estimand == "ATT") {
        matching_result <- estimate_att_ps_matching(
            data,
            outcome_var,
            treat_var,
            covariates,
            return_details = TRUE
        )
    } else {
        stop("estimand must be ATE or ATT.")
    }

    data.table(
        Estimand = estimand,
        Method = "PS matching",
        Estimate = as.numeric(matching_result$estimate),
        SE = extract_matching_se(matching_result$fit)
    )
}

estimate_all_methods <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    include_matching = TRUE
) {
    comparison <- data.table(
        Estimand = c("Crude", rep("ATE", 3), rep("ATT", 3)),
        Method = c(
            "Naive comparison",
            rep(
                c(
                    "Regression adjustment",
                    "Hajek IPW",
                    "AIPW"
                ),
                2
            )
        ),
        Estimate = c(
            estimate_naive(data, outcome_var, treat_var),
            estimate_regression_ate(data, outcome_var, treat_var, covariates),
            estimate_ate_hajek_ipw(data, outcome_var, treat_var, covariates),
            estimate_ate_aipw(data, outcome_var, treat_var, covariates),
            estimate_regression_att(data, outcome_var, treat_var, covariates),
            estimate_att_hajek_ipw(data, outcome_var, treat_var, covariates),
            estimate_att_aipw(data, outcome_var, treat_var, covariates)
        )
    )

    if (!include_matching) {
        return(comparison)
    }

    matching_rows <- rbindlist(list(
        estimate_matching_row(data, outcome_var, treat_var, covariates, "ATE"),
        estimate_matching_row(data, outcome_var, treat_var, covariates, "ATT")
    ))

    rbindlist(list(comparison, matching_rows), fill = TRUE)
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
        Estimate = c(
            estimate_regression_ate(data, outcome_var, treat_var, covariates),
            estimate_ate_hajek_ipw(data, outcome_var, treat_var, covariates),
            estimate_ate_aipw(data, outcome_var, treat_var, covariates),
            estimate_regression_att(data, outcome_var, treat_var, covariates),
            estimate_att_hajek_ipw(data, outcome_var, treat_var, covariates),
            estimate_att_aipw(data, outcome_var, treat_var, covariates)
        )
    )
}


# Estimation ----------------------------------------------------------------------------------

comparison <- estimate_all_methods(
    rhc,
    outcome_var = outcome_var,
    treat_var = treat_var,
    covariates = covariates
)

set.seed(bootstrap_seed)
bootstrap_estimates <- replicate(
    bootstrap_B,
    {
        idx <- sample(nrow(rhc), replace = TRUE)
        bs_df <- rhc[idx, ]
        estimate_bootstrap_methods(
            bs_df,
            outcome_var = outcome_var,
            treat_var = treat_var,
            covariates = covariates
        )$Estimate
    }
)

bootstrap_methods <- estimate_bootstrap_methods(
    rhc,
    outcome_var = outcome_var,
    treat_var = treat_var,
    covariates = covariates
)[, .(Estimand, Method)]
bootstrap_methods[, SE := apply(bootstrap_estimates, 1, sd)]

comparison <- merge(
    comparison,
    bootstrap_methods,
    by = c("Estimand", "Method"),
    all.x = TRUE
)
comparison[, SE := fcoalesce(SE.x, SE.y)]
comparison[, c("SE.x", "SE.y") := NULL]
comparison[
    Method == "Naive comparison",
    SE := estimate_naive_se(
        rhc,
        outcome_var,
        treat_var
    )
]

z_value <- qnorm(1 - (1 - conf_level) / 2)
comparison[, ci_low := Estimate - z_value * SE]
comparison[, ci_high := Estimate + z_value * SE]
comparison[, `95%CI` := sprintf("[%.4f, %.4f]", ci_low, ci_high)]

comparison[, Estimand := factor(Estimand, levels = c("Crude", "ATE", "ATT"))]
comparison[,
    Method := factor(
        Method,
        levels = c(
            "Naive comparison",
            "Regression adjustment",
            "Hajek IPW",
            "PS matching",
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


# Plot ----------------------------------------------------------------------------------------

comparison[, label := paste(Estimand, Method, sep = ": ")]
comparison[, label := factor(label, levels = rev(label))]

comparison_plot <- ggplot(
    comparison,
    aes(x = label, y = Estimate, color = Estimand)
) +
    geom_errorbar(
        aes(ymin = ci_low, ymax = ci_high),
        width = 0.18,
        linewidth = 0.7
    ) +
    geom_point(size = 2.8) +
    geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
    coord_flip() +
    labs(
        x = NULL,
        y = "Estimated risk difference",
        color = "Estimand",
        title = "RHC effect estimates for 30-day mortality",
        subtitle = sprintf(
            "Naive and matching use analytic SEs; other bars use %s bootstrap resamples",
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
    plot = comparison_plot,
    width = 7,
    height = 4.5,
    dpi = 300
)

rm(list = ls())
invisible(gc())
cat("...done\n")
