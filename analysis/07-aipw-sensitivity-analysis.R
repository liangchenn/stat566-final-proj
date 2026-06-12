#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-06-11
#'
#' Desc:
#' -------
#' This file implements a partial-R2 sensitivity analysis for the AIPW ATE.
#' Following the sensemakr-DML idea from Lecture 7, it first estimates a
#' cross-fitted AIPW/DML treatment effect and then asks how strong an
#' unmeasured confounder would need to be, on the sensemakr partial-R2 scale,
#' to attenuate the estimate to zero or to remove statistical significance.
#'
#' Input
#' -----
#' - data/processed/rhc-processed-data.csv
#' - data/processed/rhc-var-sets.rds
#'
#' Output
#' ------
#' - results/tables/aipw-sensitivity-summary-death30d.csv
#' - results/tables/aipw-sensitivity-benchmark-death30d.csv
#' - results/tables/aipw-sensitivity-scenarios-death30d.csv
#' - results/figures/aipw-sensitivity-contour-death30d.png

.FILE_NAME <- "07-aipw-sensitivity-analysis.R"
cat(sprintf("...Running %s ...", .FILE_NAME))


# Packages ------------------------------------------------------------------------------------
library(data.table)
library(readr)
library(ggplot2)
library(sensemakr)

source("estimators/utils.R")


# Setups --------------------------------------------------------------------------------------

input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"

output_summary_path <- "results/tables/aipw-sensitivity-summary-death30d.csv"
output_benchmark_path <- "results/tables/aipw-sensitivity-benchmark-death30d.csv"
output_scenario_path <- "results/tables/aipw-sensitivity-scenarios-death30d.csv"
output_plot_path <- "results/figures/aipw-sensitivity-contour-death30d.png"

crossfit_K <- 5
crossfit_seed <- 2026
trim <- c(0.02, 0.98)
conf_level <- 0.95
alpha <- 1 - conf_level


# Main ----------------------------------------------------------------------------------------

rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)

outcome_var <- "death30d"
treat_var <- var_sets$treatment
covariates <- var_sets$adjustment_sets$basic
required_vars <- unique(c(outcome_var, treat_var, covariates))

rhc <- rhc[, ..required_vars]
rhc <- na.omit(rhc)


# Helpers -------------------------------------------------------------------------------------

estimate_crossfit_aipw <- function(
    data,
    outcome_var,
    treat_var,
    covariates,
    K = 5,
    seed = 2026,
    trim = c(0.02, 0.98)
) {
    set.seed(seed)
    n <- nrow(data)
    folds <- sample(rep(1:K, length.out = n))

    scores <- rep(NA_real_, n)
    mu1_hat <- rep(NA_real_, n)
    mu0_hat <- rep(NA_real_, n)
    ps_hat <- rep(NA_real_, n)

    for (k in 1:K) {
        train_data <- as.data.frame(data[folds != k])
        test_data <- as.data.frame(data[folds == k])
        test_idx <- which(folds == k)

        outcome_fit <- glm(
            make_outcome_formula(outcome_var, treat_var, covariates),
            data = train_data,
            family = binomial("logit")
        )
        propensity_fit <- glm(
            make_formula(treat_var, covariates),
            data = train_data,
            family = binomial("logit")
        )

        test_treated <- test_data
        test_control <- test_data
        test_treated[[treat_var]] <- 1
        test_control[[treat_var]] <- 0

        mu1_k <- predict(outcome_fit, newdata = test_treated, type = "response")
        mu0_k <- predict(outcome_fit, newdata = test_control, type = "response")
        ps_k <- predict(propensity_fit, newdata = test_data, type = "response")
        ps_k <- pmin(pmax(ps_k, trim[1]), trim[2])

        A_k <- test_data[[treat_var]]
        Y_k <- test_data[[outcome_var]]

        scores[test_idx] <- (mu1_k - mu0_k) +
            A_k * (Y_k - mu1_k) / ps_k -
            (1 - A_k) * (Y_k - mu0_k) / (1 - ps_k)
        mu1_hat[test_idx] <- mu1_k
        mu0_hat[test_idx] <- mu0_k
        ps_hat[test_idx] <- ps_k
    }

    tau_hat <- mean(scores)
    se_hat <- sd(scores) / sqrt(n)
    z_value <- qnorm(1 - alpha / 2)

    result <- list(
        estimate = tau_hat,
        se = se_hat,
        ci_low = tau_hat - z_value * se_hat,
        ci_high = tau_hat + z_value * se_hat,
        scores = scores,
        mu1_hat = mu1_hat,
        mu0_hat = mu0_hat,
        ps_hat = ps_hat,
        riesz_score = data[[treat_var]] / ps_hat -
            (1 - data[[treat_var]]) / (1 - ps_hat)
    )

    return(result)
}

partial_r2_increment <- function(target, benchmark, adjustment_data) {
    dat <- as.data.frame(adjustment_data)
    dat$target <- target
    dat$benchmark <- benchmark

    if (length(unique(stats::na.omit(dat$benchmark))) < 2) {
        return(NA_real_)
    }

    reduced_fit <- lm(target ~ ., data = dat[, setdiff(names(dat), "benchmark")])
    full_fit <- lm(target ~ ., data = dat)

    rss_reduced <- sum(residuals(reduced_fit)^2)
    rss_full <- sum(residuals(full_fit)^2)

    if (rss_reduced <= 0) {
        return(NA_real_)
    }

    r2 <- 1 - rss_full / rss_reduced
    r2 <- pmin(pmax(r2, 0), 0.99)

    return(r2)
}

make_benchmark_table <- function(data, dml_fit, outcome_var, treat_var, covariates) {
    benchmarks <- rbindlist(lapply(covariates, function(benchmark_var) {
        adjustment_vars <- setdiff(covariates, benchmark_var)
        adjustment_data <- data[, ..adjustment_vars]

        r2_y <- partial_r2_increment(
            target = data[[outcome_var]],
            benchmark = data[[benchmark_var]],
            adjustment_data = cbind(data[, treat_var, with = FALSE], adjustment_data)
        )
        r2_d <- partial_r2_increment(
            target = dml_fit$riesz_score,
            benchmark = data[[benchmark_var]],
            adjustment_data = adjustment_data
        )

        data.table(
            covariate = benchmark_var,
            r2_riesz = r2_d,
            r2_outcome = r2_y,
            weaker_leg = pmin(r2_d, r2_y, na.rm = TRUE)
        )
    }))

    setorder(benchmarks, -weaker_leg)
    return(benchmarks)
}

make_sensitivity_scenarios <- function(summary_row, benchmarks, multipliers = c(1, 2, 3)) {
    scenarios <- rbindlist(lapply(multipliers, function(k) {
        res <- copy(benchmarks)
        res[, `:=`(
            multiplier = k,
            bound_label = sprintf("%sx %s", k, covariate),
            r2dz.x = pmin(k * r2_riesz, 0.99),
            r2yz.dx = pmin(k * r2_outcome, 0.99)
        )]
        return(res)
    }))

    scenarios[, `:=`(
        adjusted_estimate = sensemakr::adjusted_estimate(
            estimate = summary_row$estimate,
            se = summary_row$SE,
            dof = summary_row$dof,
            r2dz.x = r2dz.x,
            r2yz.dx = r2yz.dx
        ),
        adjusted_se = sensemakr::adjusted_se(
            se = summary_row$SE,
            dof = summary_row$dof,
            r2dz.x = r2dz.x,
            r2yz.dx = r2yz.dx
        ),
        adjusted_t = sensemakr::adjusted_t(
            estimate = summary_row$estimate,
            se = summary_row$SE,
            dof = summary_row$dof,
            r2dz.x = r2dz.x,
            r2yz.dx = r2yz.dx
        )
    )]
    scenarios[, adjusted_p := 2 * pt(-abs(adjusted_t), df = summary_row$dof)]
    scenarios[, adjusted_significant := adjusted_p < alpha]

    return(scenarios[])
}


# Estimation ----------------------------------------------------------------------------------

dml_fit <- estimate_crossfit_aipw(
    data = rhc,
    outcome_var = outcome_var,
    treat_var = treat_var,
    covariates = covariates,
    K = crossfit_K,
    seed = crossfit_seed,
    trim = trim
)

dof <- nrow(rhc) - length(covariates) - 2
sense_stats <- sensemakr::sensitivity_stats(
    estimate = dml_fit$estimate,
    se = dml_fit$se,
    dof = dof,
    treatment = treat_var,
    alpha = alpha
)

summary_results <- data.table(
    estimand = "ATE",
    estimator = "Cross-fitted AIPW",
    sensitivity_framework = "sensemakr-DML partial R2",
    outcome = outcome_var,
    treatment = treat_var,
    estimate = dml_fit$estimate,
    SE = dml_fit$se,
    ci_low = dml_fit$ci_low,
    ci_high = dml_fit$ci_high,
    `95%CI` = sprintf("[%.4f, %.4f]", dml_fit$ci_low, dml_fit$ci_high),
    partial_r2_treatment_outcome = sense_stats$r2yd.x,
    rv_point = as.numeric(sense_stats$rv_q),
    rv_significance = as.numeric(sense_stats$rv_qa),
    dof = dof,
    n = nrow(rhc),
    crossfit_K = crossfit_K,
    trim_lower = trim[1],
    trim_upper = trim[2]
)


# Benchmarks ----------------------------------------------------------------------------------

benchmarks <- make_benchmark_table(
    data = rhc,
    dml_fit = dml_fit,
    outcome_var = outcome_var,
    treat_var = treat_var,
    covariates = covariates
)

scenario_results <- make_sensitivity_scenarios(
    summary_row = summary_results[1],
    benchmarks = benchmarks
)

scenario_results <- scenario_results[
    ,
    .(
        bound_label,
        covariate,
        multiplier,
        r2dz.x,
        r2yz.dx,
        adjusted_estimate,
        adjusted_se,
        adjusted_t,
        adjusted_p,
        adjusted_significant
    )
]

fwrite(summary_results, output_summary_path)
fwrite(benchmarks, output_benchmark_path)
fwrite(scenario_results, output_scenario_path)


# Plot ----------------------------------------------------------------------------------------

axis_max <- max(
    summary_results$rv_point * 1.6,
    benchmarks[1:min(8, .N), max(c(r2_riesz, r2_outcome), na.rm = TRUE)] * 3,
    0.08,
    na.rm = TRUE
)
axis_max <- min(axis_max, 0.35)

grid <- CJ(
    r2dz.x = seq(0, axis_max, length.out = 100),
    r2yz.dx = seq(0, axis_max, length.out = 100)
)
grid[, adjusted_estimate := sensemakr::adjusted_estimate(
    estimate = summary_results$estimate,
    se = summary_results$SE,
    dof = summary_results$dof,
    r2dz.x = r2dz.x,
    r2yz.dx = r2yz.dx
)]

benchmark_plot_data <- copy(benchmarks[1:min(8, .N)])
benchmark_plot_data[, covariate := factor(covariate, levels = rev(covariate))]

plot <- ggplot(grid, aes(x = r2dz.x, y = r2yz.dx, z = adjusted_estimate)) +
    geom_raster(aes(fill = adjusted_estimate), interpolate = TRUE) +
    geom_contour(color = "white", linewidth = 0.3, bins = 8) +
    geom_contour(
        breaks = 0,
        color = "black",
        linewidth = 0.75
    ) +
    geom_point(
        data = benchmark_plot_data,
        aes(x = r2_riesz, y = r2_outcome),
        inherit.aes = FALSE,
        size = 2.2,
        color = "black"
    ) +
    geom_text(
        data = benchmark_plot_data,
        aes(x = r2_riesz, y = r2_outcome, label = covariate),
        inherit.aes = FALSE,
        hjust = -0.05,
        vjust = 0.4,
        size = 3.0
    ) +
    geom_abline(
        intercept = 0,
        slope = 1,
        color = "grey20",
        linetype = "dashed",
        linewidth = 0.35
    ) +
    coord_cartesian(xlim = c(0, axis_max), ylim = c(0, axis_max), expand = FALSE) +
    scale_fill_gradient2(
        low = "tomato",
        mid = "lightgrey",
        high = "steelblue",
        midpoint = 0,
        name = "Adjusted\nestimate"
    ) +
    labs(
        x = "Partial R2 with AIPW",
        y = "Partial R2 with outcome",
        title = "AIPW sensitivity analysis with sensemakr-DML"
    ) +
    theme_minimal()

ggsave(
    filename = output_plot_path,
    plot = plot,
    width = 8,
    height = 5.5,
    dpi = 300
)

rm(list = ls())
invisible(gc())
cat("...done\n")
