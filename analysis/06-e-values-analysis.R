#' README:
#' -------
#' - author: Liang-Cheng Chen, Andy Shin
#' - date: 2026-06-10
#'
#' Desc:
#' -------
#' This file computes ATE E-values for hidden-confounding sensitivity analysis
#' on the risk-ratio scale. Following Lecture 7 and Lab 7, it reports E-values
#' for the point estimate and for the 95% CI bound closest to the null, then
#' benchmarks those E-values against observed covariate strengths.
#'
#' Input
#' -----
#' - data/processed/rhc-processed-data.csv
#' - data/processed/rhc-var-sets.rds
#'
#' Output
#' ------
#' - results/tables/evalue-sensitivity-ate-death30d.csv
#' - results/tables/evalue-benchmark-basic-ate-death30d.csv
#' - results/figures/evalue-sensitivity-death30d.png

.FILE_NAME <- "06-e-values-analysis.R"
cat(sprintf("...Running %s ...", .FILE_NAME))


# Packages ------------------------------------------------------------------------------------
library(data.table)
library(readr)
library(ggplot2)
library(EValue)

source("estimators/utils.R")


# Setups --------------------------------------------------------------------------------------

input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"

output_ate_table_path <- "results/tables/evalue-sensitivity-ate-death30d.csv"
output_ate_benchmark_path <- "results/tables/evalue-benchmark-basic-ate-death30d.csv"
output_plot_path <- "results/figures/evalue-sensitivity-death30d.png"

bootstrap_B <- 100
bootstrap_seed <- 2026
conf_level <- 0.95


# Main ----------------------------------------------------------------------------------------

rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)

outcome_var <- "death30d"
treat_var <- var_sets$treatment
covariates <- var_sets$adjustment_sets$basic


# Risk-Ratio Estimators -----------------------------------------------------------------------

estimate_regression_risks <- function(data) {
    outcome_fit <- glm(
        make_outcome_formula(outcome_var, treat_var, covariates),
        data = data,
        family = binomial("logit")
    )
    
    data1 <- copy(data)
    data0 <- copy(data)
    data1[, (treat_var) := 1]
    data0[, (treat_var) := 0]
    
    p1 <- mean(predict(outcome_fit, newdata = data1, type = "response"))
    p0 <- mean(predict(outcome_fit, newdata = data0, type = "response"))
    
    res <- data.table(
        estimand = "ATE",
        estimator = "Regression adjustment",
        risk_rhc = p1,
        risk_no_rhc = p0,
        risk_difference = p1 - p0,
        risk_ratio = p1 / p0
    )
    
    return(res)
}

estimate_aipw_risks <- function(data, trim = c(0.02, 0.98)) {
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]
    
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
    
    p1 <- mean(mu1 + A * (Y - mu1) / ps)
    p0 <- mean(mu0 + (1 - A) * (Y - mu0) / (1 - ps))
    
    res <- data.table(
        estimand = "ATE",
        estimator = "AIPW",
        risk_rhc = p1,
        risk_no_rhc = p0,
        risk_difference = p1 - p0,
        risk_ratio = p1 / p0
    )
    
    return(res)
}


# Benchmark Helpers ---------------------------------------------------------------------------

# This helper is kept because the Lecture 7/Lab 7 benchmark table needs every
# covariate on a comparable binary scale before RR_AW and RR_AY are calculated.
make_binary <- function(x) {
    if (length(unique(stats::na.omit(x))) <= 2) {
        res <- as.integer(x == max(x, na.rm = TRUE))
    } else {
        res <- as.integer(x > median(x, na.rm = TRUE))
    }
    
    return(res)
}

# This helper keeps the benchmark RR calculation readable and always reports the
# association strength on the >= 1 scale used by the E-value interpretation.
safe_rr <- function(num, den) {
    if (is.na(num) || is.na(den) || num == 0 || den == 0) {
        return(NA_real_)
    }
    
    rr <- num / den
    res <- max(rr, 1 / rr)
    
    return(res)
}

# This helper mirrors the Lab 7 bench_one() block. It is abstracted because the
# exact same RR_AW, RR_AY, weaker-leg calculation is repeated for every covariate.
bench_one <- function(W, A, Y) {
    W <- as.integer(W)
    
    if (length(unique(W[!is.na(W)])) < 2) {
        res <- c(RR_AW = NA, RR_AY = NA, weaker_leg = NA)
        return(res)
    }
    
    RR_AW <- safe_rr(
        mean(W[A == 1], na.rm = TRUE),
        mean(W[A == 0], na.rm = TRUE)
    )
    RR_AY <- safe_rr(
        mean(Y[W == 1], na.rm = TRUE),
        mean(Y[W == 0], na.rm = TRUE)
    )
    
    res <- c(
        RR_AW = RR_AW,
        RR_AY = RR_AY,
        weaker_leg = min(RR_AW, RR_AY, na.rm = TRUE)
    )
    
    return(res)
}


# 1. Point Estimates --------------------------------------------------------------------------

results <- rbindlist(list(
    estimate_regression_risks(rhc),
    estimate_aipw_risks(rhc)
))


# 2. Bootstrap CI For Risk Ratios -------------------------------------------------------------

set.seed(bootstrap_seed)
bs_results <- lapply(1:bootstrap_B, function(b) {
    idx <- sample(nrow(rhc), replace = TRUE)
    bs_data <- rhc[idx, ]
    
    bs_estimates <- rbindlist(list(
        estimate_regression_risks(bs_data),
        estimate_aipw_risks(bs_data)
    ))
    bs_estimates[, bootstrap_id := b]
    
    return(bs_estimates)
}) |> rbindlist()

alpha <- 1 - conf_level
bs_summary <- bs_results[
    ,
    .(
        rr_ci_low = quantile(risk_ratio, alpha / 2, na.rm = TRUE),
        rr_ci_high = quantile(risk_ratio, 1 - alpha / 2, na.rm = TRUE),
        rr_se = sd(risk_ratio, na.rm = TRUE)
    ),
    by = .(estimand, estimator)
]

results <- merge(
    results,
    bs_summary,
    by = c("estimand", "estimator"),
    all.x = TRUE
)


# 3. E-values ---------------------------------------------------------------------------------

evalue_results <- lapply(1:nrow(results), function(i) {
    ev <- EValue::evalues.RR(
        est = results$risk_ratio[i],
        lo = results$rr_ci_low[i],
        hi = results$rr_ci_high[i]
    )
    
    if (results$rr_ci_low[i] <= 1 && results$rr_ci_high[i] >= 1) {
        evalue_ci_limit <- 1
    } else if (results$risk_ratio[i] > 1) {
        evalue_ci_limit <- as.numeric(ev["E-values", "lower"])
    } else {
        evalue_ci_limit <- as.numeric(ev["E-values", "upper"])
    }
    
    res <- data.table(
        evalue_point = as.numeric(ev["E-values", "point"]),
        evalue_ci_limit = evalue_ci_limit
    )
    
    return(res)
}) |> rbindlist()

results[, `:=`(
    evalue_point = evalue_results$evalue_point,
    evalue_ci_limit = evalue_results$evalue_ci_limit,
    n_bootstrap = bootstrap_B
)]
results[, `95%CI` := sprintf("[%.4f, %.4f]", rr_ci_low, rr_ci_high)]


# 4. Observed-Covariate Benchmark -------------------------------------------------------------

A <- rhc[[treat_var]]
Y <- rhc[[outcome_var]]

benchmarks <- t(sapply(covariates, function(v) {
    W <- make_binary(rhc[[v]])
    res <- bench_one(W, A, Y)
    return(res)
}))

benchmarks <- as.data.frame(benchmarks)
benchmarks$covariate <- rownames(benchmarks)
benchmarks <- as.data.table(benchmarks)
benchmarks <- benchmarks[order(-weaker_leg), .(
    covariate,
    RR_AW,
    RR_AY,
    weaker_leg
)]

# strongest_benchmark <- benchmarks[1, weaker_leg]
# results[, strongest_observed_weaker_leg := strongest_benchmark]
# results[, margin_vs_strongest_observed := evalue_point - strongest_benchmark]


# 5. Output Tables ----------------------------------------------------------------------------

setorder(results, estimator)

result_table <- results[
    ,
    .(
        estimand,
        estimator,
        risk_rhc,
        risk_no_rhc,
        risk_difference,
        risk_ratio,
        rr_se,
        `95%CI`,
        evalue_point,
        evalue_ci_limit,
        strongest_observed_weaker_leg,
        margin_vs_strongest_observed,
        n_bootstrap
    )
]

fwrite(result_table, output_ate_table_path)
fwrite(benchmarks, output_ate_benchmark_path)


# 6. Plot -------------------------------------------------------------------------------------

plot_data <- melt(
    results,
    id.vars = c("estimand", "estimator"),
    measure.vars = c("evalue_point", "evalue_ci_limit"),
    variable.name = "quantity",
    value.name = "evalue"
)
plot_data[,
    quantity := fifelse(
        quantity == "evalue_point",
        "Point estimate",
        "CI bound"
    )
]
plot_data[, label := factor(estimator, levels = rev(unique(estimator)))]

plot <- ggplot(
    plot_data,
    aes(x = label, y = evalue, fill = quantity)
) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62) +
    geom_hline(yintercept = 1, color = "grey45") +
    geom_hline(
        yintercept = strongest_benchmark,
        color = "grey35",
        linetype = "dashed"
    ) +
    coord_flip() +
    labs(
        x = NULL,
        y = "E-value",
        fill = NULL,
        title = "ATE E-value sensitivity analysis for hidden confounding",
        subtitle = "Dashed line: strongest observed basic covariate weaker-leg RR"
    ) +
    theme_minimal()

ggsave(
    filename = output_plot_path,
    plot = plot,
    width = 8,
    height = 4.5,
    dpi = 300
)

rm(list = ls())
invisible(gc())
cat("...done\n")
