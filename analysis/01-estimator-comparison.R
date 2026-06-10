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
# table output
library(stargazer)
library(xtable)


# Esimators -----------------------------------------------------------------------------------
source("estimators/utils.R")

source("estimators/ate/ipw.R")
source("estimators/ate/matching.R")
source("estimators/ate/AIPW.R")

source("estimators/att/ipw.R")
source("estimators/att/matching.R")
source("estimators/att/AIPW.R")


# Setups --------------------------------------------------------------------------------------

# imput
input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"

# output
output_table_path <- "results/tables/estimator-comparison-death30d.csv"
output_plot_path <- "results/figures/estimator-comparison-death30d.png"

output_tex_ate_path <- "results/tables/ate-comparison-death30d.tex"
output_tex_att_path <- "results/tables/att-comparison-death30d.tex"

# bootstrapping CI
bootstrap_B <- 200
bootstrap_seed <- 2026
conf_level <- 0.95


# Main ----------------------------------------------------------------------------------------

rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)

outcome_var <- "death30d"
treat_var <- var_sets$treatment
covariates <- var_sets$adjustment_sets$basic


# Estimator helpers ---------------------------------------------------------------------------

estimate_naive <- function(data, outcome_var, treat_var) {
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]
    tau <- mean(Y[A == 1]) - mean(Y[A == 0])
    return(tau)
}

estimate_naive_se <- function(data, outcome_var, treat_var) {
    A <- data[[treat_var]]
    Y <- data[[outcome_var]]
    Y1 <- Y[A == 1]
    Y0 <- Y[A == 0]
    se <- sqrt(var(Y1) / length(Y1) + var(Y0) / length(Y0))
    return(se)
}

estimate_regression_ate <- function(data, outcome_var, treat_var, covariates) {
    fit <- lm(
        make_outcome_formula(outcome_var, treat_var, covariates),
        data = data
    )
    # NOTE: transform will fail since treat_var is character
    data1 <- copy(data)
    data0 <- copy(data)
    data1[, (treat_var) := 1]
    data0[, (treat_var) := 0]

    mu1 <- predict(fit, newdata = data1)
    mu0 <- predict(fit, newdata = data0)
    tau <- mean(mu1 - mu0)
    return(tau)
}

estimate_regression_att <- function(data, outcome_var, treat_var, covariates) {
    fit <- lm(
        make_outcome_formula(outcome_var, treat_var, covariates),
        data = data
    )

    treated_data <- data[data[[treat_var]] == 1, ]
    treated_data1 <- copy(treated_data)
    treated_data0 <- copy(treated_data)

    # NOTE: prolly should use raw Y(1) instead of \hatY(1)
    # treated_data1[, (treat_var) := 1]
    treated_data0[, (treat_var) := 0]

    # mu1 <- predict(fit, newdata = treated_data1)
    mu1 <- treated_data1[[outcome_var]]
    mu0 <- predict(fit, newdata = treated_data0)

    tau <- mean(mu1 - mu0)
    return(tau)
}

extract_matching_se <- function(matching_fit) {
    #NOTE: handle the error
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
    }

    res <- data.table(
        Estimand = estimand,
        Method = "PS matching",
        Estimate = as.numeric(matching_result$estimate),
        SE = extract_matching_se(matching_result$fit)
    )
    return(res)
}

estimate_all_methods <- function(
    data,
    outcome_var,
    treat_var,
    covariates
) {
    comparison <- data.table(
        Estimand = c("Naive", rep("ATE", 3), rep("ATT", 3)),
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

    matching_rows <- rbindlist(list(
        estimate_matching_row(data, outcome_var, treat_var, covariates, "ATE"),
        estimate_matching_row(data, outcome_var, treat_var, covariates, "ATT")
    ))

    res <- rbindlist(list(comparison, matching_rows), fill = TRUE)

    return(res)
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

# point estimation
bootstrap_methods <- estimate_bootstrap_methods(
    rhc,
    outcome_var = outcome_var,
    treat_var = treat_var,
    covariates = covariates
)[, .(Estimand, Method)]
# se
bootstrap_methods[, SE := apply(bootstrap_estimates, 1, sd)]

comparison <- merge(
    comparison,
    bootstrap_methods,
    by = c("Estimand", "Method"),
    all.x = TRUE
)
# combine non-bs and bs results
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

comparison_table <- comparison[, .(
    Estimand,
    Method,
    Estimate,
    SE,
    `95%CI`
)]

fwrite(comparison_table, output_table_path)


# Post-processing -----------------------------------------------------------------------------
#' For tex table results

# 1. ATE results table

ATE_RES_LABEL <- "tab:ate_estimates"

comparison_df <- fread(output_table_path)

ate_df <- comparison_df[Estimand != "ATT"]
ate_df <- ate_df[order(factor(
    Method,
    levels = c(
        "Naive comparison",
        "Regression adjustment",
        "PS matching",
        "Hajek IPW",
        "AIPW"
    )
))]


x_tab <- xtable(
    ate_df,
    caption = "ATE by different causal estimators",
    label = ATE_RES_LABEL,
    align = "rllccc",
    digits = 4
)
# add notes
note_text <- "\\hline \n \\multicolumn{5}{l}{\\small Note: The SE for IPW, AIPW were obtained with bootstrapping with 200 times.} \n"
# output
print(
    x_tab,
    include.rownames = FALSE,
    booktabs = TRUE,
    caption.placement = "top",
    add.to.row = list(pos = list(nrow(x_tab)), command = note_text), # <- 核心在這
    file = output_tex_ate_path
)


# 2. ATT results table

ATT_RES_LABEL <- "tab:att_estimates"

att_df <- comparison_df[Estimand == "ATT"]
att_df <- att_df[order(factor(
    Method,
    levels = c(
        "Regression adjustment",
        "PS matching",
        "Hajek IPW",
        "AIPW"
    )
))]


x_tab <- xtable(
    att_df,
    caption = "ATT by different causal estimators",
    label = ATT_RES_LABEL,
    align = "rllccc",
    digits = 4
)
# add notes
note_text <- "\\hline \n \\multicolumn{5}{l}{\\small Note: The SE for IPW, AIPW were obtained with bootstrapping with 200 times.} \n"
# output
print(
    x_tab,
    include.rownames = FALSE,
    booktabs = TRUE,
    caption.placement = "top",
    add.to.row = list(pos = list(nrow(x_tab)), command = note_text),
    file = output_tex_att_path
)

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
            "IPW, AIPW use %s bootstrap resamples",
            bootstrap_B
        )
    ) +
    theme_minimal()
comparison_plot

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
