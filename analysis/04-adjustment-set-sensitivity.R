#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-05-28
#'
#' Desc:
#' -------
#' This file checks whether ATE and ATT estimates are sensitive to moving from
#' the minimal to basic to detailed adjustment set. The main comparison focuses
#' on Hajek IPW and AIPW, with bootstrap CIs.
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
source("estimators/ate/AIPW.R")

source("estimators/att/ipw.R")
source("estimators/att/AIPW.R")


# Setups --------------------------------------------------------------------------------------

input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"

output_table_path <- "results/tables/adjustment-set-sensitivity-death30d.csv"
output_plot_path <- "results/figures/adjustment-set-sensitivity-death30d.png"

bootstrap_B <- 50
bootstrap_seed <- 2026
conf_level <- 0.95


# Main ----------------------------------------------------------------------------------------

rhc <- fread(input_data_path)
var_sets <- read_rds(var_set_path)

outcome_var <- "death30d"
treat_var <- var_sets$treatment
adjustment_set_order <- c("minimal", "basic", "detailed")
adjustment_sets <- var_sets$adjustment_sets[adjustment_set_order]


# Helpers -------------------------------------------------------------------------------------

estimate_adjustment_methods <- function(data, covariates) {
    data.table(
        estimand = c("ATE", "ATE", "ATT", "ATT"),
        estimator = rep(c("Hajek IPW", "AIPW"), 2),
        estimate = c(
            estimate_ate_hajek_ipw(data, outcome_var, treat_var, covariates),
            estimate_ate_aipw(data, outcome_var, treat_var, covariates),
            estimate_att_hajek_ipw(data, outcome_var, treat_var, covariates),
            estimate_att_aipw(data, outcome_var, treat_var, covariates)
        )
    )
}


# Estimation ----------------------------------------------------------------------------------

# point estimate
results <- rbindlist(
    lapply(names(adjustment_sets), function(set_name) {
        covariates <- adjustment_sets[[set_name]]
        estimates <- estimate_adjustment_methods(rhc, covariates)
        estimates[,
            `:=`(
                adjustment_set = set_name,
                n_covariates = length(covariates)
            )
        ]
        return(estimates)
    })
)

# se calc.
s <- Sys.time()
set.seed(bootstrap_seed)
bs_results <- lapply(
    1:bootstrap_B,
    function(b) {
        idx <- sample(nrow(rhc), replace = TRUE)
        bs_data <- rhc[idx, ]

        res <- lapply(names(adjustment_sets), function(set_name) {
            covariates <- adjustment_sets[[set_name]]
            estimates <- estimate_adjustment_methods(bs_data, covariates)
            estimates[, `:=`(bootstrap_id = b, adjustment_set = set_name)]
            return(estimates)
        }) |>
            rbindlist()
        return(res)
    }
) |>
    rbindlist()
e <- Sys.time()


bs_se_summary <- bs_results[,
    .(
        SE = sd(estimate, na.rm = TRUE)
    ),
    by = .(adjustment_set, estimand, estimator)
]

results <- merge(
    results,
    bs_se_summary,
    by = c("adjustment_set", "estimand", "estimator"),
    all.x = TRUE
)

z_value <- qnorm(1 - (1 - conf_level) / 2)
results[, ci_low := estimate - z_value * SE]
results[, ci_high := estimate + z_value * SE]
results[, `95%CI` := sprintf("[%.4f, %.4f]", ci_low, ci_high)]

basic_reference <- results[
    adjustment_set == "basic",
    .(estimand, estimator, basic_estimate = estimate)
]
results <- merge(
    results,
    basic_reference,
    by = c("estimand", "estimator"),
    all.x = TRUE
)
results[, delta_vs_basic := estimate - basic_estimate]

results[,
    adjustment_set := factor(
        adjustment_set,
        levels = adjustment_set_order
    )
]
setorder(results, estimand, estimator, adjustment_set)

fwrite(
    results[,
        .(
            adjustment_set,
            estimand,
            estimator,
            estimate,
            SE,
            `95%CI`,
            delta_vs_basic
        )
    ],
    output_table_path
)


# Plot ----------------------------------------------------------------------------------------

plot <- ggplot(
    results,
    aes(x = adjustment_set, y = estimate, color = estimator, group = estimator)
) +
    geom_errorbar(
        aes(ymin = ci_low, ymax = ci_high),
        width = 0.16,
        linewidth = 0.65,
        position = position_dodge(width = 0.24)
    ) +
    geom_line(position = position_dodge(width = 0.24), linewidth = 0.5) +
    geom_point(position = position_dodge(width = 0.24), size = 2.0) +
    geom_hline(yintercept = 0, color = "grey45") +
    facet_wrap(~estimand) +
    labs(
        x = "Adjustment set",
        y = "Estimated risk difference",
        color = "Estimator",
        title = "Sensitivity to adjustment-set choice",
        subtitle = sprintf(
            "Minimal to basic to detailed; 95%% CIs from %s bootstrap resamples",
            bootstrap_B
        )
    ) +
    theme_minimal()

ggsave(
    filename = output_plot_path,
    plot = plot,
    width = 8,
    height = 5,
    dpi = 300
)

rm(list = ls())
invisible(gc())
cat("...done\n")
