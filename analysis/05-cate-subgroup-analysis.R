#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-05-28
#'
#' Desc:
#' -------
#' This file performs a simple exploratory CATE analysis by estimating treatment
#' effects within clinically interpretable subgroups. Estimates use the basic
#' adjustment set, with the subgroup-defining covariate removed from the nuisance
#' models. The main comparison focuses on Hajek IPW and AIPW with bootstrap CIs.
#'
#' Input
#' -----
#' - data/processed/rhc-processed-data.csv
#' - data/processed/rhc-var-sets.rds
#'
#' Output
#' ------
#' - results/tables/cate-subgroup-death30d.csv
#' - results/figures/cate-subgroup-death30d.png

.FILE_NAME <- "05-cate-subgroup-analysis.R"
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

output_table_path <- "results/tables/cate-subgroup-death30d.csv"
output_plot_path <- "results/figures/cate-subgroup-death30d.png"

bootstrap_B <- 100
bootstrap_seed <- 2026
conf_level <- 0.95


# Main ----------------------------------------------------------------------------------------


rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)

outcome_var <- "death30d"
treat_var <- var_sets$treatment
covariates <- var_sets$adjustment_sets$basic

# create sub-groups

# 1. age groups
rhc[, age_group := fifelse(age >= 65, "65+", fifelse(age < 50, "<50", "50-65"))]

# 2. day1 vital severity groups
rhc[, aps1_group := fifelse(
        aps1 >= median(aps1, na.rm = TRUE),
        "high APS",
        "low APS")
]
rhc[
    ,
    pafi1_group := fifelse(
        pafi1 >= median(pafi1, na.rm = TRUE),
        "High PaO2/FIO2",
        "Low PaO2/FIO2"
    )
]
rhc[, dnr_group := fifelse(dnr1_Yes == 1, "DNR", "No DNR")]
rhc[, chfhx_group := fifelse(chfhx == 1, "CHF history", "No CHF history")]

subgroup_specs <- list(
    age = list(group_var = "age_group", remove_covariate = "age"),
    aps1 = list(group_var = "aps1_group", remove_covariate = "aps1"),
    pafi1 = list(group_var = "pafi1_group", remove_covariate = "pafi1"),
    dnr = list(group_var = "dnr_group", remove_covariate = "dnr1_Yes"),
    chfhx = list(group_var = "chfhx_group", remove_covariate = "chfhx")
)


# Helpers -------------------------------------------------------------------------------------

estimate_subgroup_methods <- function(
    data,
    subgroup_name,
    subgroup_var,
    remove_covariate,
    level
) {
    subgroup_data <- data[get(subgroup_var) == level]
    subgroup_covariates <- setdiff(covariates, remove_covariate)

    n_treated <- sum(subgroup_data[[treat_var]] == 1)
    n_control <- sum(subgroup_data[[treat_var]] == 0)

    if (n_treated == 0 || n_control == 0) {
        return(data.table())
    }

    res <- data.table(
        subgroup = subgroup_name,
        level = level,
        estimand = c("ATE", "ATE", "ATT", "ATT"),
        estimator = rep(c("Hajek IPW", "AIPW"), 2),
        estimate = c(
            estimate_ate_hajek_ipw(
                subgroup_data,
                outcome_var,
                treat_var,
                subgroup_covariates
            ),
            estimate_ate_aipw(
                subgroup_data,
                outcome_var,
                treat_var,
                subgroup_covariates
            ),
            estimate_att_hajek_ipw(
                subgroup_data,
                outcome_var,
                treat_var,
                subgroup_covariates
            ),
            estimate_att_aipw(
                subgroup_data,
                outcome_var,
                treat_var,
                subgroup_covariates
            )
        ),
        n = nrow(subgroup_data),
        n_treated = n_treated,
        n_control = n_control,
        n_covariates = length(subgroup_covariates)
    )
    
    return(res)
}


estimate_all_subgroups <- function(data) {

    rbindlist(lapply(names(subgroup_specs), function(subgroup_name) {
        subgroup_var <- subgroup_specs[[subgroup_name]]$group_var
        remove_covariate <- subgroup_specs[[subgroup_name]]$remove_covariate

        rbindlist(lapply(unique(data[[subgroup_var]]), function(level) {
            estimate_subgroup_methods(
                data,
                subgroup_name,
                subgroup_var,
                remove_covariate,
                level
            )
        }), fill = TRUE)
    }), fill = TRUE)
}


estimate_all_subgroups <- function(df) {
    
    res <- lapply(
        names(subgroup_specs), function(subgroup_name) {
            subgroup_var <- subgroup_specs[[subgroup_name]]$group_var
            remove_covariate <- subgroup_specs[[subgroup_name]]$remove_covariate
            
            subdf <- lapply(
                unique(df[[subgroup_var]]), function(level) {
                    estimate_subgroup_methods(
                        df, subgroup_name, subgroup_var,
                        remove_covariate, level
                    )
                }
            ) |> rbindlist(fill = TRUE)
            
            return(subdf)
        }
    ) |> rbindlist(fill = TRUE)
    
    
}




# Estimation ----------------------------------------------------------------------------------

results <- estimate_all_subgroups(rhc)

set.seed(bootstrap_seed)
bootstrap_results <- lapply(1:bootstrap_B, function(b) {
    idx <- sample(nrow(rhc), replace = TRUE)
    bs_data <- rhc[idx, ]
    
    estimates <- estimate_all_subgroups(bs_data)
    estimates[, bootstrap_id := b]
    
    return(estimates)
}) |> rbindlist(fill = TRUE)

bootstrap_summary <- bootstrap_results[
    , .(SE = sd(estimate, na.rm = TRUE))
    , by = .(subgroup, level, estimand, estimator)
]

results <- merge(
    results,
    bootstrap_summary,
    by = c("subgroup", "level", "estimand", "estimator"),
    all.x = TRUE
)

z_value <- qnorm(1 - (1 - conf_level) / 2)
results[, ci_low := estimate - z_value * SE]
results[, ci_high := estimate + z_value * SE]
results[, `95%CI` := sprintf("[%.4f, %.4f]", ci_low, ci_high)]

setorder(results, subgroup, level, estimand, estimator)

fwrite(
    results[
        ,
        .(
            subgroup,
            level,
            estimand,
            estimator,
            estimate,
            SE,
            `95%CI`,
            n,
            n_treated,
            n_control
        )
    ],
    output_table_path
)


# Plot ----------------------------------------------------------------------------------------

results[, subgroup_label := paste(subgroup, level, sep = ": ")]
results[, subgroup_label := factor(subgroup_label, levels = rev(unique(subgroup_label)))]

cate_plot <- ggplot(results)+
    aes(x = subgroup_label, y = estimate, color = estimator)+
    geom_errorbar(
        aes(ymin = ci_low, ymax = ci_high),
        width = 0.18,
        linewidth = 0.6,
        position = position_dodge(width = 0.45)
    ) +
    geom_point(position = position_dodge(width = 0.45), size = 2) +
    geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
    facet_wrap(~ estimand) +
    coord_flip() +
    labs(
        x = NULL,
        y = "Estimated risk difference",
        color = "Estimator",
        title = "Subgroup Average Treatment Effect Estimates (CATE)",
        subtitle = sprintf(
            "Basic adjustment set; 95%% CIs from %s bootstrap resamples",
            bootstrap_B
        )
    ) +
    theme_minimal()

ggsave(
    filename = output_plot_path,
    plot = cate_plot,
    width = 8,
    height = 6,
    dpi = 300
)

rm(list = ls())
invisible(gc())
cat("...done\n")
