#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-05-28
#'
#' Desc:
#' -------
#' This file performs a simple exploratory CATE analysis by estimating treatment
#' effects within clinically interpretable subgroups.
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
source("estimators/att/ipw.R")


# Main ----------------------------------------------------------------------------------------



input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"
output_table_path <- "results/tables/cate-subgroup-death30d.csv"
output_plot_path <- "results/figures/cate-subgroup-death30d.png"

rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)

outcome_var <- "death30d"
treat_var <- var_sets$treatment
covariates <- var_sets$adjustment_sets$basic

rhc[, age_group := fifelse(age >= median(age, na.rm = TRUE), "Older", "Younger")]
rhc[, aps1_group := fifelse(aps1 >= median(aps1, na.rm = TRUE), "High APACHE", "Low APACHE")]
rhc[, pafi1_group := fifelse(pafi1 >= median(pafi1, na.rm = TRUE), "High PaO2/FIO2", "Low PaO2/FIO2")]
rhc[, dnr_group := fifelse(dnr1_Yes == 1, "DNR", "No DNR")]
rhc[, chfhx_group := fifelse(chfhx == 1, "CHF history", "No CHF history")]

subgroup_specs <- list(
    age = list(group_var = "age_group", remove_covariate = "age"),
    aps1 = list(group_var = "aps1_group", remove_covariate = "aps1"),
    pafi1 = list(group_var = "pafi1_group", remove_covariate = "pafi1"),
    dnr = list(group_var = "dnr_group", remove_covariate = "dnr1_Yes"),
    chfhx = list(group_var = "chfhx_group", remove_covariate = "chfhx")
)

estimate_subgroup <- function(data, subgroup_name, subgroup_var, remove_covariate, level) {
    subgroup_data <- data[get(subgroup_var) == level]
    subgroup_covariates <- setdiff(covariates, remove_covariate)

    if (length(unique(subgroup_data[[treat_var]])) < 2) {
        return(data.table())
    }

    data.table(
        subgroup = subgroup_name,
        level = level,
        estimand = c("ATE", "ATT"),
        estimator = "Hajek IPW",
        estimate = c(
            estimate_ate_hajek_ipw(subgroup_data, outcome_var, treat_var, subgroup_covariates),
            estimate_att_hajek_ipw(subgroup_data, outcome_var, treat_var, subgroup_covariates)
        ),
        n = nrow(subgroup_data),
        n_treated = sum(subgroup_data[[treat_var]] == 1),
        n_control = sum(subgroup_data[[treat_var]] == 0)
    )
}

results <- rbindlist(lapply(names(subgroup_specs), function(subgroup_name) {
    subgroup_var <- subgroup_specs[[subgroup_name]]$group_var
    remove_covariate <- subgroup_specs[[subgroup_name]]$remove_covariate
    rbindlist(lapply(unique(rhc[[subgroup_var]]), function(level) {
        estimate_subgroup(rhc, subgroup_name, subgroup_var, remove_covariate, level)
    }))
}), fill = TRUE)

fwrite(results, output_table_path)

results[, subgroup_label := paste(subgroup, level, sep = ": ")]
results[, subgroup_label := factor(subgroup_label, levels = rev(unique(subgroup_label)))]

cate_plot <- ggplot(results, aes(x = subgroup_label, y = estimate, color = estimand)) +
    geom_point(size = 2.4) +
    geom_hline(yintercept = 0, color = "grey45", linetype = "dashed") +
    coord_flip() +
    labs(
        x = NULL,
        y = "Estimated risk difference",
        color = "Estimand",
        title = "Exploratory subgroup treatment-effect estimates"
    ) +
    theme_minimal()

ggsave(
    filename = output_plot_path,
    plot = cate_plot,
    width = 7,
    height = 5,
    dpi = 300
)

rm(list = ls()); invisible(gc())
cat("...done\n")
