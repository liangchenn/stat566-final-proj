#' README:
#' -------
#' - author: ...
#' - date: 2026-05-27
#'
#' Desc:
#' -------
#' This file creates raw data absolute SMD love plot for basic adjustment variables.
#'
#' Input
#' -----
#' - data/processed/rhc-processed-data.csv
#' - data/processed/rhc-var-sets.rds
#'
#' Output
#' ------
#' - results/figures/raw-smd-love-plot.png
#' - results/tables/raw-smd-basic.csv

.FILE_NAME <- "05-smd-love-plot.R"
cat(sprintf("...Running %s ...", .FILE_NAME))

# Packages ------------------------------------------------------------------------------------
library(data.table)
library(readr)
library(ggplot2)

# Main ----------------------------------------------------------------------------------------

# setups
input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"
output_plot_path <- "results/figures/raw-smd-love-plot.png"
output_table_path <- "results/tables/raw-smd-basic.csv"

# load data
rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)
basic_vars <- var_sets$adjustment_sets$basic


# SMD plotting
smd <- function(data, covs, treat_var) {
    sapply(covs, function(v) {
        x1 <- data[get(treat_var) == 1, get(v)]
        x0 <- data[get(treat_var) == 0, get(v)]
        s1 <- var(x1)
        s0 <- var(x0)

        abs(mean(x1) - mean(x0)) / sqrt((s1 + s0) / 2)
    })
}

smd_pre <- smd(
    data = rhc,
    covs = basic_vars,
    treat_var = "RHC"
)

smd_data <- data.table(
    covariate = names(smd_pre),
    smd_pre = as.numeric(smd_pre)
)

smd_data[, covariate := factor(covariate, levels = covariate[order(smd_pre)])]

# output data
fwrite(smd_data[order(-smd_pre)], output_table_path)

# plotting
love_plot <- ggplot(smd_data) +
    geom_point(aes(x = covariate, y = smd_pre, color = "pre"), size = 2.4) +
    geom_hline(yintercept = 0.1, color = "grey", linetype = "dashed") +
    scale_color_manual(values = c("pre" = "tomato")) +
    labs(color = "raw data") +
    ylab("Absolute SMD") +
    xlab(NULL) +
    coord_flip() +
    theme_minimal()+
    ggtitle("Absolute SMD Plot for Basic Adjustments")

ggsave(
    filename = output_plot_path,
    plot = love_plot,
    width = 7,
    height = 5,
    dpi = 300
)

rm(list = ls()); invisible(gc())

cat("...done\n")
