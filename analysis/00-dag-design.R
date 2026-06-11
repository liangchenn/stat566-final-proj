#' README:
#' -------
#' - author: Andy Shin, Jian Kang
#' - date: 2026-04-27
#'
#' Desc:
#' -------
#' This file creates the DAG, representing the causal assumptions for RHC dataset.
#'
#' Output:
#' -------
#' - results/figures/rhc-dag-ggdag.png
#' - results/figures/rhc-dag-ggdag.pdf

# Packages ------------------------------------------------------------------------------------
library(data.table)
library(dagitty)
library(ggdag)
library(ggplot2)


# Main ----------------------------------------------------------------------------------------

output_png_path <- "results/figures/rhc-dag-ggdag.png"
output_pdf_path <- "results/figures/rhc-dag-ggdag.pdf"

# load data
rhc <- fread("data/processed/rhc-processed-data.csv")

# dag construction
rhc_dag <- dagitty(
    "dag {
  cat1 -> pafi1;
  cat1 -> scoma1;
  cat1 -> RHC;
  age -> chfhx;
  age -> dnr1_Yes;
  age -> meanbp1;
  chfhx -> survival;
  chfhx -> meanbp1;
  chfhx -> pafi1;
  meanbp1 -> aps1;
  meanbp1 -> survival;
  meanbp1 -> pafi1;
  meanbp1 -> scoma1;
  meanbp1 -> RHC;
  pafi1 -> aps1;
  pafi1 -> dnr1_Yes;
  pafi1 -> RHC;
  scoma1 -> aps1;
  scoma1 -> dnr1_Yes;
  aps1 -> survival;
  aps1 -> RHC;
  dnr1_Yes -> survival;
  RHC -> survival
}"
)

# plotting
coordinates(rhc_dag) <- list(
    x = c(
        cat1 = 0,
        age = 0,
        chfhx = 0,
        pafi1 = 1,
        meanbp1 = 1,
        scoma1 = 1,
        aps1 = 2,
        dnr1_Yes = 2,
        RHC = 3,
        survival = 3
    ),

    # Staggered Y-coordinates to dodge straight lines
    y = c(
        cat1 = 0,
        age = 2,
        chfhx = 4,
        pafi1 = -1,
        meanbp1 = 2,
        scoma1 = 5,
        aps1 = 0.5,
        dnr1_Yes = 3.5,
        RHC = 0,
        survival = 2.5
    )
)

# 1. Base dagitty plot
# plot(rhc_dag)

# 2. Presentation-ready plot using ggdag
dag_plot <- ggdag(rhc_dag, text = FALSE, use_labels = "name") +
    theme_dag() +
    ggtitle("DAG model of RHC study") +
    theme(plot.title = element_text(hjust = 0.5))


ggsave(
    filename = output_png_path,
    plot = dag_plot,
    width = 8,
    height = 5.5,
    dpi = 300
)

ggsave(
    filename = output_pdf_path,
    plot = dag_plot,
    width = 8,
    height = 5.5
)

# print(dag_plot)
