#' README:
#' -------
#' - author: ...
#' - date: 2026-05-27
#' 
#' Desc:
#' -------
#' This file is the main entry point of the project.
#' Run this file and the codes will generate all results we include in our write-up.
#' 
#' Output
#' --------
#' -

setwd(".")

# Project Preparations ------------------------------------------------------------------------

cat("Preparing project data: \n")
source("scripts/00-setups.R")
source("scripts/01-load-raw-data.R")
source("scripts/02-define-variables.R")
source("scripts/03-construct-data.R")
source("scripts/04-descriptive-stat.R")
source("scripts/05-smd-love-plot.R")
cat("finished \n\n")


# Project Analysis ----------------------------------------------------------------------------

cat("Preparing project analysis: \n")
source("analysis/00-dag-design.R")
source("analysis/00-estimator-comparison.R")
source("analysis/01-balance-diagnostics.R")
source("analysis/02-odds-ratio-estimator-comparison.R")
source("analysis/03-trimming-sensitivity.R")
source("analysis/04-adjustment-set-sensitivity.R")
source("analysis/05-cate-subgroup-analysis.R")
cat("finished \n\n")
