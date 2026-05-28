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
source("scripts/01-load-raw-data.R")
source("scripts/02-define-variables.R")
source("scripts/03-construct-data.R")
source("scripts/04-descriptive-stat.R")
source("scripts/05-smd-love-plot.R")
cat("finished \n\n")


# Project Analysis ----------------------------------------------------------------------------

cat("Preparing project analysis: \n\n")

