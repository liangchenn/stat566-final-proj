#' README:
#' -------
#' - author: ...
#' - date: 2026-05-26
#' 
#' Desc:
#' -------
#' This file does the following things:
#' 1. load raw data with corrected outcome
#' 2. keep important variables (adjustment sets)
#' 3. save analysis-ready data for estimation
#' 
#' Output
#' --------
#' - data/processed/rhc-data.csv
#' - data/processed/rhc-adjustment-sets.rds

.FILE_NAME <- "03-construct-data.R"
cat(sprintf("...Running %s ...", .FILE_NAME))


# Packages ------------------------------------------------------------------------------------
library(data.table)
library(readr)


# Main ----------------------------------------------------------------------------------------

# output file
output_data_path <- "data/processed/rhc-processed-data.csv"

raw_data_path <- "data/raw/rhc-raw-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"

rhc <- fread(raw_data_path)
var_sets <- readr::read_rds(var_set_path)


ordered_vars <- unique(
    c(
        var_sets$outcomes, 
        var_sets$treatment, 
        var_sets$adjustment_sets$all
    )
)

rhc_processed <- rhc[, ..ordered_vars]

fwrite(rhc_processed, output_data_path)


rm(list = ls()); invisible(gc())
cat("...done\n")
