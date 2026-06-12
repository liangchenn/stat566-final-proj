#' README:
#' -------
#' - author: ...
#' - date: 2026-05-25
#'
#' - Desc:
#'  - This file prepares relevant setups for the project.

.FILE_NAME <- "00-setups.R"
cat(sprintf("...Running %s ...", .FILE_NAME))

# Packages ------------------------------------------------------------------------------------

.required_pkgs <- c(
    # pkg dependencies
    "rlang",
    "devtools"

    # RHC data
    "Hmisc",
    "ATbounds",

    # data wrangling
    "data.table",
    "dplyr",
    "readr",
    "tidyverse",

    # estimation, analysis
    "fixest",
    "sensemakr",
    "stargazer",
    "xtable",
    "EValue",

    # causal
    "dagitty",
    "MatchIt",
    "Matching",

    # DAG
    "dagitty",
    "ggdag",
    "ggplot2"
)


for (.pkg in .required_pkgs) {
    tryCatch(
        {
            install.packages(.pkg)
        },
        error = function(err) {
            cat(sprintf("got error when installing pkg=%s", .pkg))
        }
    )
}

invisible(lapply(.required_pkgs, library, character.only = TRUE))

# for dml.sensemakr
if (!"dml.sensemakr" %in% installed.packages()) {
    devtools::install_github("carloscinelli/dml.sensemakr")
}

# Folder Structure ----------------------------------------------------------------------------

.paths <- c(
    "data/raw",
    "data/processed",
    "results",
    "results/tables",
    "results/figures",
    "estimators"
)

for (.path in .paths) {
    if (!dir.exists(.path)) {
        dir.create(.path, recursive = TRUE)
    }
}


cat("...done\n")
