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
    # RHC data
    "Hmisc",
    "ATbounds",
    
    # data wrangling
    "data.table",
    "dplyr",
    "tidyverse",
    
    # estimation, analysis
    "fixest",
    "sensemakr",
    
    # causal
    "dagitty",
    "MatchIt",
    "Matching",
    
    # DAG
    "dagitty",
    "ggdag"
)


for (.pkg in .required_pkgs) {
    if (.pkg %in% installed.packages()) {
        next
    }
    tryCatch({
        install.packages(.pkg)  
    }, error = function(err) {
        cat(sprintf("got error when installing pkg=%s", .pkg))
    }
    )
}

# invisible(lapply(required_packages, library, character.only = TRUE))


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


cat("...done")
