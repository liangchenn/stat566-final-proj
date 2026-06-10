#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-05-27
#'
#' Desc:
#' -------
#' This file creates a descriptive statistics table for the constructed RHC
#' dataset, comparing No RHC vs RHC groups using the basic adjustment variables.
#'
#' Input
#' -----
#' - data/processed/rhc-processed-data.csv
#' - data/processed/rhc-var-sets.rds
#'
#' Output
#' ------
#' - results/tables/descriptive-stat-basic-by-rhc.tex
#' - results/tables/descriptive-stat-basic-by-rhc.csv

.FILE_NAME <- "04-descriptive-stat.R"
cat(sprintf("...Running %s ...", .FILE_NAME))


# Packages ------------------------------------------------------------------------------------
library(data.table)
library(readr)
library(stargazer)


# Main ----------------------------------------------------------------------------------------

output_tex_path <- "results/tables/descriptive-stat-basic-by-rhc.tex"
output_csv_path <- "results/tables/descriptive-stat-basic-by-rhc.csv"

input_data_path <- "data/processed/rhc-processed-data.csv"
var_set_path <- "data/processed/rhc-var-sets.rds"
var_desc_path <- "data/RHC-variable-desc.csv"


rhc <- fread(input_data_path)
var_sets <- readr::read_rds(var_set_path)
basic_vars <- c(
    var_sets$adjustment_sets$basic,
    var_sets$outcomes[1:2]
)



# Variable labels -----------------------------------------------------------------------------

variable_desc <- if (file.exists(var_desc_path)) {
    fread(var_desc_path)
} else {
    data.table(variable_name = character(), desc = character())
}

label_lookup <- setNames(variable_desc$desc, variable_desc$variable_name)

custom_labels <- c(
    age = "Age",
    edu = "Years of education",
    sex_Female = "Female",
    race_black = "Black race",
    race_other = "Other race",
    income1 = "Income $11k-$25k",
    income2 = "Income $25k-$50k",
    income3 = "Income >$50k",
    aps1 = "APACHE score",
    pafi1 = "PaO2/FIO2 ratio",
    meanbp1 = "Mean blood pressure",
    scoma1 = "Glasgow coma score",
    dnr1_Yes = "DNR on day 1",
    cardiohx = "Cardiovascular symptoms",
    chfhx = "Congestive heart failure",
    amihx = "Myocardial infarction",
    death30d = "Death (30D)",
    death180d = "Death (180D)"
)

make_label <- function(var) {
    if (var %in% names(custom_labels)) {
        return(unname(custom_labels[var]))
    }
    if (var %in% names(label_lookup) && !is.na(label_lookup[[var]]) && nzchar(label_lookup[[var]])) {
        label <- label_lookup[[var]]
        label <- sub("^1 if ", "", label)
        label <- sub(", and 0 otherwise.*$", "", label)
        label <- sub(" \\(Omitted category.*$", "", label)
        return(label)
    }
    var
}


# Summary helpers -----------------------------------------------------------------------------

is_binary <- function(x) {
    vals <- na.omit(unique(x))
    length(vals) <= 2 && all(vals %in% c(0, 1))
}

fmt <- function(x, digits = 3) {
    formatC(x, format = "f", digits = digits)
}

smd_continuous <- function(x, a) {
    x1 <- x[a == 1]
    x0 <- x[a == 0]
    pooled_sd <- sqrt((stats::var(x1, na.rm = TRUE) + stats::var(x0, na.rm = TRUE)) / 2)
    if (is.na(pooled_sd) || pooled_sd == 0) return(0)
    (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / pooled_sd
}

smd_binary <- function(x, a) {
    p1 <- mean(x[a == 1] == 1, na.rm = TRUE)
    p0 <- mean(x[a == 0] == 1, na.rm = TRUE)
    pooled_sd <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
    if (is.na(pooled_sd) || pooled_sd == 0) return(0)
    (p1 - p0) / pooled_sd
}

summarize_one <- function(var) {
    x <- rhc[[var]]
    a <- rhc$RHC
    binary <- is_binary(x)

    if (binary) {
        no_rhc_mean <- mean(x[a == 0] == 1, na.rm = TRUE)
        rhc_mean <- mean(x[a == 1] == 1, na.rm = TRUE)
        no_rhc_se <- sqrt(no_rhc_mean * (1 - no_rhc_mean) / sum(!is.na(x[a == 0])))
        rhc_se <- sqrt(rhc_mean * (1 - rhc_mean) / sum(!is.na(x[a == 1])))
        smd <- smd_binary(x, a)
    } else {
        no_rhc_mean <- mean(x[a == 0], na.rm = TRUE)
        rhc_mean <- mean(x[a == 1], na.rm = TRUE)
        no_rhc_se <- stats::sd(x[a == 0], na.rm = TRUE) / sqrt(sum(!is.na(x[a == 0])))
        rhc_se <- stats::sd(x[a == 1], na.rm = TRUE) / sqrt(sum(!is.na(x[a == 1])))
        smd <- smd_continuous(x, a)
    }

    data.table(
        Variable = make_label(var),
        `No RHC` = sprintf("%s (%s)", fmt(no_rhc_mean), fmt(no_rhc_se)),
        RHC = sprintf("%s (%s)", fmt(rhc_mean), fmt(rhc_se))
    )
}

n_control <- sum(rhc$RHC == 0)
n_treated <- sum(rhc$RHC == 1)

desc_table <- rbind(
    data.table(
        Variable = "Observations",
        `No RHC` = as.character(n_control),
        RHC = as.character(n_treated)
    ),
    rbindlist(lapply(basic_vars, summarize_one))
)

# save table
fwrite(desc_table, output_csv_path)


# LaTeX output ----------------------------------------------------------------------

# output tex table
stargazer(
    as.data.frame(desc_table),
    type = "latex",
    summary = FALSE,
    rownames = FALSE,
    title = paste0(
        "Descriptive statistics \\ Basic adjustment variables by RHC treatment group. "
    ),
    label = "tab:rhc-desc-stat",
    notes = paste0(
        "NOTE: Cells report mean or proportion with standard error in parentheses."
    ),
    out = output_tex_path
) |> capture.output()


rm(list = ls()); invisible(gc())
# cat(sprintf("\nSaved descriptive CSV to %s\n", output_csv_path))
# cat(sprintf("Saved stargazer LaTeX table to %s\n", output_tex_path))
cat("...done\n")
