#' README:
#' -------
#' - author: ...
#' - date: 2026-05-26
#' 
#' Desc:
#' -------
#' This file defines important variables and adjustment sets for estimation.
#' 
#' Output
#' --------
#' - data/processed/rhc-adjustment-sets.rds


.FILE_NAME <- "02-define-variables.R"
cat(sprintf("...Running %s ...", .FILE_NAME))


# Main ----------------------------------------------------------------------------------------

# output file
var_set_path <- "./data/processed/rhc-var-sets.rds"

# Variable groups -----------------------------------------------------------------------------

# Outcomes
outcome_vars <- c(
    "death30d", "death180d",
    "survival30d", "survival180d"
)

# treatment
treatment_var <- c("RHC")

# Demographics
basic_demographic_vars <- c(
    "age", "edu", "sex_Female",
    "race_black", "race_other",
    "income1", "income2", "income3" # baseline group: <$11k
)

demographic_vars <- c(
    basic_demographic_vars,
    
    "ninsclas_Medicaid", "ninsclas_Medicare",
    "ninsclas_Medicare_and_Medicaid", "ninsclas_No_insurance",
    "ninsclas_Private_and_Medicare"
) |> unique()


# Clinical history before RHC assignment.
basic_comorbidity_vars <- c(
    "cardiohx", # cardiovascular symptons
    "chfhx",    # congestive heart failure
    "amihx"     # myocardial infarction
)

comorbidity_vars <- c(
    basic_comorbidity_vars,
    "cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx",
    "renalhx", "liverhx", "gibledhx", "malighx", "immunhx",
    "transhx", "amihx"
) |> unique()


# Day 1 physiology diagnosis

basic_physio_vars <- c(
    "aps1",    # APACHE score, higher->worse
    "pafi1",   # PaO2/FIO2 ratio, blood O2 level (?)
    "meanbp1", # mean blood pressure
    "scoma1",  # coma score
    "dnr1_Yes" # Do Not Resuscitate status
)

physio_vars <- c(
    basic_physio_vars,
    "aps1", "scoma1", "meanbp1", "resp1", "hrt1", "pafi1",
    "paco21", "ph1", "wblc1", "hema1", "sod1", "pot1",
    "crea1", "bili1", "alb1", "temp1", "wtkilo1", "wt0"
) |> unique()

# Disease category / diagnosis indicators. 
diagnosis_vars <- c(
    "cat1_CHF", "cat1_Cirrhosis", "cat1_Colon_Cancer", "cat1_Coma",
    "cat1_COPD", "cat1_Lung_Cancer", "cat1_MOSF_Malignancy",
    "cat1_MOSF_Sepsis",
    "cat2_Cirrhosis", "cat2_Colon_Cancer", "cat2_Coma",
    "cat2_Lung_Cancer", "cat2_MOSF_Malignancy", "cat2_MOSF_Sepsis",
    "resp_Yes", "card_Yes", "neuro_Yes", "gastr_Yes", "renal_Yes",
    "meta_Yes", "hema_Yes", "seps_Yes", "trauma_Yes", "ortho_Yes",
    "ca_Metastatic", "ca_Yes"
)

# other variables
other_vars <- c("surv2md1", "das2d3pc")



# Three adjustment sets definitions:
# - minimal: demographics only, mininmal baseline.
# - basic: core backdoor set
# - detailed: detailed backdoor

minimal_adjustment_vars <- basic_demographic_vars

basic_adjustment_vars <- c(
    minimal_adjustment_vars,
    basic_comorbidity_vars,
    basic_physio_vars
) |> unique()


detailed_adjustment_vars <- c(
    basic_adjustment_vars,
    demographic_vars,
    comorbidity_vars,
    physio_vars,
    diagnosis_vars,
    other_vars
) |> unique()

adjustment_sets <- list(
    minimal = minimal_adjustment_vars,
    basic = basic_adjustment_vars,
    detailed = detailed_adjustment_vars
)

all_adjustment_vars <- unique(unlist(adjustment_sets, use.names = FALSE))


rhc_var_sets <- list(
    # outcomes
    outcomes=outcome_vars,
    # treatment
    treatment=treatment_var,
    # adjustment variable sets
    adjustment_sets=list(
        minimal=minimal_adjustment_vars,
        basic=basic_adjustment_vars,
        detailed=detailed_adjustment_vars,
        all=all_adjustment_vars
    ),
    # demographics
    demographics=list(
        basic=basic_demographic_vars,
        all=demographic_vars
    ),
    # comorbidity
    comorbidity=list(
        basic=basic_comorbidity_vars,
        all=comorbidity_vars
    ),
    # day1 physio
    physio=list(
        basic=basic_physio_vars,
        all=physio_vars
    ),
    # diagnosis
    diagnosis=diagnosis_vars,
    # others
    others=other_vars
)

readr::write_rds(rhc_var_sets, var_set_path)

rm(list = ls()); gc() |> invisible()
cat("...done\n")
