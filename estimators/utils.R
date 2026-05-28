#' README:
#' -------
#' - author: Liang-Cheng Chen
#' - date: 2026-05-27
#'
#' Desc:
#' -------
#' utility functions for causal estimators.
#' This file contains nuiance helpers for 
#' outcome regression, IPW, AIPW, and matching estimators.


# Packages ------------------------------------------------------------------------------------
library(data.table)


# Formula helpers -----------------------------------------------------------------------------

make_rhs <- function(covariates) {
    if (length(covariates) == 0) return("1")
    return(paste(covariates, collapse = " + "))
}

make_formula <- function(lhs, covariates) {
    fmla <- as.formula(paste(lhs, "~", make_rhs(covariates)))
    return(fmla)
}

make_outcome_formula <- function(outcome_var, treat_var, covariates) {
    if (length(covariates) == 0) {
        return(as.formula(paste(outcome_var, "~", treat_var)))
    }
    
    fmla <- as.formula(
        paste0(outcome_var, " ~ ", treat_var, " + ", make_rhs(covariates))
    )
    return(fmla)
}

make_interacted_outcome_formula <- function(outcome_var, treat_var, covariates) {
    if (length(covariates) == 0) {
        return(as.formula(paste(outcome_var, "~", treat_var)))
    }

    fmla <- as.formula(
        paste0(outcome_var, " ~ ", treat_var, " * (", make_rhs(covariates), ")")
    )
    return(fmla)
}


# Shared nuisance models ----------------------------------------------------------------------

fit_propensity_score <- function(data,
                                 treat_var,
                                 covariates,
                                 trim = c(0.02, 0.98),
                                 return_fit = FALSE) {
    ps_formula <- make_formula(treat_var, covariates)
    ps_fit <- glm(ps_formula, data = data, family = binomial("logit"))
    ps <- predict(ps_fit, type = "response")
    ps <- pmin(pmax(ps, trim[1]), trim[2])

    if (return_fit) {
        return(list(fit = ps_fit, ps = ps))
    }

    return(ps)
}

fit_outcome_model <- function(data,
                              outcome_var,
                              treat_var,
                              covariates,
                              formula_func=make_outcome_formula
) {
    outcome_formula <- formula_func(outcome_var, treat_var, covariates)
    fit <- lm(outcome_formula, data = data)
    return(fit)
}
