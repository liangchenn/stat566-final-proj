# Analysis TODO and Report Completion Plan

Updated: 2026-06-03

## Purpose

This document defines the remaining work needed to turn the current RHC analysis into a complete, reproducible causal-inference report. The priorities are:

1. Establish one clearly defined primary analysis with uncertainty estimates.
2. Diagnose the main identifying assumptions.
3. Evaluate robustness to modeling and design choices.
4. Separate confirmatory results from exploratory analyses.

## File Organization Rules

- Put reusable estimation methods in `estimators/`.
- Put data construction, shared configuration, and data-quality checks in `scripts/`.
- Put each report-producing analysis in `analysis/`.
- Save numerical outputs in `results/tables/`.
- Save report-ready plots in `results/figures/`.
- Add every final analysis script to `main.R`.
- Each analysis table should include the outcome, estimand, estimator, adjustment set, point estimate, uncertainty method, standard error when available, confidence interval, and sample size.

## Priority 0: Required for the Final Report

### 1. Define a Reproducible Analysis Specification

**Proposed file:** `scripts/06-define-analysis-specification.R`

**Goal**

Create a single configuration object defining the primary outcome, secondary outcome, primary estimand, adjustment set, propensity-score trimming rule, confidence level, bootstrap repetitions, cross-fitting folds, and random seeds.

**Why this is needed**

These choices are currently distributed across multiple scripts and notebooks. A centralized specification prevents accidental differences between analyses and makes the final report reproducible.

**Expected outputs**

- `data/processed/rhc-analysis-spec.rds`
- A short console summary of the primary and secondary analysis definitions.

**Recommended initial specification**

- Primary outcome: `death30d`
- Secondary outcome: `death180d`
- Primary estimand: ATE
- Secondary estimand: ATT
- Primary adjustment set: `basic`
- Primary propensity-score trim: `[0.02, 0.98]`
- Confidence level: 95%
- Bootstrap repetitions: at least 500 for final tables
- Cross-fitting folds: 5
- Random seed: 2026

### 2. Produce the Primary Effect-Estimate Table with Confidence Intervals

**Proposed file:** `analysis/06-primary-effect-estimates.R`

**Supporting estimator files**

- `estimators/ate/outcome-regression.R`
- `estimators/att/outcome-regression.R`
- `estimators/att/AIPW.R`
- Existing `estimators/ate/AIPW.R`
- Existing IPW and matching estimator files

**Goal**

Estimate the effect of RHC on 30-day mortality using regression adjustment, Hajek IPW, propensity-score matching, and AIPW for the supported estimands. Report point estimates and confidence intervals in one consistent table.

**Why this is needed**

The current estimator-comparison table only reports point estimates. A final report needs uncertainty estimates to distinguish stable findings from sampling variation. AIPW also provides a doubly robust comparison with the existing IPW result.

**Expected outputs**

- `results/tables/primary-effect-estimates-death30d.csv`
- `results/figures/primary-effect-estimates-death30d.png`

**Required table fields**

- Outcome and estimand
- Estimator
- Risk difference
- Standard error
- 95% confidence interval
- Uncertainty method
- Effective sample size when relevant
- Number of observations

**Uncertainty plan**

- AIPW: influence-function asymptotic confidence interval
- Regression adjustment, IPW, and matching: bootstrap confidence intervals
- Clearly label the uncertainty method instead of presenting all intervals as if they were calculated identically

### 3. Formalize Positivity and Weight Diagnostics

**Proposed file:** `analysis/07-positivity-weight-diagnostics.R`

**Goal**

Evaluate whether treated and control patients have adequate propensity-score overlap and whether IPW estimates depend on a small number of extreme weights.

**Why this is needed**

Good weighted covariate balance does not by itself establish positivity. Extreme propensity scores or weights can produce unstable estimates even when standardized mean differences are small.

**Expected outputs**

- `results/tables/positivity-diagnostics-basic.csv`
- `results/tables/ipw-weight-summary-basic.csv`
- `results/figures/propensity-score-overlap-basic.png`
- `results/figures/ipw-weight-distribution-basic.png`

**Required diagnostics**

- Propensity-score quantiles by treatment group
- Minimum and maximum propensity scores
- Percentage of observations affected by trimming
- Weight quantiles and maximum weight
- Effective sample size for ATE and ATT
- Number of observations with weak overlap

### 4. Add a Reproducible E-Value Analysis

**Proposed file:** `analysis/08-evalue-sensitivity.R`

**Goal**

Reproduce the existing E-value results from a formal analysis script and calculate E-values for both the point estimate and the confidence-limit estimate.

**Why this is needed**

The E-value results currently exist in `results/`, but there is no corresponding analysis script called by `main.R`. E-values provide a compact assessment of how strong an unmeasured confounder would need to be to explain away the estimated effect.

**Expected outputs**

- `results/tables/evalue-sensitivity-death30d.csv`
- `results/figures/evalue-sensitivity-death30d.png`

**Required reporting**

- Risk ratio and confidence interval used for the calculation
- Point-estimate E-value
- Confidence-limit E-value
- Plain-language interpretation

### 5. Add Partial R-Squared Sensitivity Analysis

**Proposed file:** `analysis/09-partial-r2-sensitivity.R`

**Goal**

Benchmark the strength of potential unmeasured confounding against the observed covariates using partial R-squared measures.

**Why this is needed**

E-values operate on the risk-ratio scale but do not directly compare an unmeasured confounder with observed confounders. Partial R-squared benchmarks help assess whether an omitted variable would need to be stronger than clinically important measured variables.

**Expected outputs**

- `results/tables/partial-r2-sensitivity-death30d.csv`
- `results/figures/partial-r2-sensitivity-death30d.png`

**Required reporting**

- Partial R-squared of each benchmark covariate with treatment
- Partial R-squared of each benchmark covariate with the outcome
- Treatment-outcome partial R-squared
- Identification of the strongest observed benchmark covariates

**Interpretation limitation**

This analysis will likely use a linear probability model approximation. The report must state that this is a sensitivity benchmark rather than the primary binary-outcome estimator.

### 6. Analyze the Secondary Outcome

**Proposed file:** `analysis/10-secondary-outcome-death180d.R`

**Goal**

Estimate the effect of RHC on 180-day mortality using the same primary adjustment set and the main estimators.

**Why this is needed**

The current report contains odds-ratio estimates for 180-day mortality but does not provide a complete risk-difference analysis with uncertainty. Comparing 30-day and 180-day effects clarifies whether the estimated harm is concentrated in the short term or persists over time.

**Expected outputs**

- `results/tables/secondary-effect-estimates-death180d.csv`
- `results/figures/secondary-effect-estimates-death180d.png`

**Required reporting**

- ATE and ATT estimates
- Risk differences and confidence intervals
- A comparison with the corresponding 30-day estimates

## Priority 1: Strong Robustness Checks

### 7. Evaluate Propensity-Score and Outcome-Model Specification

**Proposed file:** `analysis/11-model-specification-sensitivity.R`

**Goal**

Compare the primary estimates under alternative nuisance-model specifications.

**Why this is needed**

The current propensity-score and outcome models assume linear covariate effects and a specific interaction structure. Misspecification can bias regression adjustment and IPW, while AIPW is only doubly robust if at least one nuisance model is adequately specified.

**Expected outputs**

- `results/tables/model-specification-sensitivity-death30d.csv`
- `results/figures/model-specification-sensitivity-death30d.png`

**Recommended specifications**

- Current linear main-effects propensity-score model
- Restricted cubic splines or natural splines for age, APACHE score, PaO2/FIO2, mean blood pressure, and coma score
- Clinically defensible interactions
- A reduced model and a more flexible model

**Required reporting**

- Estimate under each specification
- Covariate balance under each propensity-score specification
- Effective sample size and extreme-weight diagnostics

### 8. Evaluate Cross-Fitting Stability for AIPW

**Proposed file:** `analysis/12-aipw-crossfit-sensitivity.R`

**Goal**

Assess whether the AIPW result is sensitive to the number of folds and random fold assignment.

**Why this is needed**

Cross-fitted AIPW estimates can vary with sample splitting, especially with moderately sized samples or flexible nuisance models. The final report should demonstrate that the AIPW result is not an artifact of one random split.

**Expected outputs**

- `results/tables/aipw-crossfit-sensitivity-death30d.csv`
- `results/figures/aipw-crossfit-sensitivity-death30d.png`

**Recommended design**

- Compare `K = 2`, `5`, and `10`
- Repeat each setting across multiple seeds
- Report the distribution, mean, and range of estimates

### 9. Add Inference to the Subgroup Analysis

**Proposed file:** `analysis/13-subgroup-heterogeneity-inference.R`

**Goal**

Add confidence intervals and formal interaction tests to the existing exploratory subgroup estimates.

**Why this is needed**

The current subgroup analysis reports point estimates only. Large differences between subgroup estimates may reflect sampling variation, especially in smaller groups such as patients with DNR status.

**Expected outputs**

- `results/tables/subgroup-heterogeneity-death30d.csv`
- `results/figures/subgroup-heterogeneity-death30d.png`

**Required reporting**

- Subgroup-specific estimates and confidence intervals
- Sample size and treatment counts by subgroup
- Interaction-test p-values
- A clear exploratory label and a multiple-comparisons caution

### 10. Evaluate Matching-Design Choices

**Proposed file:** `analysis/14-matching-design-sensitivity.R`

**Goal**

Assess whether matching estimates change under different matching specifications.

**Why this is needed**

The current matching analysis uses one-nearest-neighbor propensity-score matching with replacement. Matching results can depend strongly on the caliper, replacement rule, number of matches, and common-support restriction.

**Expected outputs**

- `results/tables/matching-design-sensitivity-death30d.csv`
- `results/tables/matching-balance-diagnostics-death30d.csv`
- `results/figures/matching-design-sensitivity-death30d.png`
- `results/figures/matching-balance-diagnostics-death30d.png`

**Recommended specifications**

- Matching with and without replacement
- One versus multiple nearest neighbors
- Several defensible caliper widths
- Common-support restriction

**Required reporting**

- Treatment-effect estimate
- Matched sample size
- Number of discarded observations
- Post-matching covariate balance

### 11. Compare Weighting Schemes

**Proposed estimator file:** `estimators/ate/overlap-weighting.R`

**Proposed analysis file:** `analysis/15-weighting-scheme-sensitivity.R`

**Goal**

Compare Hajek IPW with stabilized IPW and overlap weighting.

**Why this is needed**

Alternative weighting schemes target different populations and respond differently to weak overlap. Overlap weighting can provide a useful robustness check when extreme propensity scores are a concern.

**Expected outputs**

- `results/tables/weighting-scheme-sensitivity-death30d.csv`
- `results/tables/weighting-scheme-balance-death30d.csv`
- `results/figures/weighting-scheme-sensitivity-death30d.png`

**Required reporting**

- Target estimand for each weighting method
- Effect estimate and confidence interval
- Effective sample size
- Maximum weight and post-weighting balance

## Priority 2: Optional Extensions

### 12. Rosenbaum-Bounds Sensitivity for the Matched Analysis

**Proposed file:** `analysis/16-rosenbaum-bounds-matching.R`

**Goal**

Assess how strong hidden bias would need to be to change the matched-analysis conclusion.

**Why this is needed**

Rosenbaum bounds provide an unmeasured-confounding sensitivity analysis tailored to matched observational studies and complement the E-value analysis.

**Expected outputs**

- `results/tables/rosenbaum-bounds-death30d.csv`
- `results/figures/rosenbaum-bounds-death30d.png`

**Implementation condition**

Only perform this analysis after selecting and validating a final matched design. The interpretation must be tied specifically to the matched estimand and test statistic.

### 13. Consolidate Results Across Effect Scales

**Proposed file:** `analysis/17-effect-scale-comparison.R`

**Goal**

Present risk differences, risk ratios, and odds ratios together for the main outcomes.

**Why this is needed**

Odds ratios can exaggerate the apparent magnitude of common outcomes and are less directly interpretable than risk differences. A consolidated table improves clinical interpretation and prevents conclusions from depending on one effect scale.

**Expected outputs**

- `results/tables/effect-scale-comparison.csv`
- `results/figures/effect-scale-comparison.png`

**Required reporting**

- Explicit effect-scale labels
- Confidence intervals
- Separate panels or columns for 30-day and 180-day mortality

## Pipeline and Documentation Tasks

### 14. Update the Main Reproduction Pipeline

**File to update:** `main.R`

**Goal**

Ensure that one command regenerates every table and figure used in the final report.

**Why this is needed**

Several current outputs, including the primary IPW and E-value tables, do not have corresponding scripts called by `main.R`. Results that cannot be regenerated should not be treated as final.

**Expected outputs**

- A successful end-to-end `Rscript main.R` run
- All final tables and figures regenerated in `results/`

**Required changes**

- Add `analysis/00-estimator-comparison.R` to the current pipeline until it is replaced by `analysis/06-primary-effect-estimates.R`.
- Add all completed Priority 0 scripts.
- Use the shared analysis specification.
- Set seeds explicitly for every stochastic analysis.

### 15. Add a Result Manifest

**Proposed file:** `scripts/07-build-result-manifest.R`

**Goal**

Create an inventory linking every final table and figure to its generating script.

**Why this is needed**

A result manifest makes it easy to verify that the report only cites reproducible outputs and helps identify stale files.

**Expected output**

- `results/result-manifest.csv`

**Required fields**

- Result file
- Generating script
- Analysis purpose
- Primary or secondary status
- Date generated

## Recommended Execution Order

1. `scripts/06-define-analysis-specification.R`
2. Complete reusable outcome-regression and ATT AIPW estimators.
3. `analysis/06-primary-effect-estimates.R`
4. `analysis/07-positivity-weight-diagnostics.R`
5. `analysis/08-evalue-sensitivity.R`
6. `analysis/09-partial-r2-sensitivity.R`
7. `analysis/10-secondary-outcome-death180d.R`
8. Update and verify `main.R`.
9. Complete Priority 1 robustness checks.
10. Add optional Priority 2 analyses only if time permits.

## Final Report Minimum

The report should not be considered complete until it contains:

- A clearly defined primary estimand and outcome
- A DAG and justified adjustment set
- Descriptive statistics and raw imbalance
- Positivity and weighted-balance diagnostics
- A primary estimate with a confidence interval
- A doubly robust AIPW comparison
- At least one modeling robustness check
- At least one unmeasured-confounding sensitivity analysis
- A secondary-outcome analysis
- Clearly labeled exploratory subgroup results
- A fully reproducible `main.R` pipeline
