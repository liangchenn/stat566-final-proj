# Project Progress

Updated: 2026-06-03

## Current Study Setup

- Research question: Estimate the causal effect of receiving right heart catheterization (RHC) on mortality risk.
- Analysis sample: 5,735 patients, including 2,184 who received RHC and 3,551 who did not receive RHC.
- Primary outcome: 30-day mortality; 180-day mortality is also available.
- Main estimands: ATE and ATT.
- Primary adjustment set: The basic adjustment set, containing 16 demographic, comorbidity, and day-1 physiological covariates.
- A DAG representing the causal assumptions has been created: `results/figures/rhc-dag-ggdag.png`.

## Completed Analyses and Results

### 1. Descriptive Statistics and Unadjusted Group Differences

- The 30-day mortality rate was 38.0% in the RHC group and 30.6% in the no-RHC group, corresponding to an unadjusted risk difference of 7.36 percentage points.
- The 180-day mortality rate was 68.0% in the RHC group and 63.0% in the no-RHC group.
- Patients who received RHC were generally more severely ill. For example:
  - APACHE score: 60.7 vs. 50.9.
  - PaO2/FIO2 ratio: 192.4 vs. 240.6.
  - Mean blood pressure: 68.2 vs. 84.9.
- In the unweighted data, 7 of the 16 basic covariates had an absolute SMD greater than 0.1. The largest imbalances were observed for APACHE score (0.501), mean blood pressure (0.455), and PaO2/FIO2 ratio (0.433).

Related outputs:

- `results/tables/descriptive-stat-basic-by-rhc.csv`
- `results/tables/descriptive-stat-table.csv`
- `results/tables/raw-smd-basic.csv`
- `results/figures/raw-smd-love-plot.png`

### 2. Covariate-Balance Diagnostics After Propensity-Score Weighting

- After ATE IPW and ATT IPW, all basic covariates had absolute SMDs below 0.1.
- The maximum absolute SMD was 0.020 after ATE IPW and 0.048 after ATT IPW.
- These results indicate that the basic propensity-score model substantially improved balance in the observed covariates.

Related outputs:

- `results/tables/balance-diagnostics-basic.csv`
- `results/figures/balance-diagnostics-basic.png`
- `results/figures/propensity-score-overlap-basic.png`

### 3. Primary IPW Causal-Effect Estimates

Using 30-day mortality as the outcome, the basic adjustment set, and Hajek IPW:

| Estimand | RHC risk | No-RHC counterfactual risk | Risk difference | 95% bootstrap CI | Risk ratio | 95% bootstrap CI |
|---|---:|---:|---:|---:|---:|---:|
| ATE | 37.31% | 31.67% | 5.64 pp | [2.91, 7.94] pp | 1.178 | [1.090, 1.256] |
| ATT | 38.00% | 33.32% | 4.68 pp | [1.94, 7.32] pp | 1.140 | [1.055, 1.228] |

- The primary estimates indicate that receiving RHC was associated with a higher risk of 30-day mortality. The bootstrap confidence intervals for both the ATE and ATT excluded zero.
- The weighted effective sample sizes were approximately 4,411 for the ATE and 4,020 for the ATT.

Related output:

- `results/tables/primary-ipw-estimate-death30d.csv`

### 4. Comparison Across Estimation Methods

Estimated risk differences for 30-day mortality:

| Estimand | Regression adjustment | Hajek IPW | PS matching |
|---|---:|---:|---:|
| ATE | 6.15 pp | 5.64 pp | 5.15 pp |
| ATT | 6.19 pp | 4.68 pp | 4.07 pp |

- All adjusted estimates were smaller than the unadjusted risk difference of 7.36 percentage points, but their directions were consistent.
- Across the ATE, ATT, and different estimation methods, the estimates ranged from approximately 4.1 to 6.2 percentage points. This suggests that the main conclusion is not driven entirely by a single estimation method.

Related outputs:

- `results/tables/estimator-comparison-death30d.csv`
- `results/figures/estimator-comparison-death30d.png`

### 5. Odds-Ratio Results

- The unadjusted odds ratio for 30-day mortality was 1.388. Adjusted odds ratios were approximately:
  - ATE: 1.256 to 1.307.
  - ATT: 1.194 to 1.305.
- The unadjusted odds ratio for 180-day mortality was 1.252. Adjusted odds ratios were approximately:
  - ATE: 1.154 to 1.190.
  - ATT: 1.041 to 1.143.
- The estimated effects were generally larger for 30-day mortality than for 180-day mortality.

Related outputs:

- `results/tables/odds-ratio-estimator-comparison.csv`
- `results/figures/odds-ratio-estimator-comparison.png`

### 6. Robustness Checks

#### Propensity-Score Trimming

- Changing the trimming threshold from `[0.001, 0.999]` to `[0.100, 0.900]` changed the ATE estimate only from 5.64 to 5.59 percentage points and the ATT estimate from 4.68 to 4.66 percentage points.
- The results were highly stable across propensity-score trimming thresholds.

#### Adjustment-Set Sensitivity

- Estimates remained positive when using the minimal, basic, and detailed adjustment sets.
- Across adjustment sets and estimators, the estimated ranges were:
  - ATE: 3.79 to 8.04 percentage points.
  - ATT: 4.07 to 7.60 percentage points.
- Estimates using the minimal adjustment set were generally larger, suggesting that adjustment for demographic variables alone may leave residual confounding. Results using the basic and detailed adjustment sets were more similar.

Related outputs:

- `results/tables/trimming-sensitivity-death30d.csv`
- `results/figures/trimming-sensitivity-death30d.png`
- `results/tables/adjustment-set-sensitivity-death30d.csv`
- `results/figures/adjustment-set-sensitivity-death30d.png`

### 7. Sensitivity to Unmeasured Confounding

- For the ATE IPW risk ratio, the point-estimate E-value was 1.64 and the confidence-limit E-value was 1.40.
- For the ATT IPW risk ratio, the point-estimate E-value was 1.54 and the confidence-limit E-value was 1.29.
- In other words, after conditioning on the measured covariates, an unmeasured confounder would need to have risk-ratio associations of approximately 1.5 to 1.6 with both RHC receipt and 30-day mortality to fully explain away the point estimate.

Related outputs:

- `results/tables/evalue-sensitivity-death30d.csv`
- `results/figures/evalue-sensitivity-death30d.png`

### 8. Exploratory Subgroup Analysis

- All currently examined subgroup estimates were positive, but their magnitudes varied substantially.
- Larger ATE estimates were observed among:
  - Patients with DNR status: 21.1 percentage points.
  - The high-APACHE group: 10.5 percentage points.
  - The older group: 9.2 percentage points.
- In comparison, the ATE was 3.7 percentage points for patients without DNR status, 3.5 percentage points for the low-APACHE group, and 2.0 percentage points for the younger group.
- These results may be treated as exploratory evidence of treatment-effect heterogeneity. However, there are currently no subgroup-specific confidence intervals or formal interaction tests, so the results should not be interpreted as definitive evidence of effect modification.

Related outputs:

- `results/tables/cate-subgroup-death30d.csv`
- `results/figures/cate-subgroup-death30d.png`

## Methods Implemented but Not Yet Available as Formal Reportable Results

- AIPW:
  - `notebooks/AIPW.Rmd` contains a cross-fitted AIPW/DML implementation.
  - `notebooks/att.Rmd` contains ATT AIPW and bootstrap code.
  - There is currently no formal result table that consistently uses the processed data and basic adjustment set, and the analysis has not yet been added to `main.R`.
- Partial R-squared sensitivity benchmark:
  - `notebooks/simple-sensitivity-analysis.Rmd` contains the analysis code, but the results have not yet been exported to `results/`.
- ATT notebook PDF:
  - `notebooks/fp.pdf` contains an older analysis and warnings. Its numerical results should not be cited directly.

## Reproducibility Gaps

- `main.R` currently does not execute `analysis/00-estimator-comparison.R`, so the estimator-comparison table and figure are not regenerated by the main pipeline.
- `primary-ipw-estimate-death30d.csv` and the E-value results currently exist in `results/`, but the repository does not contain corresponding formal analysis scripts called by `main.R`.
- Before producing the next version of the formal report, AIPW, primary IPW uncertainty estimation, E-value analysis, and partial R-squared analysis should be incorporated into the reproducible pipeline.
