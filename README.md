# STAT566 Causal Modeling Final Project

## Project Structure

### Overview:
```
├── analysis/
├── data/
│   ├── processed/
│   └── raw/
├── estimators/
│   ├── ate/
│   ├── att/
│   └── cate/
├── notebooks/
├── results/
│   ├── figures/
│   └── tables/
└── scripts/
```
### Folders
- `scripts`: reproducible data constructing scripts
- `analysis`: reproducible analysis results scripts
- `data/`
    - `raw`: raw data location
    - `processed`: processed data location
- `results/` : analysis results location
- `estimators/`: estimator functions definitions and utils by types
- `notebooks/`: analysis Rmd, pdf notebook location

## RHC data



## Progress

- May 12, 2026
    - progress report
        - SMD table, love plot, propensity score distribution
        - new DAG on overleaf


- May 19, 2026
    - Estimator implementation
        - ATT (Andy)
        - AIPW (Jian)
        - CATE, matching (LC)
