# STAT566 Causal Modeling Final Project

## Project Location

- Github Repo: [Link](https://github.com/liangchenn/stat566-final-proj)

- Overleaf Project: [Link](https://www.overleaf.com/read/yqswbkvkcgqz#f30bf9)



## Project Structure

### Usage

- Execute `main.R` will: 
    - create raw data in `data/raw/` from relevant packages
    - create analysis dataset `data/processed/` from raw data on the Internet
    - generate all figures and tables we put in the report in `results/{figures, tables}/`
    
- Or run `$ make run` in the terminal to regenerate all the results we use in the final report.
- `$ make clean` to clean up all results, processed data and intermediate files.

- `$ make`: to reproduce all results.

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

- `final-report.pdf`: The final report
- `presentation-slides.pdf`: final presentation slides

- `scripts`: reproducible data constructing scripts
- `analysis`: reproducible analysis results scripts
- `data/`
    - `raw`: raw data location
    - `processed`: processed data location
- `results/` : analysis results location
- `estimators/`: estimator functions definitions and utils by types
- `notebooks/`: legacy analysis Rmd, pdf notebook location
- `results/`: all the tables and figures we used in the final report

## RHC data

- CHECK the `scripts/02-defined-variables.R` for the outcomes, treatment, and different adjustment sets
