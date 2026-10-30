# DOAC versus warfarin initiation and one-year mortality in antiphospholipid syndrome

Analysis code for:

> Odabas O, Sheth P, Ehwerhemuepha L, Yaghmaei E, Rakovski C. *Direct oral anticoagulant versus warfarin initiation and one-year mortality in antiphospholipid syndrome: a target trial emulation using electronic health record data.* PLOS Medicine (under review).

The study emulates a target trial in de-identified electronic health record data and estimates the effect of initiating a direct oral anticoagulant (DOAC) versus warfarin on all-cause mortality using targeted maximum likelihood estimation (TMLE) with Super Learner.

## Data access

The analysis uses de-identified patient-level data from the Oracle Health Real-World Data network (formerly Cerner Real-World Data), licensed under a data use agreement that prohibits redistribution. **No data are included in this repository.** Researchers can request access from Oracle Health (https://www.oracle.com/health/). The expected input files and their columns are described in [`docs/data_dictionary.md`](docs/data_dictionary.md).

## Repository structure

```
apls_prim_sens_analysis/
├── 01_preprocess_tmle.R                  # cleans horizon-specific cohort files; derives covariates
├── 02_tmle_primary_1y.R                  # primary analysis: 1-year mortality, laboratory-adjusted
├── 03_tmle_sensitivity_2y_3y_4y_5y.R     # 1-year without laboratory covariates; 2- to 5-year horizons
├── diss_tables_figs/                     # scripts that format manuscript tables and figures
├── input_data/raw_data/                  # place licensed input files here (not tracked)
├── input_data/tmle_ready/                # written by 01_preprocess_tmle.R (not tracked)
└── output_data/                          # results written by the analysis scripts (not tracked)
docs/data_dictionary.md
```

## How to run

Run all scripts from the repository root (paths are relative to it), in this order:

1. `Rscript apls_prim_sens_analysis/01_preprocess_tmle.R`
2. `Rscript apls_prim_sens_analysis/02_tmle_primary_1y.R`
3. `Rscript apls_prim_sens_analysis/03_tmle_sensitivity_2y_3y_4y_5y.R`
4. Scripts in `apls_prim_sens_analysis/diss_tables_figs/` (any order)

Random seeds are fixed (`set.seed(42)`) for imputation and estimation.

## Scripts and manuscript items

| Script | Produces | Manuscript item |
| --- | --- | --- |
| `01_preprocess_tmle.R` | Analysis-ready files for each horizon, with and without laboratory values | All analyses |
| `02_tmle_primary_1y.R` | Primary TMLE estimates at gbound 0.025 and 0.010; crude risks | Table 4; Table 5 (1-year, with laboratory) |
| `03_tmle_sensitivity_2y_3y_4y_5y.R` | TMLE estimates without laboratory covariates (1 year) and at 2–5 years, with E-values | Table 5; S4 Table |
| `diss_tables_figs/figure_4_1.R` | Cohort flow diagram | Fig 2 |
| `diss_tables_figs/table_4_1.R` | Baseline characteristics | Tables 1–2 |
| `diss_tables_figs/table_4_3.R` | Unadjusted mortality risks and rates | Table 3 |
| `diss_tables_figs/table_4_4.R`, `table_4_5.R` | Formatted primary TMLE results | Table 4 |
| `diss_tables_figs/table_4_7.R` | Formatted horizon results | S4 Table |
| `diss_tables_figs/figure_4_2.R` | Propensity score distributions (main-effects logistic regression, diagnostic only) | S2 Fig |
| `diss_tables_figs/fig_4_3.R` | Covariate differences before and after weighting | S3 Fig |

## Key analytic settings

- Treatment: A = 1 for DOAC, 0 for warfarin, by first recorded class after the APS diagnosis date.
- Outcome: Y = 1 for recorded death within the horizon.
- Nuisance models: Super Learner (10-fold cross-validation, `method.NNLS`) with `SL.glmnet`, `SL.mean` and `SL.ranger`; full-sample predictions supplied to `tmle()` (no cross-fitting).
- Propensity score truncation: `gbound` = 0.025 (primary) and 0.010.
- Missing INR and creatinine: single imputation with `mice` (predictive mean matching, 10 iterations, seed 42). Platelet count is excluded (about 97% missing).

## Not included

- **Cohort extraction.** Construction of the horizon-specific cohort files from the Oracle Health database was performed in PySpark within the Oracle Health environment. [TO ADD if available]
- **Super Learner benchmarking.** The comparison of 63 library configurations summarized in S1 Fig. [TO ADD if available]

## Software

R 4.5.3 with packages: data.table, stringr, dplyr, readr, mice, SuperLearner, glmnet, ranger, tmle, EValue, cobalt, ggplot2, gt, gtsummary, scales, glue, DiagrammeR.

## License

MIT (see `LICENSE`). If you use this code, please cite the article above.
