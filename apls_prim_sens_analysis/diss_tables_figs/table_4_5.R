library(readr)
library(dplyr)
library(gt)
library(scales)

fmt_ci <- function(est, lo, hi, digits = 3) {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
    est, lo, hi
  )
}

res <- read_csv("apls_prim_sens_analysis/output_data/final_dissertation_data/primary_v1/tmle_primary_1y_summary.csv") 

tbl_tmle <- res |>
  transmute(
    `Propensity score truncation (gbound)` = gbound,
    `Counterfactual risk under DOAC` = percent(risk_doac, accuracy = 0.1),
    `Counterfactual risk under warfarin` = percent(risk_warfarin, accuracy = 0.1),
    `Risk difference (95% CI)` = fmt_ci(RD, RD_ci_l, RD_ci_u),
    `Risk ratio (95% CI)` = fmt_ci(RR, RR_ci_l, RR_ci_u)
  ) |>
  gt() |>
  tab_header(
    title = md("**Table 4.4. Primary 1-year TMLE estimates using the laboratory-adjusted covariate set**"),
    subtitle = "Fixed-horizon TMLE analysis of 1-year all-cause mortality"
  ) |>
  cols_align(
    align = "center",
    columns = everything()
  ) |>
  tab_source_note(
    source_note = md(
      "TMLE = targeted maximum likelihood estimation; CI = confidence interval."
    )
  ) |>
  tab_source_note(
    source_note = md(
      "The adjustment set included demographic characteristics, comorbidities, APS phenotype indicators, and laboratory covariates measured prior to treatment initiation."
    )
  )

tbl_tmle
