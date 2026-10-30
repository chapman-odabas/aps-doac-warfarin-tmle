library(readr)
library(dplyr)
library(gt)

fmt_ci <- function(est, lo, hi, digits = 3) {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
    est, lo, hi
  )
}

# -----------------------------
# Read files
# -----------------------------


res_sens <- read_csv(
  "apls_prim_sens_analysis/output_data/final_dissertation_data/sensitivity_analysis_v1/tmle_sensitivity_summaries.csv"
)


# -----------------------------
# 2) Remaining rows from sensitivity file
#    (exclude sens_1y_labs from this table)
# -----------------------------
rows_sens <- res_sens %>%
  filter(gbound == 0.01, analysis != "sens_1y_labs") %>%
  mutate(
    `Follow-up horizon` = case_when(
      analysis %in% c("sens_1y") ~ "1 year",
      analysis %in% c("sens_2y", "sens_2y_labs") ~ "2 years",
      analysis %in% c("sens_3y", "sens_3y_labs") ~ "3 years",
      analysis %in% c("sens_4y", "sens_4y_labs") ~ "4 years",
      analysis %in% c("sens_5y", "sens_5y_labs") ~ "5 years",
      TRUE ~ NA_character_
    ),
    `Covariate specification` = case_when(
      analysis %in% c("sens_2y_labs", "sens_3y_labs", "sens_4y_labs", "sens_5y_labs") ~
        "Demographic, clinical, and laboratory covariates",
      analysis %in% c("sens_1y", "sens_2y", "sens_3y", "sens_4y", "sens_5y") ~
        "Demographic and clinical covariates only",
      TRUE ~ NA_character_
    ),
    `Risk difference (95% CI)` = fmt_ci(RD, RD_ci_l, RD_ci_u),
    `Risk ratio (95% CI)` = fmt_ci(RR, RR_ci_l, RR_ci_u, digits = 2)
  ) %>%
  select(
    `Follow-up horizon`,
    `Covariate specification`,
    `Risk difference (95% CI)`,
    `Risk ratio (95% CI)`
  )

# -----------------------------
# Combine and order row_primary_1y, 
# -----------------------------
table_4_7_data <- bind_rows(rows_sens) %>%
  mutate(
    `Follow-up horizon` = factor(
      `Follow-up horizon`,
      levels = c("1 year", "2 years", "3 years", "4 years", "5 years")
    ),
    `Covariate specification` = factor(
      `Covariate specification`,
      levels = c(
        "Demographic and clinical covariates only",
        "Demographic, clinical, and laboratory covariates"
      )
    )
  ) %>%
  arrange(`Follow-up horizon`, `Covariate specification`)

# -----------------------------
# Build table
# -----------------------------
tbl_4_7 <- table_4_7_data %>%
  gt() %>%
  tab_header(
    title = md("**Table 4.7. Summary of sensitivity analyses across follow-up horizons and covariate specifications**"),
    subtitle = "TMLE estimates using propensity score truncation level gbound = 0.01"
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_options(
    table.font.size = px(13),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(12),
    data_row.padding = px(6)
  ) %>%
  tab_source_note(
    source_note = md("TMLE = targeted maximum likelihood estimation; CI = confidence interval.")
  ) %>%
  tab_source_note(
    source_note = md("The 1-year row reports the alternative model specification excluding laboratory covariates. Longer-horizon rows report sensitivity analyses with and without laboratory covariates.")
  )

tbl_4_7
