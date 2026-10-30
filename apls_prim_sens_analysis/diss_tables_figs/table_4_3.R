library(readr)
library(dplyr)
library(gt)
library(scales)

dat <- read_csv("apls_prim_sens_analysis/input_data/tmle_ready/tmle_primary_1y_final_dataset.csv") |>
  mutate(
    treatment_group = factor(A, levels = c(0, 1), labels = c("Warfarin", "DOAC"))
  )

dat <- dat[dat$death_before_index_flag == 0, ]


# treatment-specific summaries
by_trt <- dat |>
  group_by(treatment_group) |>
  summarise(
    `Number of patients` = n(),
    `Deaths within 1 year` = sum(Y, na.rm = TRUE),
    `Mortality risk` = mean(Y, na.rm = TRUE),
    `Total patient-years` = sum(pt_years, na.rm = TRUE),
    `Mortality rate per 1,000 patient-years` = 1000 * sum(Y, na.rm = TRUE) / sum(pt_years, na.rm = TRUE),
    .groups = "drop"
  )

# overall summary
overall <- dat |>
  summarise(
    treatment_group = "Overall",
    `Number of patients` = n(),
    `Deaths within 1 year` = sum(Y, na.rm = TRUE),
    `Mortality risk` = mean(Y, na.rm = TRUE),
    `Total patient-years` = sum(pt_years, na.rm = TRUE),
    `Mortality rate per 1,000 patient-years` = 1000 * sum(Y, na.rm = TRUE) / sum(pt_years, na.rm = TRUE)
  )

# combine and reshape for display
tbl_crude <- bind_rows(by_trt, overall) |>
  mutate(
    `Mortality risk` = percent(`Mortality risk`, accuracy = 0.1),
    `Total patient-years` = round(`Total patient-years`, 1),
    `Mortality rate per 1,000 patient-years` = round(`Mortality rate per 1,000 patient-years`, 2)
  ) |>
  rename(`Treatment group` = treatment_group) |>
  gt() |>
  tab_header(
    title = md("**Table 4.3. Unadjusted 1-year mortality risk and mortality rates by treatment group**"),
    subtitle = "Primary 1-year analytic cohort"
  ) |>
  cols_label(
    `Treatment group` = "Treatment Group",
    `Number of patients` = "N",
    `Deaths within 1 year` = "Deaths",
    `Mortality risk` = "Mortality Risk",
    `Total patient-years` = "Patient-Years",
    `Mortality rate per 1,000 patient-years` = "Mortality Rate per 1,000 PY"
  ) |>
  tab_source_note(
    source_note = md(
      "Mortality risk is the proportion of patients who died within 365 days of treatment initiation."
    )
  ) |>
  tab_source_note(
    source_note = md(
      "Mortality rate per 1,000 patient-years was calculated as deaths divided by total follow-up time, multiplied by 1,000."
    )
  ) |>
  cols_align(
    align = "center",
    columns = everything()
  )

tbl_crude
