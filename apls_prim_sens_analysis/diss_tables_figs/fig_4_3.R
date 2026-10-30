
library(readr)
library(dplyr)
library(cobalt)
library(ggplot2)

# ============================
# Load data
# ============================
dat <- read_csv("apls_prim_sens_analysis/input_data/tmle_ready/tmle_primary_1y_final_dataset.csv") %>%
  mutate(
    A = as.integer(A),
    gender = factor(gender),
    race_clean = factor(race_clean),
    ethnicity_clean = factor(ethnicity_clean)
  )

# Keep analytic cohort only
if ("death_before_index_flag" %in% names(dat)) {
  dat <- dat %>% filter(death_before_index_flag == 0)
}

dat <- dat %>% filter(A %in% c(0L, 1L))

# ============================
# Covariate specification
# ============================
bal_formula <- A ~ age_at_index + gender + race_clean + ethnicity_clean +
  htn + dm + dyslipidemia + obesity + tobacco +
  stroke_tia + mi_cad + pad + prior_vte + ckd + liver_dz +
  cancer + copd + heart_failure + atrial_fib + valvular + sle +
  prior_major_bleed + pud + anemia + b2gpi_pos + la_pos +
  arterial_thrombosis + venous_thrombosis + high_risk_phenotype

# ============================
# Balance diagnostics before adjustment
# ============================
bal_unadj <- bal.tab(
  bal_formula,
  data = dat,
  estimand = "ATE"
)

# ============================
# Pretty variable labels
# ============================
nice_labels <- c(
  age_at_index = "Age at treatment initiation",
  htn = "Hypertension",
  dm = "Diabetes mellitus",
  dyslipidemia = "Dyslipidemia",
  obesity = "Obesity",
  tobacco = "Tobacco use",
  stroke_tia = "Prior stroke or TIA",
  mi_cad = "Coronary artery disease or myocardial infarction",
  pad = "Peripheral artery disease",
  prior_vte = "History of venous thromboembolism",
  ckd = "Chronic kidney disease",
  liver_dz = "Liver disease",
  cancer = "Malignancy",
  copd = "Chronic obstructive pulmonary disease",
  heart_failure = "Heart failure",
  atrial_fib = "Atrial fibrillation",
  valvular = "Valvular heart disease",
  sle = "Systemic lupus erythematosus",
  prior_major_bleed = "Prior major bleeding",
  pud = "Peptic ulcer disease",
  anemia = "Anemia",
  b2gpi_pos = "Anti-β2 glycoprotein I antibody positive",
  la_pos = "Lupus anticoagulant positive",
  arterial_thrombosis = "Prior arterial thrombosis",
  venous_thrombosis = "Prior venous thrombosis",
  high_risk_phenotype = "High-risk APS phenotype",
  "gender_Female" = "Female sex",
  "gender_Male" = "Male sex",
  "gender_Unknown" = "Sex unknown",
  "race_clean_White" = "White race",
  "race_clean_Black or African American" = "Black or African American race",
  "race_clean_Asian" = "Asian race",
  "race_clean_American Indian or Alaska Native" = "American Indian or Alaska Native race",
  "race_clean_Native Hawaiian or Other Pacific Islander" = "Native Hawaiian or Pacific Islander race",
  "race_clean_Other" = "Other race",
  "race_clean_Unknown" = "Race unknown",
  "ethnicity_clean_Hispanic or Latino" = "Hispanic or Latino ethnicity",
  "ethnicity_clean_Non-Hispanic" = "Non-Hispanic ethnicity",
  "ethnicity_clean_Unknown" = "Ethnicity unknown"
)

# ============================
# Plot
# ============================
p_before <- love.plot(
  bal_unadj,
  stat = "mean.diffs",
  abs = TRUE,
  thresholds = c(m = 0.1),
  var.order = "unadjusted",
  var.names = nice_labels,
  sample.names = c("Before adjustment"),
  title = "Standardized Mean Differences of Baseline Covariates Before Adjustment",
  xlab = "Absolute standardized mean difference"
) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

print(p_before)

ggsave(
  "apls_prim_sens_analysis/output_data/final_dissertation_data/tables_figures/figure_5_smd_before_adjustment.png",
  plot = p_before,
  width = 10,
  height = 8,
  dpi = 600
)




######################################################
# After-weighting SMD plot script 
######################################################

library(readr)
library(dplyr)
library(cobalt)
library(ggplot2)

# ============================
# Load data
# ============================
dat <- read_csv("apls_prim_sens_analysis/input_data/tmle_ready/tmle_primary_1y_final_dataset.csv") %>%
  mutate(
    A = as.integer(A),
    gender = factor(gender),
    race_clean = factor(race_clean),
    ethnicity_clean = factor(ethnicity_clean)
  )

if ("death_before_index_flag" %in% names(dat)) {
  dat <- dat %>% filter(death_before_index_flag == 0)
}

dat <- dat %>% filter(A %in% c(0L, 1L))

# ============================
# Covariate specification
# ============================
bal_formula <- A ~ age_at_index + gender + race_clean + ethnicity_clean +
  htn + dm + dyslipidemia + obesity + tobacco +
  stroke_tia + mi_cad + pad + prior_vte + ckd + liver_dz +
  cancer + copd + heart_failure + atrial_fib + valvular + sle +
  prior_major_bleed + pud + anemia + b2gpi_pos + la_pos +
  arterial_thrombosis + venous_thrombosis + high_risk_phenotype

# ============================
# Fit propensity score model
# ============================
ps_model <- glm(
  bal_formula,
  data = dat,
  family = binomial()
)

dat <- dat %>%
  mutate(
    ps = predict(ps_model, type = "response"),
    ps = pmin(pmax(ps, 0.01), 0.99),   # aligns with gbound = 0.01
    wt_ate = if_else(A == 1L, 1 / ps, 1 / (1 - ps))
  )

# ============================
# Balance diagnostics after weighting
# ============================
bal_adj <- bal.tab(
  bal_formula,
  data = dat,
  weights = dat$wt_ate,
  method = "weighting",
  estimand = "ATE",
  un = TRUE
)

print(bal_adj)

# ============================
# Pretty variable labels
# ============================
nice_labels <- c(
  age_at_index = "Age at treatment initiation",
  htn = "Hypertension",
  dm = "Diabetes mellitus",
  dyslipidemia = "Dyslipidemia",
  obesity = "Obesity",
  tobacco = "Tobacco use",
  stroke_tia = "Prior stroke or TIA",
  mi_cad = "Coronary artery disease or myocardial infarction",
  pad = "Peripheral artery disease",
  prior_vte = "History of venous thromboembolism",
  ckd = "Chronic kidney disease",
  liver_dz = "Liver disease",
  cancer = "Malignancy",
  copd = "Chronic obstructive pulmonary disease",
  heart_failure = "Heart failure",
  atrial_fib = "Atrial fibrillation",
  valvular = "Valvular heart disease",
  sle = "Systemic lupus erythematosus",
  prior_major_bleed = "Prior major bleeding",
  pud = "Peptic ulcer disease",
  anemia = "Anemia",
  b2gpi_pos = "Anti-β2 glycoprotein I antibody positive",
  la_pos = "Lupus anticoagulant positive",
  arterial_thrombosis = "Prior arterial thrombosis",
  venous_thrombosis = "Prior venous thrombosis",
  high_risk_phenotype = "High-risk APS phenotype",
  "gender_Female" = "Female sex",
  "gender_Male" = "Male sex",
  "gender_Unknown" = "Sex unknown",
  "race_clean_White" = "White race",
  "race_clean_Black or African American" = "Black or African American race",
  "race_clean_Asian" = "Asian race",
  "race_clean_American Indian or Alaska Native" = "American Indian or Alaska Native race",
  "race_clean_Native Hawaiian or Other Pacific Islander" = "Native Hawaiian or Pacific Islander race",
  "race_clean_Other" = "Other race",
  "race_clean_Unknown" = "Race unknown",
  "ethnicity_clean_Hispanic or Latino" = "Hispanic or Latino ethnicity",
  "ethnicity_clean_Non-Hispanic" = "Non-Hispanic ethnicity",
  "ethnicity_clean_Unknown" = "Ethnicity unknown"
)

# ============================
# Adjusted-only plot
# ============================
p_after <- love.plot(
  bal_adj,
  stat = "mean.diffs",
  abs = TRUE,
  thresholds = c(m = 0.1),
  var.order = "unadjusted",
  var.names = nice_labels,
  sample.names = c("Before adjustment", "After weighting"),
  title = "Standardized Mean Differences of Baseline Covariates After Weighting",
  xlab = "Absolute standardized mean difference",
  which = "adjusted"
) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.85, 0.5)
  )

print(p_after)

ggsave(
  "apls_prim_sens_analysis/output_data/final_dissertation_data/tables_figures/figure_5_smd_before_after_weighting.png",
  plot = p_after,
  width = 10,
  height = 12,
  dpi = 600
)

