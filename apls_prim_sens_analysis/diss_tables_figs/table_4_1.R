library(dplyr)
library(gtsummary)
library(readr)

dat <- read_csv("apls_prim_sens_analysis/input_data/tmle_ready/tmle_primary_1y_labs.csv") |>
  mutate(
    treatment_group = factor(A, levels = c(0, 1), labels = c("Warfarin", "DOAC")),
    gender = factor(gender),
    race_clean = factor(race_clean),
    ethnicity_clean = factor(ethnicity_clean)
  )

dat <- dat[dat$death_before_index_flag == 0, ]


tbl_demog <- dat |>
  select(treatment_group, age_at_index, gender, race_clean, ethnicity_clean) |>
  tbl_summary(
    by = treatment_group,
    statistic = list(
      age_at_index ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany"
  ) |>
  add_overall() |>
  add_p() |>
  bold_labels() |>
  modify_caption("**Table 4.1. Baseline demographic characteristics by treatment group**")

tbl_demog

############################################
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(gtsummary)
  library(gt)
  library(scales)
  library(glue)
})

# ============================
# Config
# ============================
IN_FILE <- "apls_prim_sens_analysis/input_data/tmle_ready/tmle_primary_1y_final_dataset.csv"
OUT_DIR <- "apls_prim_sens_analysis/output_data/final_dissertation_data/diss_tables_figs/"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DATASET_LABEL <- "primary_1y_labs"

NA_STRINGS <- c("NaN", "nan", "NA", "na", "null", "NULL", "", " ")

COMORB_CANDIDATES <- c(
  "htn","dm","dyslipidemia","obesity","tobacco",
  "stroke_tia","mi_cad","pad","prior_vte",
  "ckd","liver_dz","cancer","copd","heart_failure",
  "atrial_fib","valvular","sle","prior_major_bleed",
  "pud","anemia"
)

label_comorb <- list(
  htn = "Hypertension (I10–I15)",
  dm = "Diabetes mellitus (E10–E13)",
  dyslipidemia = "Dyslipidemia (E78)",
  obesity = "Obesity (E66)",
  tobacco = "Tobacco use (F17, Z72.0, Z87.891)",
  stroke_tia = "Stroke/TIA (I63–I64, G45)",
  mi_cad = "MI/CAD (I21–I22, I25)",
  pad = "Peripheral artery disease (I70.2, I73.9, I74)",
  prior_vte = "Prior VTE (I26, I80–I82)",
  ckd = "Chronic kidney disease (N18–N19, Z99.2)",
  liver_dz = "Liver disease (K70–K74, K76)",
  cancer = "Cancer (C00–C97)",
  copd = "COPD (J44)",
  heart_failure = "Heart failure (I50)",
  atrial_fib = "Atrial fibrillation (I48)",
  valvular = "Valvular disease (I34–I37)",
  sle = "Systemic lupus erythematosus (M32)",
  prior_major_bleed = "Prior major bleed (I60–I62, K92.2)",
  pud = "Peptic ulcer disease (K25–K28)",
  anemia = "Anemia (D50–D53, D64)"
)

# ============================
# Helpers
# ============================
safe_factor <- function(x, unknown_level = "Unknown") {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- unknown_level
  factor(x)
}

make_age_group <- function(age) {
  age <- suppressWarnings(as.numeric(age))
  cut(
    age,
    breaks = c(-Inf, 39, 49, 59, 69, Inf),
    labels = c("<40", "40–49", "50–59", "60–69", "≥70"),
    right = TRUE
  )
}

style_gt <- function(gt_obj) {
  gt_obj %>%
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#F2F2F2"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_column_labels(everything())
    ) %>%
    gt::opt_row_striping()
}

# ============================
# Load data
# ============================
if (!file.exists(IN_FILE)) stop("Missing input file: ", IN_FILE)

apls <- fread(IN_FILE, na.strings = NA_STRINGS, showProgress = FALSE)
setDT(apls)

apls <- apls[apls$death_before_index_flag == 0, ]


# Standardize A/Y if needed
if (!"A" %in% names(apls) && "treatment_binary" %in% names(apls)) {
  apls[, A := as.integer(treatment_binary)]
}
if (!"Y" %in% names(apls) && "event_flag" %in% names(apls)) {
  apls[, Y := as.integer(event_flag)]
}

stopifnot(all(c("A", "Y") %in% names(apls)))

# Treatment labels
apls[, treatment_group := factor(
  ifelse(A == 1L, "DOAC", "Warfarin"),
  levels = c("Warfarin", "DOAC")
)]

# Age group
if ("age_at_index" %in% names(apls)) {
  apls[, age_group := make_age_group(age_at_index)]
}

# Key categoricals
if ("gender" %in% names(apls))          apls[, gender := safe_factor(gender)]
if ("race_clean" %in% names(apls))      apls[, race_clean := safe_factor(race_clean)]
if ("ethnicity_clean" %in% names(apls)) apls[, ethnicity_clean := safe_factor(ethnicity_clean)]

# Comorbidity vars present in this file
COMORB_VARS <- intersect(COMORB_CANDIDATES, names(apls))

# Convert comorbidities to Absence/Presence
for (v in COMORB_VARS) {
  apls[, (v) := suppressWarnings(as.integer(as.numeric(get(v))))]
  apls[is.na(get(v)), (v) := 0L]
  apls[, (v) := factor(
    ifelse(get(v) == 1L, "Presence", "Absence"),
    levels = c("Absence", "Presence")
  )]
}

# ============================
# Table 1: Demographics
# ============================
demo_vars <- c("age_group", "gender", "race_clean", "ethnicity_clean")
demo_vars <- intersect(demo_vars, names(apls))

if (length(demo_vars) > 0) {
  demo_label <- list(
    age_group ~ "Age at index (years)",
    gender ~ "Gender",
    race_clean ~ "Race",
    ethnicity_clean ~ "Ethnicity"
  )
  
  tbl1 <- apls %>%
    select(treatment_group, any_of(demo_vars)) %>%
    tbl_summary(
      by = treatment_group,
      statistic = all_categorical() ~ "{n} ({p})",
      digits = all_categorical() ~ c(0, 2),
      missing = "no",
      percent = "column",
      label = demo_label
    ) %>%
    add_overall(last = TRUE) %>%
    #add_p(test = all_categorical() ~ "chisq.test") %>%
    add_stat_label(label = all_categorical() ~ "N (%)") %>%
    bold_labels() %>%
    modify_header(label ~ "**Baseline characteristic**",
                  all_stat_cols() ~ "**{level}**") %>%
    modify_spanning_header(all_stat_cols() ~ NA) %>%
    modify_caption("**Table 1. Demographics (primary 1-year, labs file)**") %>%
    as_gt() %>%
    style_gt()
  
  gt::gtsave(tbl1, file.path(OUT_DIR, paste0(DATASET_LABEL, "_Table1_Demographics.html")))
  tryCatch(
    gt::gtsave(tbl1, file.path(OUT_DIR, paste0(DATASET_LABEL, "_Table1_Demographics.docx"))),
    error = function(e) message("DOCX export skipped for Table 1: ", e$message)
  )
} else {
  message("Skipping Table 1 (no demographics columns found).")
}

# ============================
# Table 2: Comorbidities
# ============================
if (length(COMORB_VARS) > 0) {
  tbl2 <- apls %>%
    select(treatment_group, all_of(COMORB_VARS)) %>%
    tbl_summary(
      by = treatment_group,
      statistic = all_categorical() ~ "{n} ({p}%)",
      missing = "no",
      percent = "column",
      label = label_comorb
    ) %>%
    add_overall() %>%
    add_p(test = all_categorical() ~ "fisher.test") %>%
    bold_labels() %>%
    modify_header(label ~ "**Comorbidities**") %>%
    modify_footnote(label ~ "Designations in parenthesis are ICD-10 codes.") %>%
    modify_caption("**Table 2. Comorbidities (primary 1-year, labs file)**") %>%
    as_gt() %>%
    style_gt()
  
  gt::gtsave(tbl2, file.path(OUT_DIR, paste0(DATASET_LABEL, "_Table2_Comorbidities.html")))
  tryCatch(
    gt::gtsave(tbl2, file.path(OUT_DIR, paste0(DATASET_LABEL, "_Table2_Comorbidities.docx"))),
    error = function(e) message("DOCX export skipped for Table 2: ", e$message)
  )
} else {
  message("Skipping Table 2 (no comorbidity columns found).")
}

# ============================
# Table 3: Mortality risk + rate
# ============================
if (!"pt_years" %in% names(apls)) {
  if ("time_to_event_days" %in% names(apls)) {
    apls[, pt_years := suppressWarnings(as.numeric(time_to_event_days) / 365.25)]
  } else {
    apls[, pt_years := NA_real_]
  }
}

apls[, pt_years := suppressWarnings(as.numeric(pt_years))]
apls[is.na(pt_years), pt_years := 0]

tbl3_dt <- apls[, .(
  N = .N,
  deaths = sum(Y == 1L, na.rm = TRUE),
  pt_years = sum(pt_years, na.rm = TRUE),
  risk = mean(Y == 1L, na.rm = TRUE)
), by = treatment_group]

tbl3_overall <- apls[, .(
  treatment_group = factor("Overall", levels = c("Overall", "Warfarin", "DOAC")),
  N = .N,
  deaths = sum(Y == 1L, na.rm = TRUE),
  pt_years = sum(pt_years, na.rm = TRUE),
  risk = mean(Y == 1L, na.rm = TRUE)
)]

tbl3_dt[, treatment_group := factor(as.character(treatment_group), levels = c("Warfarin", "DOAC"))]
tbl3_dt <- rbindlist(list(tbl3_overall, tbl3_dt), use.names = TRUE, fill = TRUE)

tbl3_dt[, rate_per_1000_py := fifelse(pt_years > 0, (deaths / pt_years) * 1000, NA_real_)]

tbl3_out <- tbl3_dt[, .(
  Treatment = as.character(treatment_group),
  N,
  `Mortality risk` = percent(risk, accuracy = 0.1),
  Deaths = deaths,
  `Patient-years` = round(pt_years, 1),
  `Mortality rate per 1,000 PY` = round(rate_per_1000_py, 2)
)]

gt_tbl3 <- tbl3_out %>%
  as.data.frame() %>%
  gt() %>%
  tab_caption("**Table 3. Mortality risk and rates (primary 1-year, labs file)**") %>%
  style_gt()

gt::gtsave(gt_tbl3, file.path(OUT_DIR, paste0(DATASET_LABEL, "_Table3_MortalityRisk_Rates.html")))
tryCatch(
  gt::gtsave(gt_tbl3, file.path(OUT_DIR, paste0(DATASET_LABEL, "_Table3_MortalityRisk_Rates.docx"))),
  error = function(e) message("DOCX export skipped for Table 3: ", e$message)
)

cat("\nSaved tables to: ", OUT_DIR, "\n", sep = "")
