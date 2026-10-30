#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# ============================
# Config
# ============================
RAW_DIR  <- "apls_prim_sens_analysis/input_data/raw_data/"
OUT_DIR  <- "apls_prim_sens_analysis/input_data/tmle_ready"
#dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

INPUTS <- list(
  primary_1y      = file.path(RAW_DIR, "df_1y_df.csv"),
  primary_1y_labs = file.path(RAW_DIR, "df_1y_labs.csv"),
  sens_2y         = file.path(RAW_DIR, "df_2y_df.csv"),
  sens_2y_labs    = file.path(RAW_DIR, "df_2y_labs.csv"),
  sens_3y         = file.path(RAW_DIR, "df_3y_df.csv"),
  sens_3y_labs    = file.path(RAW_DIR, "df_3y_labs.csv"),
  sens_4y         = file.path(RAW_DIR, "df_4y_df.csv"),
  sens_4y_labs    = file.path(RAW_DIR, "df_4y_labs.csv"),
  sens_5y         = file.path(RAW_DIR, "df_5y_df.csv"),
  sens_5y_labs    = file.path(RAW_DIR, "df_5y_labs_df.csv")
)

OUTPUTS <- list(
  primary_1y      = file.path(OUT_DIR, "tmle_primary_1y.csv"),
  primary_1y_labs = file.path(OUT_DIR, "tmle_sens_1y_labs.csv"),
  sens_2y         = file.path(OUT_DIR, "tmle_sens_2y.csv"),
  sens_2y_labs    = file.path(OUT_DIR, "tmle_sens_2y_labs.csv"),
  sens_3y         = file.path(OUT_DIR, "tmle_sens_3y.csv"),
  sens_3y_labs    = file.path(OUT_DIR, "tmle_sens_3y_labs.csv"),
  sens_4y         = file.path(OUT_DIR, "tmle_sens_4y.csv"),
  sens_4y_labs    = file.path(OUT_DIR, "tmle_sens_4y_labs.csv"),
  sens_5y         = file.path(OUT_DIR, "tmle_sens_5y.csv"),
  sens_5y_labs    = file.path(OUT_DIR, "tmle_sens_5y_labs.csv")
)

NA_STRINGS <- c("NaN","nan","NA","na","null","NULL",""," ")

# Expected binary flags you commonly have
COMORB_VARS <- c(
  "htn","dm","dyslipidemia","obesity","tobacco",
  "stroke_tia","mi_cad","pad","prior_vte",
  "ckd","liver_dz","cancer","copd","heart_failure",
  "atrial_fib","valvular","sle","prior_major_bleed",
  "pud","anemia"
)

AB_FLAGS <- c("acl_pos","b2gpi_pos","la_pos","triple_positive")

LAB_VARS <- c(
  "inr_value","platelet_value","creatinine_value",
  "inr_date","platelet_date","creatinine_date"
)

DATE_VARS <- c("apls_first_date","index_date","first_anticoag_date","death_date","censor_date",
               "inr_date","platelet_date","creatinine_date")

# ============================
# Helpers
# ============================

clean_list_string <- function(x) {
  # Turns "['White', 'White']" -> "White"
  x <- as.character(x)
  x[x %in% NA_STRINGS] <- NA_character_
  x <- gsub("\\[|\\]|'", "", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x <- sub(",.*$", "", x)
  x[x == ""] <- NA_character_
  x
}

drop_constant_covariates <- function(dt, vars) {
  vars <- intersect(vars, names(dt))
  dropped <- character(0)
  
  for (v in vars) {
    x <- dt[[v]]
    # treat NA as its own category; but in your case these are 0/1 already
    ux <- unique(x)
    ux <- ux[!is.na(ux)]
    if (length(ux) <= 1) dropped <- c(dropped, v)
  }
  
  if (length(dropped) > 0) {
    dt[, (dropped) := NULL]
    cat("Dropped constant covariates:", paste(dropped, collapse = ", "), "\n")
  } else {
    cat("No constant covariates dropped.\n")
  }
  
  invisible(dropped)
}


parse_date_col <- function(dt, col) {
  if (!col %in% names(dt)) return(invisible(NULL))
  dt[, (col) := as.character(get(col))]
  dt[get(col) %in% NA_STRINGS, (col) := NA_character_]
  # Many of your files have YYYY-MM-DD; as.IDate handles that
  suppressWarnings(dt[, (col) := as.IDate(get(col))])
  invisible(NULL)
}

median_impute <- function(x) {
  x <- as.numeric(x)
  med <- suppressWarnings(median(x, na.rm = TRUE))
  if (!is.finite(med)) return(x)
  x[is.na(x)] <- med
  x
}

standardize_gender <- function(dt) {
  # If one-hot exists, use it
  if (all(c("gender_female","gender_male","gender_unknown") %in% names(dt))) {
    dt[, gender := fifelse(gender_female == 1, "Female",
                           fifelse(gender_male == 1, "Male", "Unknown"))]
    dt[, gender := factor(gender, levels = c("Female","Male","Unknown"))]
    return(invisible(NULL))
  }
  
  if (!"gender" %in% names(dt)) {
    dt[, gender := factor("Unknown", levels = c("Female","Male","Unknown"))]
    return(invisible(NULL))
  }
  
  g <- as.character(dt$gender)
  g[g %in% NA_STRINGS] <- NA_character_
  g <- str_trim(g)
  gl <- str_to_lower(g)
  
  # IMPORTANT: fcase conditions must all be length N; use logical vectors
  is_female <- !is.na(gl) & gl %in% c("female","f")
  is_male   <- !is.na(gl) & gl %in% c("male","m")
  is_unk    <- is.na(gl) | gl %in% c("other","unknown","other or unknown")
  
  out <- rep("Unknown", length(gl))
  out[is_female] <- "Female"
  out[is_male]   <- "Male"
  out[is_unk]    <- "Unknown"
  
  dt[, gender := factor(out, levels = c("Female","Male","Unknown"))]
  invisible(NULL)
}

coerce_binary_01 <- function(dt, vars, na_to_zero = TRUE) {
  vars <- intersect(vars, names(dt))
  for (v in vars) {
    dt[, (v) := suppressWarnings(as.numeric(as.character(get(v))))]
    # map weird values defensively
    if (na_to_zero) dt[is.na(get(v)), (v) := 0]
    dt[, (v) := as.integer(get(v) != 0)]
  }
  invisible(NULL)
}

fix_time_to_event <- function(dt) {
  if (!("time_to_event_days" %in% names(dt))) return(invisible(NULL))
  dt[, time_to_event_days := suppressWarnings(as.integer(as.numeric(time_to_event_days)))]
  dt[is.na(time_to_event_days), time_to_event_days := 0L]
  
  # If negative, set to 0 (usually date rounding / source weirdness)
  dt[time_to_event_days < 0, time_to_event_days := 0L]
  
  invisible(NULL)
}

derive_race_ethnicity <- function(dt) {
  if ("races" %in% names(dt) && !"race_clean" %in% names(dt)) {
    dt[, race_clean := clean_list_string(races)]
  }
  if ("ethnicities" %in% names(dt) && !"ethnicity_clean" %in% names(dt)) {
    dt[, ethnicity_clean := clean_list_string(ethnicities)]
  }
  if ("race_clean" %in% names(dt)) {
    rc <- as.character(dt$race_clean)
    rc[is.na(rc) | rc == ""] <- "Unknown"
    dt[, race_clean := factor(rc)]
  }
  if ("ethnicity_clean" %in% names(dt)) {
    ec <- as.character(dt$ethnicity_clean)
    ec[is.na(ec) | ec == ""] <- "Unknown"
    dt[, ethnicity_clean := factor(ec)]
  }
  invisible(NULL)
}

ensure_core_fields <- function(dt) {
  # Treatment
  if (!"treatment_binary" %in% names(dt) && "A" %in% names(dt)) {
    dt[, treatment_binary := as.integer(A)]
  }
  dt[, treatment_binary := as.integer(treatment_binary)]
  dt <- dt[treatment_binary %in% c(0L, 1L)]
  
  # Outcome
  if (!"event_flag" %in% names(dt) && "Y" %in% names(dt)) {
    dt[, event_flag := as.integer(Y)]
  }
  dt[, event_flag := as.integer(event_flag)]
  dt[is.na(event_flag), event_flag := 0L]
  
  # Standard names A and Y (TMLE scripts will use these)
  dt[, A := as.integer(treatment_binary)]
  dt[, Y := as.integer(event_flag)]
  
  dt
}

restrict_adults <- function(dt) {
  if (!"age_at_index" %in% names(dt)) return(dt)
  dt[, age_at_index := suppressWarnings(as.numeric(age_at_index))]
  # median impute then restrict adults
  dt[, age_at_index := median_impute(age_at_index)]
  dt <- dt[age_at_index >= 18]
  dt
}

add_person_time <- function(dt) {
  if ("time_to_event_days" %in% names(dt)) {
    dt[, pt_years := as.numeric(time_to_event_days) / 365.25]
  } else {
    dt[, pt_years := NA_real_]
  }
  dt
}

# ============================
# APS phenotype / severity covariates
# ============================

add_severity_covariates <- function(dt) {
  # Helper: safe 0/1 column fetch (missing -> 0)
  v01 <- function(v) {
    if (!v %in% names(dt)) return(rep.int(0L, nrow(dt)))
    x <- dt[[v]]
    x <- suppressWarnings(as.integer(as.numeric(x)))
    x[is.na(x)] <- 0L
    as.integer(x != 0L)
  }
  
  # --- Core phenotype flags (APS-relevant) ---
  # Arterial: stroke/TIA, MI/CAD, PAD
  dt[, arterial_thrombosis := as.integer(
    v01("stroke_tia") == 1L |
      v01("mi_cad") == 1L |
      v01("pad") == 1L
  )]
  
  # Venous: prior VTE
  dt[, venous_thrombosis := as.integer(v01("prior_vte") == 1L)]
  
  # Any thrombosis: arterial OR venous
  dt[, any_thrombosis := as.integer(arterial_thrombosis == 1L | venous_thrombosis == 1L)]
  
  # OPTIONAL: a simple “high-risk phenotype” proxy
  # (can be handy as a single binary interaction term)
  dt[, high_risk_phenotype := as.integer(
    arterial_thrombosis == 1L | (venous_thrombosis == 1L & v01("sle") == 1L)
  )]
  
  # Quick sanity coercion
  dt[, `:=`(
    arterial_thrombosis = as.integer(arterial_thrombosis),
    venous_thrombosis   = as.integer(venous_thrombosis),
    any_thrombosis      = as.integer(any_thrombosis),
    high_risk_phenotype = as.integer(high_risk_phenotype)
  )]
  
  invisible(NULL)
}



preprocess_one <- function(path_in, path_out) {
  cat("\n=== Reading:", path_in, "===\n")
  dt <- fread(path_in, na.strings = NA_STRINGS, showProgress = FALSE)
  setDT(dt)
  
  # Dates
  for (d in intersect(DATE_VARS, names(dt))) parse_date_col(dt, d)
  
  # Race/ethnicity
  derive_race_ethnicity(dt)
  
  # Gender
  standardize_gender(dt)
  
  # Age + restrict adults
  dt <- restrict_adults(dt)
  
  # Comorbidities NA->0; antibodies NA->0
  coerce_binary_01(dt, COMORB_VARS, na_to_zero = TRUE)
  coerce_binary_01(dt, AB_FLAGS,    na_to_zero = TRUE)
  add_severity_covariates(dt)
  
  # Labs (leave numeric; do not impute here unless you explicitly want to)
  for (lv in intersect(c("inr_value","platelet_value","creatinine_value"), names(dt))) {
    dt[, (lv) := suppressWarnings(as.numeric(as.character(get(lv))))]
  }
  
  # Make sure treatment/outcome exist and are 0/1
  dt <- ensure_core_fields(dt)
  
  # Drop constant binary flags (no variance)
  drop_constant_covariates(
    dt,
    c(COMORB_VARS, AB_FLAGS, "arterial_thrombosis", "venous_thrombosis", "any_thrombosis", "high_risk_phenotype")
  )
  
  # Fix time_to_event_days negatives
  fix_time_to_event(dt)
  
  # Add person-time
  dt <- add_person_time(dt)
  
  # Optional: ensure comorb_count exists
  if ("comorb_count" %in% names(dt)) {
    dt[, comorb_count := suppressWarnings(as.integer(as.numeric(comorb_count)))]
    dt[is.na(comorb_count), comorb_count := 0L]
  }
  
  # Last QC
  stopifnot(all(dt$A %in% c(0L,1L)))
  stopifnot(all(dt$Y %in% c(0L,1L)))
  
  cat("Final dims:", paste(dim(dt), collapse = " x "), "\n")
  cat("Deaths (Y=1):", sum(dt$Y == 1L, na.rm = TRUE), "\n")
  cat("Warfarin (A=0):", sum(dt$A == 0L), " DOAC (A=1):", sum(dt$A == 1L), "\n")
  cat("time_to_event_days < 0:", sum(dt$time_to_event_days < 0, na.rm = TRUE), "\n")
  cat("Arterial thrombosis = 1:", sum(dt$arterial_thrombosis == 1L, na.rm = TRUE), "\n")
  cat("Venous thrombosis = 1:", sum(dt$venous_thrombosis == 1L, na.rm = TRUE), "\n")
  
  # ============================
  # APS phenotype QC
  # ============================
  
  if (all(c("arterial_thrombosis","venous_thrombosis") %in% names(dt))) {
    
    cat("\n--- APS Phenotype Distribution ---\n")
    
    # 2x2 overlap table
    tab <- dt[, .N, by = .(arterial_thrombosis, venous_thrombosis)]
    print(tab)
    
    # Percent distribution
    tab[, pct := round(100 * N / sum(N), 2)]
    cat("\nPercent distribution:\n")
    print(tab)
    
    # Any thrombosis %
    if ("any_thrombosis" %in% names(dt)) {
      cat("\nAny thrombosis = 1:",
          sum(dt$any_thrombosis == 1L),
          " (",
          round(100 * mean(dt$any_thrombosis == 1L), 2),
          "%)\n", sep="")
    }
    
    # High-risk phenotype %
    if ("high_risk_phenotype" %in% names(dt)) {
      cat("High-risk phenotype = 1:",
          sum(dt$high_risk_phenotype == 1L),
          " (",
          round(100 * mean(dt$high_risk_phenotype == 1L), 2),
          "%)\n", sep="")
    }
    
    # Overlap proportion (both arterial & venous)
    both <- dt[, sum(arterial_thrombosis == 1L & venous_thrombosis == 1L)]
    cat("Both arterial & venous:",
        both,
        " (",
        round(100 * both / nrow(dt), 2),
        "%)\n\n", sep="")
  }
  
  
  fwrite(dt, path_out)
  cat("Wrote:", path_out, "\n")
  invisible(dt)
}

# ============================
# Run
# ============================
for (nm in names(INPUTS)) {
  if (!file.exists(INPUTS[[nm]])) stop("Missing input file: ", INPUTS[[nm]])
  preprocess_one(INPUTS[[nm]], OUTPUTS[[nm]])
}

cat("\nAll preprocessing complete.\nOutputs in: ", OUT_DIR, "\n")
