#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tmle)
  library(SuperLearner)
  library(ranger)
  library(dplyr)
  library(mice)
})

# ============================
# Config
# ============================
IN_FILE  <- "apls_prim_sens_analysis/input_data/tmle_ready/tmle_primary_1y_labs.csv"
OUT_DIR  <- "apls_prim_sens_analysis/output_data/final_dissertation_data/primary_v1"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TRUNC_LEVELS <- c(0.01, 0.025)

# Pick SuperLearner library
# (you can add/remove learners based on what you have installed)
SL_LIB <- c("SL.glmnet", "SL.mean", "SL.ranger")

# Optional learners if installed
#if (requireNamespace("randomForest", quietly = TRUE)) SL_LIB <- c(SL_LIB, "SL.randomForest")
#if (requireNamespace("xgboost", quietly = TRUE))      SL_LIB <- c(SL_LIB, "SL.xgboost")
#if (requireNamespace("gam", quietly = TRUE))          SL_LIB <- c(SL_LIB, "SL.gam")

# ============================
# Helper: choose covariates W
# ============================
build_W <- function(dt) {
  drop_cols <- c(
    "personid",
    "apls_first_date","index_date","first_anticoag_date",
    "death_date","censor_date",
    "event_flag","time_to_event_days","pt_years",
    "treatment_group","treatment_binary",
    "A","Y",
    "races","ethnicities", "comorb_count", 
    "inr_date" , "platelet_date", "creatinine_date",
    "inr_missing", "creatinine_missing",
    #"arterial_thrombosis", "venous_thrombosis", 
    "any_thrombosis","high_risk_phenotype",
    #dropping platelet_value due to 97% missingness
    "platelet_value",
    #dropping sensitivity, dis sev 
    "aps_severity_score_sens", "aps_disease_severity_sens"
  )
  
  
  # Keep race/ethnicity/gender if present (as factors)
  # Make sure any date columns remaining are dropped
  W_names <- setdiff(names(dt), drop_cols)
  
  # Remove any residual date-like columns
  date_like <- W_names[sapply(dt[, ..W_names], inherits, what = "IDate") |
                         grepl("_date$", W_names)]
  W_names <- setdiff(W_names, date_like)
  
  # Ensure W is not empty
  if (length(W_names) == 0) stop("No covariates left for W after exclusions.")
  
  # Convert characters to factors (tmle will one-hot via model.matrix inside learners)
  for (v in W_names) {
    if (is.character(dt[[v]])) dt[[v]] <- as.factor(dt[[v]])
  }
  list(dt = dt, W_names = W_names)
}

extract_tmle_summary <- function(fit, analysis_label, gbound) {
  # tmle() often returns estimates for ATE, RR, OR
  get_est <- function(path, default = NA_real_) {
    tryCatch({
      x <- fit$estimates
      for (p in path) x <- x[[p]]
      x
    }, error = function(e) default)
  }
  
  out <- data.table(
    analysis = analysis_label,
    gbound = gbound,
    risk_doac    = get_est(c("EY1","psi")),
    risk_warfarin= get_est(c("EY0","psi")),
    RD           = get_est(c("ATE","psi")),
    RD_ci_l      = get_est(c("ATE","CI"))[1],
    RD_ci_u      = get_est(c("ATE","CI"))[2],
    RR           = get_est(c("RR","psi")),
    RR_ci_l      = get_est(c("RR","CI"))[1],
    RR_ci_u      = get_est(c("RR","CI"))[2],
    OR           = get_est(c("OR","psi")),
    OR_ci_l      = get_est(c("OR","CI"))[1],
    OR_ci_u      = get_est(c("OR","CI"))[2]
  )
  out
}

# ============================
# Run
# ============================
dt <- fread(IN_FILE)
setDT(dt)

# remove patients that have death date before index date
dt <- dt[dt$death_before_index_flag == 0, ]

# handle NA values in inr_value and creatinine_value
#dt$inr_value[is.na(dt$inr_value)] <- median(dt$inr_value, na.rm = TRUE)
#dt$creatinine_value[is.na(dt$creatinine_value)] <- median(dt$creatinine_value, na.rm = TRUE)

############### IMPUTATION STARTS ################

# handle missing baseline labs using mice + predictive mean matching (PMM)
# create missingness indicators so the analysis model can retain missingness information
dt[, inr_missing := as.integer(is.na(inr_value))]
dt[, creatinine_missing := as.integer(is.na(creatinine_value))]

# variables used to impute the lab values
imp_vars <- c(
  "inr_value", "creatinine_value",
  "age_at_index", "gender", "race_clean", "ethnicity_clean",
  "htn", "dm", "dyslipidemia", "obesity", "tobacco",
  "stroke_tia", "mi_cad", "pad", "prior_vte", "ckd",
  "liver_dz", "cancer", "copd", "heart_failure",
  "atrial_fib", "valvular", "sle", "prior_major_bleed",
  "pud", "anemia", "arterial_thrombosis", "venous_thrombosis",
  "la_pos", "b2gpi_pos",
  "A", "Y"
)

imp_data <- as.data.frame(dt[, ..imp_vars])

# make character variables explicit factors for mice
for (v in names(imp_data)) {
  if (is.character(imp_data[[v]])) imp_data[[v]] <- factor(imp_data[[v]])
}

# initialize mice setup
meth <- make.method(imp_data)
pred <- make.predictorMatrix(imp_data)

# only impute the two lab variables
meth[] <- ""
meth[c("inr_value", "creatinine_value")] <- "pmm"

# do not use each variable to predict itself
diag(pred) <- 0

set.seed(42)
imp <- mice(
  imp_data,
  m = 1,
  maxit = 10,
  method = meth,
  predictorMatrix = pred,
  printFlag = FALSE
)

imp_complete <- complete(imp, 1)

dt[, inr_value := imp_complete$inr_value]
dt[, creatinine_value := imp_complete$creatinine_value]

############### IMPUTATION ENDS ################

############### DISEASE SEVERITY SCORE ADDITION BEGINS ################
add_aps_severity <- function(dt, include_la_sensitivity = TRUE) {
  # Works with data.table or data.frame
  # Adds:
  #   renal_dysfunction
  #   aps_severity_score
  #   aps_disease_severity
  # Optionally adds:
  #   aps_severity_score_sens
  #   aps_disease_severity_sens
  
  # ----------------------------
  # Helpers
  # ----------------------------
  v01 <- function(v) {
    if (!v %in% names(dt)) return(rep.int(0L, nrow(dt)))
    x <- dt[[v]]
    x <- suppressWarnings(as.integer(as.numeric(x)))
    x[is.na(x)] <- 0L
    as.integer(x != 0L)
  }
  
  vnum <- function(v) {
    if (!v %in% names(dt)) return(rep(NA_real_, nrow(dt)))
    suppressWarnings(as.numeric(dt[[v]]))
  }
  
  # ----------------------------
  # Core APS-relevant components
  # ----------------------------
  arterial <- v01("arterial_thrombosis")
  venous   <- v01("venous_thrombosis")
  sle      <- v01("sle")
  ckd      <- v01("ckd")
  b2gpi    <- v01("b2gpi_pos")
  la       <- v01("la_pos")
  bleed    <- v01("prior_major_bleed")
  creat    <- vnum("creatinine_value")
  
  # Renal dysfunction:
  # CKD flag OR elevated creatinine
  renal_dysfunction <- as.integer(
    ckd == 1L | (!is.na(creat) & creat >= 1.5)
  )
  
  # ----------------------------
  # Primary APS severity score
  # ----------------------------
  # Recommended formula:
  #   2 * arterial_thrombosis
  # + 1 * venous_thrombosis
  # + 1 * sle
  # + 1 * renal_dysfunction
  # + 1 * b2gpi_pos
  # + 1 * prior_major_bleed
  aps_severity_score <-
    2L * arterial +
    1L * venous +
    1L * sle +
    1L * renal_dysfunction +
    1L * b2gpi +
    1L * bleed
  
  # Categorical version for tables / subgroup analyses
  aps_disease_severity <- ifelse(
    aps_severity_score <= 1L, "low",
    ifelse(aps_severity_score <= 3L, "moderate", "high")
  )
  
  aps_disease_severity <- factor(
    aps_disease_severity,
    levels = c("low", "moderate", "high")
  )
  
  # ----------------------------
  # Optional sensitivity version
  # ----------------------------
  # Because la_pos is very sparse, do not use it in the primary score.
  # If desired, add it only in a sensitivity version.
  if (isTRUE(include_la_sensitivity)) {
    aps_severity_score_sens <- aps_severity_score + 1L * la
    
    aps_disease_severity_sens <- ifelse(
      aps_severity_score_sens <= 1L, "low",
      ifelse(aps_severity_score_sens <= 3L, "moderate", "high")
    )
    
    aps_disease_severity_sens <- factor(
      aps_disease_severity_sens,
      levels = c("low", "moderate", "high")
    )
  }
  
  # ----------------------------
  # Assign back
  # ----------------------------
  if (inherits(dt, "data.table")) {
    dt[, renal_dysfunction := as.integer(renal_dysfunction)]
    dt[, aps_severity_score := as.integer(aps_severity_score)]
    dt[, aps_disease_severity := aps_disease_severity]
    
    if (isTRUE(include_la_sensitivity)) {
      dt[, aps_severity_score_sens := as.integer(aps_severity_score_sens)]
      dt[, aps_disease_severity_sens := aps_disease_severity_sens]
    }
    
    return(invisible(dt))
  } else {
    dt$renal_dysfunction    <- as.integer(renal_dysfunction)
    dt$aps_severity_score   <- as.integer(aps_severity_score)
    dt$aps_disease_severity <- aps_disease_severity
    
    if (isTRUE(include_la_sensitivity)) {
      dt$aps_severity_score_sens   <- as.integer(aps_severity_score_sens)
      dt$aps_disease_severity_sens <- aps_disease_severity_sens
    }
    
    return(dt)
  }
}

# data.table
add_aps_severity(dt)

# check distribution
table(dt$aps_disease_severity, useNA = "ifany")
summary(dt$aps_severity_score)

# optional sensitivity version
add_aps_severity(dt, include_la_sensitivity = TRUE)

############### DISEASE SEVERITY SCORE ADDITION ENDS ################


stopifnot(all(c("A","Y") %in% names(dt)))
dt[, A := as.integer(A)]
dt[, Y := as.integer(Y)]
dt <- dt[A %in% c(0L,1L)]

############### SAve the final datafile Begins ################
dt_final <- copy(dt)
FINAL_DATA_CSV <- file.path("apls_prim_sens_analysis/input_data/tmle_ready/tmle_primary_1y_final_dataset.csv")
fwrite(dt_final, FINAL_DATA_CSV)

cat("\nSaved final processed dataset to:", FINAL_DATA_CSV, "\n")
############### SAve the final datafile ENDS ################


# Build W
tmp <- build_W(dt)
dt <- tmp$dt
W_names <- tmp$W_names
W <- dt[, ..W_names]
W_df <- as.data.frame(W)

cat("N =", nrow(dt), "\n")
cat("A=1 (DOAC) =", sum(dt$A == 1L), " A=0 (Warfarin) =", sum(dt$A == 0L), "\n")
cat("Deaths (Y=1) =", sum(dt$Y == 1L), "\n")
cat("W variables =", length(W_names), "\n")
cat("SuperLearner library:\n"); print(SL_LIB)

# 1) Propensity model g(W) = P(A=1 | W)
g_fit <- SuperLearner(
  Y = dt$A,
  X = W_df,
  family = binomial(),
  SL.library = SL_LIB,
  method = "method.NNLS",
  verbose=TRUE
)
g_fit$coef
g1W <- as.numeric(g_fit$SL.predict)  # P(A=1|W)

# 2) Outcome model Q(A,W) = P(Y=1 | A, W) using one model with A included
X_Q <- data.frame(A = dt$A, W_df)

Q_fit <- SuperLearner(
  Y = dt$Y,
  X = X_Q,
  family = binomial(),
  SL.library = SL_LIB,
  method = "method.NNLS",
  verbose=TRUE
)

# Predict counterfactual risks under A=0 and A=1
X_Q0 <- X_Q; X_Q0$A <- 0L
X_Q1 <- X_Q; X_Q1$A <- 1L

Q0W <- as.numeric(predict(Q_fit, newdata = X_Q0)$pred)
Q1W <- as.numeric(predict(Q_fit, newdata = X_Q1)$pred)

# Defensive bounds (optional but recommended for stability)
eps <- 1e-6
g1W <- pmin(pmax(g1W, eps), 1 - eps)
Q0W <- pmin(pmax(Q0W, eps), 1 - eps)
Q1W <- pmin(pmax(Q1W, eps), 1 - eps)

results <- rbindlist(lapply(TRUNC_LEVELS, function(gb) {
  
  fit <- tmle(
    Y = dt$Y,
    A = dt$A,
    W = W_df,                 # still pass W for clever covariate construction
    Q = cbind(Q0W, Q1W),       # <- your explicit Q predictions
    g1W = g1W,                 # <- your explicit g predictions (P(A=1|W))
    family = "binomial",
    gbound = gb,
    verbose = TRUE
  )
  
  extract_tmle_summary(fit, analysis_label = "primary_1y_labs", gbound = gb)
}))

out_csv <- file.path(OUT_DIR, "tmle_primary_1y_summary.csv")
fwrite(results, out_csv)

out_rds <- file.path(OUT_DIR, "tmle_primary_1y_fit_objects.rds")
saveRDS(list(summary = results,
             W_names = W_names,
             SL_LIB = SL_LIB),
        out_rds)

cat("\nSaved:\n- ", out_csv, "\n- ", out_rds, "\n", sep = "")

# ============================
# CRUDE Analysis
# ============================

OUT_CSV <- file.path(OUT_DIR, "crude_primary_1y.csv")

# If your treatment is already in A (0=warfarin, 1=DOAC), use that.
# Otherwise set TREAT_COL to "treatment_binary" or similar.
TREAT_COL <- "A"
OUTCOME_COL <- "Y"

# ============================
# Load + basic cleaning
# ============================
dt <- fread(IN_FILE)
setDT(dt)

stopifnot(all(c(TREAT_COL, OUTCOME_COL) %in% names(dt)))

# Coerce to integer 0/1
dt[, (TREAT_COL) := as.integer(get(TREAT_COL))]
dt[, (OUTCOME_COL) := as.integer(get(OUTCOME_COL))]

# Keep only valid 0/1 rows for A and Y
dt <- dt[get(TREAT_COL) %in% c(0L, 1L) & get(OUTCOME_COL) %in% c(0L, 1L)]

A <- dt[[TREAT_COL]]
Y <- dt[[OUTCOME_COL]]

# ============================
# Crude estimates
# ============================
n_total <- nrow(dt)
n_a1 <- sum(A == 1L)
n_a0 <- sum(A == 0L)

events_a1 <- sum(Y[A == 1L] == 1L)
events_a0 <- sum(Y[A == 0L] == 1L)

risk_a1 <- mean(Y[A == 1L] == 1L)  # DOAC crude risk
risk_a0 <- mean(Y[A == 0L] == 1L)  # Warfarin crude risk

rd <- risk_a1 - risk_a0
rr <- risk_a1 / risk_a0

# Odds ratio (with guard for zeros)
odds_a1 <- risk_a1 / (1 - risk_a1)
odds_a0 <- risk_a0 / (1 - risk_a0)
or <- odds_a1 / odds_a0

# ============================
# Simple Wald 95% CI for RD
# (risk difference CI; OK as a quick crude check)
# ============================
se_rd <- sqrt((risk_a1 * (1 - risk_a1)) / n_a1 + (risk_a0 * (1 - risk_a0)) / n_a0)
rd_ci <- rd + c(-1, 1) * 1.96 * se_rd

# ============================
# Log-scale Wald 95% CI for RR and OR (quick crude)
# ============================
# RR: use standard log(RR) SE using counts; add small continuity correction if needed
cc <- 0.5
a <- events_a1
b <- n_a1 - events_a1
c_ <- events_a0
d <- n_a0 - events_a0

if (any(c(a, b, c_, d) == 0)) {
  a <- a + cc; b <- b + cc; c_ <- c_ + cc; d <- d + cc
}

rr_hat <- (a / (a + b)) / (c_ / (c_ + d))
se_log_rr <- sqrt((1 / a) - (1 / (a + b)) + (1 / c_) - (1 / (c_ + d)))
rr_ci <- exp(log(rr_hat) + c(-1, 1) * 1.96 * se_log_rr)

or_hat <- (a * d) / (b * c_)
se_log_or <- sqrt(1 / a + 1 / b + 1 / c_ + 1 / d)
or_ci <- exp(log(or_hat) + c(-1, 1) * 1.96 * se_log_or)

# ============================
# Output
# ============================
out <- data.table(
  analysis = "crude_primary_1y",
  n_total = n_total,
  n_doac = n_a1,
  n_warfarin = n_a0,
  deaths_doac = events_a1,
  deaths_warfarin = events_a0,
  risk_doac = risk_a1,
  risk_warfarin = risk_a0,
  RD = rd,
  RD_ci_l = rd_ci[1],
  RD_ci_u = rd_ci[2],
  RR = rr_hat,
  RR_ci_l = rr_ci[1],
  RR_ci_u = rr_ci[2],
  OR = or_hat,
  OR_ci_l = or_ci[1],
  OR_ci_u = or_ci[2]
)

print(out)

dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
fwrite(out, OUT_CSV)
cat("\nSaved crude results to:", OUT_CSV, "\n")

