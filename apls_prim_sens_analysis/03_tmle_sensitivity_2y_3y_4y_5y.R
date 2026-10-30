#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tmle)
  library(SuperLearner)
  library(EValue)
  library(dplyr)
  library(mice)
})

set.seed(42)


# ============================
# Config
# ============================
IN_FILES <- list(
  # 1-year
  sens_1y  = "apls_prim_sens_analysis/input_data/tmle_ready/tmle_sens_1y.csv",
  # 2-year
  sens_2y       = "apls_prim_sens_analysis/input_data/tmle_ready/tmle_sens_2y.csv",
  sens_2y_labs  = "apls_prim_sens_analysis/input_data/tmle_ready/tmle_sens_2y_labs.csv",
  # 3-year
  sens_3y       = "apls_prim_sens_analysis/input_data/tmle_ready/tmle_sens_3y.csv",
  sens_3y_labs  = "apls_prim_sens_analysis/input_data/tmle_ready/tmle_sens_3y_labs.csv",
  # 4-year
  sens_4y       = "apls_prim_sens_analysis/input_data/tmle_ready/tmle_sens_4y.csv",
  sens_4y_labs  = "apls_prim_sens_analysis/input_data/tmle_ready/tmle_sens_4y_labs.csv",
  # 5-year
  sens_5y       = "apls_prim_sens_analysis/input_data/tmle_ready/tmle_sens_5y.csv",
  sens_5y_labs  = "apls_prim_sens_analysis/input_data/tmle_ready/tmle_sens_5y_labs.csv"
)

OUT_DIR <- "apls_prim_sens_analysis/output_data/final_dissertation_data/sensitivity_analysis_v1"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TRUNC_LEVELS <- c(0.01, 0.025)

SL_LIB <- c("SL.ranger", "SL.glmnet", "SL.mean")
#if (requireNamespace("randomForest", quietly = TRUE)) SL_LIB <- c(SL_LIB, "SL.randomForest")
#if (requireNamespace("xgboost", quietly = TRUE))      SL_LIB <- c(SL_LIB, "SL.xgboost")
#if (requireNamespace("gam", quietly = TRUE))          SL_LIB <- c(SL_LIB, "SL.gam")

# ============================
# Helpers
# ============================



is_binary01 <- function(x) {
  ux <- unique(x[!is.na(x)])
  (is.numeric(x) || is.integer(x)) && length(ux) <= 2 && all(ux %in% c(0, 1))
}


prep_for_SL <- function(dt, W_names) {
  stopifnot(all(c("A", "Y") %in% names(dt)))
  dt[, A := as.integer(A)]
  dt[, Y := as.integer(Y)]
  dt <- dt[!is.na(A) & !is.na(Y)]
  dt <- dt[A %in% c(0L, 1L) & Y %in% c(0L, 1L)]
  
  for (v in W_names) {
    x <- dt[[v]]
    
    if (is.character(x)) {
      x <- as.character(x)
      x[is.na(x) | trimws(x) == ""] <- "Unknown"
      dt[[v]] <- factor(x)
      next
    }
    
    if (is.factor(x)) {
      x <- as.character(x)
      x[is.na(x) | trimws(x) == ""] <- "Unknown"
      dt[[v]] <- factor(x)
      next
    }
    
    if (is.logical(x)) x <- as.integer(x)
    
    if (is.numeric(x) || is.integer(x)) {
      if (is_binary01(x)) {
        x[is.na(x)] <- 0L
        dt[[v]] <- as.integer(x)
      } else {
        miss_name <- paste0(v, "_missing")
        dt[[miss_name]] <- as.integer(is.na(x))
        med <- suppressWarnings(median(x, na.rm = TRUE))
        if (!is.finite(med)) med <- 0
        x[is.na(x)] <- med
        dt[[v]] <- as.numeric(x)
      }
      next
    }
    
    x <- as.character(x)
    x[is.na(x) | trimws(x) == ""] <- "Unknown"
    dt[[v]] <- factor(x)
  }
  
  dt
}



build_W <- function(dt) {
  drop_cols <- c(
    "personid",
    "apls_first_date", "index_date", "first_anticoag_date",
    "death_date", "censor_date",
    "event_flag", "treatment_binary", "treatment_group",
    "time_to_event_days", "pt_years",
    "races", "ethnicities", "comorb_count", 
    "death_before_index_flag", "any_thrombosis", "high_risk_phenotype",
    "A","Y",
    "inr_date" , "platelet_date", "creatinine_date",
    "inr_missing", "creatinine_missing",
    # dropping platelet_value due to 97% missingness
    "platelet_value",
    #dropping sensitivity, dis sev 
    "aps_severity_score", "aps_disease_severity"
  )
  
  
  W_names <- setdiff(names(dt), c(drop_cols, "A", "Y"))
  W_names <- setdiff(W_names, W_names[grepl("_date$", W_names)])
  
  if (length(W_names) == 0) stop("No covariates left for W after exclusions.")
  
  cat("Number of W covariates:", length(W_names), "\n")
  cat("W_names:", paste(W_names, collapse = ", "), "\n")
  
  return(list(W_names = W_names))
}



extract_tmle_summary <- function(fit, analysis_label, gbound) {
  get_est <- function(path, default = NA_real_) {
    tryCatch({
      x <- fit$estimates
      for (p in path) x <- x[[p]]
      x
    }, error = function(e) default)
  }
  
  rd_ci <- get_est(c("ATE", "CI"), default = c(NA_real_, NA_real_))
  rr_ci <- get_est(c("RR",  "CI"), default = c(NA_real_, NA_real_))
  or_ci <- get_est(c("OR",  "CI"), default = c(NA_real_, NA_real_))
  rr_val <- get_est(c("RR", "psi"))
  ev_point <- NA_real_
  ev_lower <- NA_real_
  
  if (!is.na(rr_val) && !is.na(rr_ci[1]) && !is.na(rr_ci[2])) {
    # We wrap this in try() just in case the EValue package 
    # encounters a math error (like log of zero)
    try({
      res <- evalues.RR(est = rr_val, lo = rr_ci[1], hi = rr_ci[2])
      
      # Using as.numeric() on the second row ignores the names/labels 
      # that are causing the "out of bounds" error.
      # The EValue matrix structure: 
      # Row 1: Original RR, Lower, Upper
      # Row 2: E-Value Point, E-Value CI bound
      
      ev_point <- res[2, 1] # Row 2, Column 1 (Point Estimate E-value)
      ev_lower <- res[2, 2] # Row 2, Column 2 (CI E-value)
    }, silent = TRUE)
  }
  
  data.table(
    analysis = analysis_label,
    gbound = gbound,
    n = fit$n,
    risk_doac     = get_est(c("EY1", "psi")),
    risk_warfarin = get_est(c("EY0", "psi")),
    RD            = get_est(c("ATE", "psi")),
    RD_ci_l       = rd_ci[1],
    RD_ci_u       = rd_ci[2],
    RR            = rr_val,
    RR_ci_l       = rr_ci[1],
    RR_ci_u       = rr_ci[2],
    OR            = get_est(c("OR", "psi")),
    OR_ci_l       = or_ci[1],
    OR_ci_u       = or_ci[2],
    e_value_est   = ev_point,
    e_value_lower = ev_lower
  )
}

fit_g_q_predictions <- function(dt, W_df, SL_LIB, label = NULL) {
  g_fit <- SuperLearner(
    Y = dt$A,
    X = W_df,
    family = binomial(),
    SL.library = SL_LIB,
    method = "method.NNLS", 
    verbose = TRUE
  )
  cat("Number of W_df covariates:", length(W_df), "\n")
  cat("W_df:", paste(colnames(W_df), collapse = ", "), "\n")
  cat("\n--- SuperLearner Weights: g-model (", label, ") ---\n", sep = "")
  print(round(g_fit$coef, 4))
  g1W <- as.numeric(g_fit$SL.predict)
  
  X_Q <- data.frame(A = dt$A, W_df)
  Q_fit <- SuperLearner(
    Y = dt$Y,
    X = X_Q,
    family = binomial(),
    SL.library = SL_LIB,
    method = "method.NNLS",
    verbose = TRUE
  )
  cat("\n--- SuperLearner Weights: Q-model (", label, ") ---\n", sep = "")
  print(round(Q_fit$coef, 4))
  X_Q0 <- X_Q; X_Q0$A <- 0L
  X_Q1 <- X_Q; X_Q1$A <- 1L
  Q0W <- as.numeric(predict(Q_fit, newdata = X_Q0)$pred)
  Q1W <- as.numeric(predict(Q_fit, newdata = X_Q1)$pred)
  
  eps <- 1e-6
  g1W <- pmin(pmax(g1W, eps), 1 - eps)
  Q0W <- pmin(pmax(Q0W, eps), 1 - eps)
  Q1W <- pmin(pmax(Q1W, eps), 1 - eps)
  
  list(g1W = g1W, Q0W = Q0W, Q1W = Q1W)
}

run_tmle_with_explicit_gq <- function(dt, W_names, label, TRUNC_LEVELS, SL_LIB) {
  W_df <- as.data.frame(dt[, ..W_names])
  
  preds <- fit_g_q_predictions(dt = dt, W_df = W_df, SL_LIB = SL_LIB)
  
  # optional diagnostics (kept, but not printed unless you want)
  # print(summary(preds$g1W))
  # print(quantile(preds$g1W, c(0.001,0.01,0.05,0.5,0.95,0.99,0.999)))
  
  rbindlist(lapply(TRUNC_LEVELS, function(gb) {
    fit <- tmle(
      Y = dt$Y,
      A = dt$A,
      W = W_df,
      Q = cbind(preds$Q0W, preds$Q1W),
      g1W = preds$g1W,
      family = "binomial",
      gbound = gb,
      verbose = TRUE
    )
    extract_tmle_summary(fit, label, gb)
  }))
}

run_one <- function(label, path) {
  if (!file.exists(path)) stop("Missing TMLE-ready input file: ", path)
  
  dt <- fread(path)
  setDT(dt)
  # remove patients that have death date before index date
  dt <- dt[dt$death_before_index_flag == 0, ]
  
  ############### IMPUTATION STARTS ################
  
  # identify which lab variables are actually present in this file
  lab_targets <- intersect(c("inr_value", "creatinine_value"), names(dt))
  
  if (length(lab_targets) > 0) {
    
    # create missingness indicators only for labs that exist
    if ("inr_value" %in% lab_targets) {
      dt[, inr_missing := as.integer(is.na(inr_value))]
    }
    if ("creatinine_value" %in% lab_targets) {
      dt[, creatinine_missing := as.integer(is.na(creatinine_value))]
    }
    
    # candidate predictors for imputation
    predictor_candidates <- c(
      "age_at_index", "gender", "race_clean", "ethnicity_clean",
      "htn", "dm", "dyslipidemia", "obesity", "tobacco",
      "stroke_tia", "mi_cad", "pad", "prior_vte", "ckd",
      "liver_dz", "cancer", "copd", "heart_failure",
      "atrial_fib", "valvular", "sle", "prior_major_bleed",
      "pud", "anemia", "arterial_thrombosis", "venous_thrombosis",
      "la_pos", "b2gpi_pos",
      "A", "Y"
    )
    
    # only keep predictors that are actually in this dataset
    predictor_vars <- intersect(predictor_candidates, names(dt))
    
    # final imputation dataset = labs present + predictors present
    imp_vars <- unique(c(lab_targets, predictor_vars))
    imp_data <- as.data.frame(dt[, ..imp_vars])
    
    # make character variables factors for mice
    for (v in names(imp_data)) {
      if (is.character(imp_data[[v]])) imp_data[[v]] <- factor(imp_data[[v]])
    }
    
    # initialize mice
    meth <- make.method(imp_data)
    pred <- make.predictorMatrix(imp_data)
    
    # default: do not impute anything
    meth[] <- ""
    
    # impute only the lab variables that exist
    for (v in lab_targets) {
      meth[v] <- "pmm"
    }
    
    # do not use a variable to predict itself
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
    
    # write imputed values back only for labs that were imputed
    for (v in lab_targets) {
      dt[, (v) := imp_complete[[v]]]
    }
    
    cat("Imputed lab variables for ", label, ": ", paste(lab_targets, collapse = ", "), "\n", sep = "")
    
  } else {
    cat("No INR/creatinine lab columns found for ", label, "; skipping imputation.\n", sep = "")
  }
  
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
  
  if (!("A" %in% names(dt)) && ("treatment_binary" %in% names(dt))) dt[, A := as.integer(treatment_binary)]
  if (!("Y" %in% names(dt)) && ("event_flag" %in% names(dt)))       dt[, Y := as.integer(event_flag)]
  stopifnot(all(c("A", "Y") %in% names(dt)))
  
  w_info <- build_W(dt)
  W_names <- w_info$W_names
  
  dt <- prep_for_SL(dt, W_names)
  const_w <- W_names[sapply(W_names, function(v) length(unique(dt[[v]])) <= 1)]
  if (length(const_w) > 0) {
    cat("Dropped constant W covariates:", paste(const_w, collapse = ", "), "\n")
    W_names <- setdiff(W_names, const_w)
  }
  if (length(W_names) == 0) stop("All W covariates became constant after prep.")
  
  
  cat("\n=== ", label, " ===\n", sep = "")
  cat("File:", path, "\n")
  cat("N =", nrow(dt), "Deaths(Y=1) =", sum(dt$Y == 1L), "\n")
  cat("A=1 =", sum(dt$A == 1L), "A=0 =", sum(dt$A == 0L), "\n")
  cat("W covariates =", length(W_names), "\n")
  cat("Total NA in W after prep =", sum(is.na(dt[, ..W_names])), "\n")
  
  run_tmle_with_explicit_gq(dt, W_names, label, TRUNC_LEVELS, SL_LIB)
}

# ============================
# Run all
# ============================
all_results <- rbindlist(
  Map(run_one, names(IN_FILES), IN_FILES),
  fill = TRUE
)

out_csv <- file.path(OUT_DIR, "tmle_sensitivity_summaries.csv")
fwrite(all_results, out_csv)

out_rds <- file.path(OUT_DIR, "tmle_sensitivity_run_meta.rds")
saveRDS(
  list(
    summary = all_results,
    SL_LIB = SL_LIB,
    trunc_levels = TRUNC_LEVELS,
    inputs = IN_FILES
  ),
  out_rds
)

cat("\nSaved:\n- ", out_csv, "\n- ", out_rds, "\n", sep = "")
