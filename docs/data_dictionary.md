# Data dictionary: expected input files

Input files are placed in `apls_prim_sens_analysis/input_data/raw_data/`. Each file is one row per patient for one follow-up horizon, with or without baseline laboratory values:

`df_1y_df.csv`, `df_1y_labs.csv`, `df_2y_df.csv`, `df_2y_labs.csv`, `df_3y_df.csv`, `df_3y_labs.csv`, `df_4y_df.csv`, `df_4y_labs.csv`, `df_5y_df.csv`, `df_5y_labs_df.csv`

| Column | Description |
| --- | --- |
| `personid` | De-identified patient identifier (dropped before modeling) |
| `apls_first_date` | Date of first recorded APS diagnosis (ICD-10-CM D68.61) |
| `index_date`, `first_anticoag_date` | Date of first qualifying warfarin or DOAC record after diagnosis (time zero) |
| `death_date` | Recorded death date, if any |
| `censor_date` | End of horizon-specific follow-up |
| `event_flag` | 1 if death recorded within the horizon |
| `time_to_event_days` | Days from index to death or end of horizon |
| `age_at_index` | Age in years at index (median-imputed if missing, before adult restriction) |
| `gender`, `races`, `ethnicities` | Recorded sex, race and ethnicity |
| `treatment_group`, `treatment_binary` | Initial anticoagulant class (1 = DOAC, 0 = warfarin) |
| `htn` … `anemia` | Coded comorbidity indicators recorded on or before index (see manuscript S2 Table) |
| `comorb_count` | Number of comorbidities (descriptive only) |
| `death_before_index_flag` | 1 if recorded death on or before index (excluded) |
| `acl_pos`, `b2gpi_pos`, `la_pos`, `triple_positive` | Documented antibody positivity on or before index (untested coded 0) |
| `inr_value`, `inr_date` | Baseline international normalized ratio (laboratory files only) |
| `platelet_value`, `platelet_date` | Baseline platelet count (not used; about 97% missing) |
| `creatinine_value`, `creatinine_date` | Baseline serum creatinine (laboratory files only) |
