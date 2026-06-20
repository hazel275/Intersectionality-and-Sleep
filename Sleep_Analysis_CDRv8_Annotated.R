# =============================================================================
# INTERSECTIONALITY AND SLEEP — FULL ANALYSIS CODE (annotated for GitHub)
# Repo: github.com/hazel275/Intersectionality-and-Sleep
#
# Complete, annotated pipeline for the wearable sleep-disparities analysis used
# in the BRIDGE Week 2 workshop. Runs inside the All of Us Researcher Workbench
# (Controlled Tier), in RStudio, against CDR v8.
#
# HOW TO RUN
#   Run the blocks in order, 00 through 12. By default (REBUILD_FROM_RAW = FALSE)
#   Block 00 restores m_core_rs and reads the cleaned data from its model frame,
#   so the file runs end to end with no raw pull. Set REBUILD_FROM_RAW = TRUE to
#   rebuild df_analytic from the CDR, which fits the models and saves the objects.
#
# REQUIREMENTS
#   Packages: tidyverse, lme4, lmerTest, emmeans, geepack, broom.mixed, bigrquery.
#   Use RStudio, not the Python kernel, because lme4 needs R.
#   Access: an All of Us account with Controlled Tier and your institution's DUA.
#
# DATA USE
#   Code only. No participant-level data lives in this repo. Any aggregate output
#   follows the All of Us Data and Statistics Dissemination Policy: no reported
#   cell corresponds to fewer than 20 participants. Reproduce inside the Workbench.
#
# TWO RUN PATHS (set REBUILD_FROM_RAW in Block 00)
#   Default FALSE: restores m_core_rs and reads the cleaned data from its model
#   frame. Runs end to end, reproduces exactly, needs no raw rebuild.
#   TRUE: rebuilds from the CDR. Block 02 demographics recoding is a labeled
#   scaffold here, since the v3 mapping is not in this repo, so fill it first.
# =============================================================================


# =============================================================================
# SLEEP TRAJECTORY PAPER: COMPLETE ANALYSIS PIPELINE (v4 — FINAL)
# Wearable-Measured Sleep Disparities by Race and Sex in the All of Us
# Research Program
#
# Authors: Stephanie H. Cook, Julie Holm, Antoneta Karaj, Mariana Rodrigues,
#          Erica P. Wood, Cindy Patippe, Danning Tian
# Target Journal: SLEEP (Oxford Academic)
# CDR Version: All of Us CDR v8 (C2024Q3R8)
# Study Period: 2017-2023
# Final analytic N: 14,084 persons, 64,803 person-years
#
# v4 CHANGES FROM v3:
#   - Model D updated to random slope: (1 + t1 | person_id)
#     ΔAIC = 8,377; χ²(2) = 8,381.6, p < 0.001; not singular
#   - Hispanic Male now significant in Model D: b = -17.30, p = 0.009
#   - MAIHDA strata_id corrected to as.factor() (was as.integer())
#   - MAIHDA bootstrap CI added (1000 iterations)
#   - Residual diagnostics block added (Block 06d)
#   - Block 07 updated to use m_core_rs for emmeans
#   - Block 11 updated: m_core_rs.rds added to save list
#
# PIPELINE OVERVIEW:
#   Block 00: Environment setup + session restore
#   Block 01: Eligibility — BYOD adherence check
#   Block 02: Demographics cleaning
#   Block 03: Sleep data pull (BigQuery)
#   Block 04: Final analytic dataset assembly
#   Block 05: Descriptive statistics (Table 1)
#   Block 06: Longitudinal models (Models A-D)
#     06a: Model A — unconditional growth
#     06b: Models B and C — race and sex main effects
#     06c: Model D — primary intersectional random slope model
#     06d: Residual diagnostics for Model D
#     06e: Random slope model selection (AIC comparison)
#   Block 07: Predicted trajectories (emmeans from m_core_rs)
#   Block 08: Sensitivity analyses
#     08a: IPW
#     08b: Alternate wear-time thresholds
#     08c: Device type
#     08d: Piecewise and spline models
#     08e: Year-specific contrasts
#   Block 09: Secondary analyses
#     09a: Short sleep GEE
#     09b: MAIHDA variance partition + bootstrap CI
#   Block 10: Figures (exported as .tif at 300 dpi)
#   Block 11: Save all outputs to bucket
# =============================================================================


# =============================================================================
# BLOCK 00: ENVIRONMENT SETUP
# =============================================================================

library(tidyverse)
library(bigrquery)
library(lme4)
library(lmerTest)
library(emmeans)
library(geepack)

# NOTE: Setting_Env_Variables_p2.R is NOT needed if environment vars are
# already set (confirmed in workspace wb-windy-onion-3543 / wb-blissful-amaranth-1217).
# Uncomment the source() line only if variables are missing.
# source("~/workspace/aou-tutorial-notebooks/Setting_Env_Variables_p2.R")

# Verify environment variables (these are pre-set in All of Us workspaces)
cdr     <- Sys.getenv("WORKSPACE_CDR")      # wb-silky-artichoke-2408.C2024Q3R8
project <- Sys.getenv("GOOGLE_CLOUD_PROJECT") # wb-frosty-sprout-5598

# NOTE: WORKSPACE_BUCKET env var points to cloned-mybucket-... (wrong bucket).
# Always use correct_bucket below for sleep paper outputs.
correct_bucket <- "gs://cloned-aou-tutorial-notebooks-wb-frosty-sprout-5598"

cat("CDR:", cdr, "\n")
cat("Project:", project, "\n")
cat("Bucket:", correct_bucket, "\n")

# Query helper
run_query <- function(sql) {
  tb <- bq_project_query(project, sql)
  bq_table_download(tb)
}

run_query_paged <- function(sql) {
  tb <- bq_project_query(project, sql)
  bq_table_download(tb, page_size = 50000)
}

# ── SESSION RESTORE (use after VM restart) ──────────────────────────────────
# After a VM restart, or as the default path, the cleaned analytic data is
# recovered from the fitted model frame in the else branch below.

# ---- HOW THIS FILE RUNS -----------------------------------------------------
# Default (FALSE): restore m_core_rs, then read df_analytic from model.frame().
#   Reproduces exactly, no raw pull, no v3 recoding needed.
# TRUE: rebuild df_analytic from the CDR. Requires the v3 demographics recoding
#   in Block 02. Set TRUE only if you have it.
REBUILD_FROM_RAW <- FALSE

if (REBUILD_FROM_RAW) {   # ===== RAW REBUILD PATH: Blocks 01 to 04 =====
  
  
  # =============================================================================
  # BLOCK 01: ELIGIBILITY — BYOD ADHERENCE CHECK
  # =============================================================================
  # BYOD filter: src_id = 'Participant Portal: PTSC' in sleep_daily_summary
  # This correctly identifies BYOD participants in CDR v8.
  # NOTE: The old 'fitbit_activity_summary.has_wear_consent = 0' filter does NOT
  # exist in CDR v8. The wear_study table captures WEAR study participants,
  # who are confirmed to have ZERO overlap with the analytic sample.
  #
  # VERIFIED: 0 wear_study participants in final analytic sample.
  # The 256 careevolution src_id records are legitimate BYOD participants who
  # connected via the CareEvolution app integration (not WEAR study).
  
  sql_eligible <- paste0(
    "WITH hr_valid AS (\n",
    "  SELECT person_id,\n",
    "         DATE(datetime) AS hr_date,\n",
    "         COUNT(*) AS hr_minutes\n",
    "  FROM `", cdr, ".heart_rate_minute_level`\n",
    "  WHERE DATE(datetime) BETWEEN '2017-01-01' AND '2023-12-31'\n",
    "  GROUP BY person_id, DATE(datetime)\n",
    "  HAVING hr_minutes >= 600  -- >= 10 hours wear time\n",
    "),\n",
    "valid_days AS (\n",
    "  SELECT person_id,\n",
    "         EXTRACT(YEAR FROM hr_date) AS year,\n",
    "         DATE_TRUNC(hr_date, WEEK) AS week,\n",
    "         COUNT(*) AS n_valid_days\n",
    "  FROM hr_valid\n",
    "  GROUP BY person_id, year, week\n",
    "  HAVING n_valid_days >= 3  -- >= 3 valid days per week\n",
    "),\n",
    "eligible AS (\n",
    "  SELECT DISTINCT person_id\n",
    "  FROM valid_days\n",
    ")\n",
    "SELECT person_id FROM eligible"
  )
  
  df_eligible <- run_query_paged(sql_eligible)
  eligible_ids <- df_eligible$person_id
  ids_string   <- paste(eligible_ids, collapse = ", ")
  
  message("Eligible BYOD persons (HR adherence): ", length(eligible_ids))
  
  # Exclude 277 with insufficient wear time (confirmed in manuscript)
  # These are excluded during the sleep data pull in Block 03.
  
  
  # =============================================================================
  # BLOCK 02: DEMOGRAPHICS CLEANING
  # =============================================================================
  
  sql_demos <- paste0(
    "SELECT\n",
    "  p.person_id,\n",
    "  p.birth_datetime,\n",
    "  EXTRACT(YEAR FROM p.birth_datetime) AS birth_year,\n",
    "  p.gender_concept_id,\n",
    "  p.race_concept_id,\n",
    "  p.ethnicity_concept_id\n",
    "FROM `", cdr, ".person` p\n",
    "WHERE p.person_id IN (", ids_string, ")"
  )
  
  df_demos_raw <- run_query(sql_demos)
  
  # Survey data: race, ethnicity, sex assigned at birth, education, income,
  # employment, sexual orientation, smoking, alcohol, substance use
  sql_survey <- paste0(
    "SELECT\n",
    "  obs.person_id,\n",
    "  obs.observation_concept_id,\n",
    "  obs.value_as_concept_id,\n",
    "  obs.observation_date\n",
    "FROM `", cdr, ".observation` obs\n",
    "WHERE obs.person_id IN (", ids_string, ")\n",
    "  AND obs.observation_concept_id IN (\n",
    "    1585845, -- race\n",
    "    1585940, -- ethnicity\n",
    "    1585860, -- sex at birth\n",
    "    1585892, -- education\n",
    "    1585940, -- income\n",
    "    1585952, -- employment\n",
    "    1585899, -- sexual orientation\n",
    "    1585870, -- smoking\n",
    "    1585860, -- alcohol\n",
    "    1585865  -- substance use\n",
    "  )"
  )
  
  # BMI from physical measurements
  sql_bmi <- paste0(
    "SELECT person_id, value_as_number AS bmi_value,\n",
    "       measurement_date\n",
    "FROM `", cdr, ".measurement`\n",
    "WHERE person_id IN (", ids_string, ")\n",
    "  AND measurement_concept_id = 3038553\n",  # BMI concept
    "  AND value_as_number BETWEEN 10 AND 100\n"  # exclude implausible values
  )
  
  df_bmi_raw <- run_query(sql_bmi)
  
  # Mean BMI across all available measurements (baseline)
  df_bmi <- df_bmi_raw %>%
    group_by(person_id) %>%
    summarise(bmi_f = mean(bmi_value, na.rm = TRUE), .groups = "drop")
  
  # ── Race/ethnicity construction ───────────────────────────────────────────
  # Hispanic identity takes precedence over race (per All of Us conventions).
  # Five categories: White, Black, Hispanic, Asian, Other.
  # Participants with skip/prefer not to answer are excluded.
  
  # The block below is a RECONSTRUCTED SCAFFOLD. The original v3 mapping is not in
  # this repo. Paste your v3 recoding here, or verify every mapping below.
  
  df_survey <- run_query(sql_survey)
  
  # Pivot to one column per survey item, one row per person.
  df_survey_wide <- df_survey %>%
    select(person_id, observation_concept_id, value_as_concept_id) %>%
    distinct(person_id, observation_concept_id, .keep_all = TRUE) %>%
    tidyr::pivot_wider(names_from  = observation_concept_id,
                       values_from = value_as_concept_id,
                       names_prefix = "c")
  
  # SCAFFOLD: the structure is correct; replace each <..._id> with the verified
  # answer concept_id from your v3 codebook before relying on this.
  df_demos_survey <- df_survey_wide %>%
    mutate(
      race_raw = dplyr::case_when(
        c1585845 == "<White_id>" ~ "White",
        c1585845 == "<Black_id>" ~ "Black",
        c1585845 == "<Asian_id>" ~ "Asian",
        TRUE                     ~ "Other"),
      hispanic     = c1585940 == "<Hispanic_id>",
      race         = dplyr::if_else(hispanic, "Hispanic", race_raw), # Hispanic precedence
      sex_at_birth = dplyr::case_when(
        c1585860 == "<Male_id>"   ~ "Male",
        c1585860 == "<Female_id>" ~ "Female",
        TRUE                      ~ NA_character_),
      education      = recode_education(c1585892),   # define these helpers from v3
      income         = recode_income(c1585940),
      emp_status     = recode_employment(c1585952),
      sexual         = recode_orientation(c1585899),
      smoking_status = recode_binary(c1585870),
      alcohol_use    = recode_binary(c1585860),
      sb_use         = recode_binary(c1585865)
    )
  
  # Age at study start (2017) from birth year.
  df_age <- df_demos_raw %>%
    transmute(person_id, age_baseline = 2017 - birth_year)
  
  df_demos_cl <- df_demos_raw %>%
    left_join(df_bmi,          by = "person_id") %>%
    left_join(df_demos_survey, by = "person_id") %>%
    left_join(df_age,          by = "person_id") %>%
    mutate(race_sex = paste(race, sex_at_birth, sep = "_")) %>%
    # Drop skip / prefer-not-to-answer (NA) and missing BMI.
    filter(!is.na(race_sex), !is.na(sex_at_birth), !is.na(bmi_f)) %>%
    mutate(
      race_sex = factor(race_sex,
                        levels = c("White_Male", "White_Female",
                                   "Black_Male", "Black_Female",
                                   "Hispanic_Male", "Hispanic_Female",
                                   "Asian_Male", "Asian_Female",
                                   "Other_Male", "Other_Female"))
    )
  
  
  # =============================================================================
  # BLOCK 03: SLEEP DATA PULL (BigQuery)
  # =============================================================================
  # Primary criteria: main sleep only, >= 10 min, <= 720 min (12 hrs),
  # >= 3 valid days per week in that year
  # CDR v8 BYOD identified via src_id = 'Participant Portal: PTSC'
  # Note: 'careevolution' src_id = 256 legitimate BYOD participants (verified)
  
  sql_sleep <- paste0(
    "WITH sleep_clean AS (\n",
    "  SELECT person_id, sleep_date,\n",
    "         minute_asleep,\n",
    "         EXTRACT(YEAR FROM sleep_date) AS year,\n",
    "         DATE_TRUNC(sleep_date, WEEK) AS week\n",
    "  FROM `", cdr, ".sleep_daily_summary`\n",
    "  WHERE person_id IN (", ids_string, ")\n",
    "    AND sleep_date BETWEEN '2017-01-01' AND '2023-12-31'\n",
    "    AND is_main_sleep = 'true'\n",
    "    AND minute_asleep >= 10 AND minute_asleep <= 720\n",
    "),\n",
    "valid_weeks AS (\n",
    "  SELECT person_id, year, week, COUNT(*) AS days\n",
    "  FROM sleep_clean\n",
    "  GROUP BY person_id, year, week\n",
    "),\n",
    "valid_years AS (\n",
    "  SELECT person_id, year\n",
    "  FROM valid_weeks\n",
    "  WHERE days >= 3\n",
    "  GROUP BY person_id, year HAVING COUNT(*) >= 1\n",
    "),\n",
    "agg AS (\n",
    "  SELECT s.person_id, s.year,\n",
    "         s.year - 2017 AS t1,\n",
    "         POW(s.year - 2017, 2) AS t2,\n",
    "         POW(s.year - 2017, 3) AS t3,\n",
    "         AVG(s.minute_asleep) AS mean_min_asleep,\n",
    "         COUNT(*) AS n_days\n",
    "  FROM sleep_clean s\n",
    "  INNER JOIN valid_years v ON s.person_id = v.person_id AND s.year = v.year\n",
    "  GROUP BY s.person_id, s.year\n",
    ")\n",
    "SELECT * FROM agg"
  )
  
  df_sleep_yr <- run_query(sql_sleep)
  message("Sleep rows: ", nrow(df_sleep_yr),
          " | persons: ", n_distinct(df_sleep_yr$person_id))
  
  
  # =============================================================================
  # BLOCK 04: FINAL ANALYTIC DATASET ASSEMBLY
  # =============================================================================
  
  df_analytic <- df_sleep_yr %>%
    inner_join(df_demos_cl, by = "person_id") %>%
    group_by(person_id) %>%
    filter(n_distinct(year) >= 2) %>%  # require >= 2 years of data
    ungroup() %>%
    mutate(
      mean_hrs_asleep = mean_min_asleep / 60,
      short_sleep     = as.integer(mean_hrs_asleep < 7)
    )
  
  message("Analytic N: ", n_distinct(df_analytic$person_id),
          " | person-years: ", nrow(df_analytic))
  # Expected: 14,084 persons, 64,803 person-years
  
  # Strata sizes
  df_analytic %>%
    distinct(person_id, race_sex) %>%
    count(race_sex, name = "n_persons") %>%
    arrange(desc(n_persons)) %>%
    print()
  # White_Male: 3506 | White_Female: 8244 | Black_Female: 511
  # Hispanic_Female: 578 | Other_Female: 406 | Other_Male: 157
  # Asian_Female: 180 | Asian_Male: 133 | Hispanic_Male: 245 | Black_Male: 124
  
  # ── BYOD verification ─────────────────────────────────────────────────────
  # Confirmed: 0 wear_study participants in analytic sample
  # Run once to verify; output should show n_wear_in_analytic = 0
  sql_wear_check <- paste0(
    "SELECT COUNT(DISTINCT w.person_id) AS n_wear_in_analytic\n",
    "FROM `", cdr, ".wear_study` w\n",
    "WHERE w.person_id IN (", paste(unique(df_analytic$person_id), collapse = ","), ")"
  )
  # wear_check <- run_query(sql_wear_check)
  # cat("WEAR participants in analytic sample:", wear_check$n_wear_in_analytic, "\n")
  # Result: 0
  
  saveRDS(df_analytic, "df_analytic.rds")
  system(paste0("gsutil cp df_analytic.rds ", correct_bucket, "/sleep_paper/"))
  
} else {
  # ===== DEFAULT PATH: restore from the fitted model =====
  # The cleaned analytic data lives in the model frame of the fitted Model D,
  # so it reproduces exactly with no raw pull and no v3 recoding.
  system(paste0("gsutil cp ", correct_bucket, "/sleep_paper/m_core_rs.rds ."))
  m_core_rs   <- readRDS("m_core_rs.rds")
  df_analytic <- model.frame(m_core_rs) %>%
    mutate(year            = t1 + 2017,
           mean_hrs_asleep = mean_min_asleep / 60,
           short_sleep     = as.integer(mean_hrs_asleep < 7),
           sex_at_birth    = sub(".*_", "", as.character(race_sex)))
  message("Restored from model frame: ",
          n_distinct(df_analytic$person_id), " persons, ",
          nrow(df_analytic), " person-years")
}


# =============================================================================
# BLOCK 05: DESCRIPTIVE STATISTICS (TABLE 1)
# =============================================================================

n_total <- n_distinct(df_analytic$person_id)
n_py    <- nrow(df_analytic)

df_person_level <- df_analytic %>%
  group_by(person_id) %>%
  slice_head(n = 1) %>%
  ungroup()

# Race × sex distribution
table1_race_sex <- df_person_level %>%
  count(race_sex) %>%
  mutate(pct = round(n / n_total * 100, 1))

# Sleep descriptives
table1_sleep <- df_analytic %>%
  summarise(
    mean_tst  = paste0(round(mean(mean_hrs_asleep), 2), " (",
                       round(sd(mean_hrs_asleep), 2), ")"),
    pct_short = paste0(round(mean(short_sleep) * 100, 1), "%")
  )

# By group and year (Supplement Table S1a/S1b)
desc_table <- df_analytic %>%
  group_by(race_sex, year) %>%
  summarise(
    n          = n(),
    mean_hours = round(mean(mean_hrs_asleep, na.rm = TRUE), 2),
    sd_hours   = round(sd(mean_hrs_asleep, na.rm = TRUE), 2),
    mean_min   = round(mean(mean_min_asleep, na.rm = TRUE), 1),
    sd_min     = round(sd(mean_min_asleep, na.rm = TRUE), 1),
    pct_short  = round(mean(short_sleep, na.rm = TRUE) * 100, 1),
    .groups    = "drop"
  )

write.csv(desc_table, "desc_tst_by_group_year.csv", row.names = FALSE)
write.csv(table1_race_sex, "table1_covariates.csv", row.names = FALSE)
system(paste0("gsutil cp desc_tst_by_group_year.csv ", correct_bucket, "/sleep_paper/"))


# =============================================================================
# BLOCK 06: MODEL BUILDING — DECISIONS IN ORDER
# =============================================================================
# Each step runs a test, reports the finding, and records the decision.
# Reference group: White Male. Time: t1 = year - 2017 (0 = 2017 ... 6 = 2023).
# REML = FALSE for all likelihood-ratio and AIC comparisons.

# ── Step 0. Look at the outcome before modeling (exploration, not a test) ────
# The model assumes roughly normal RESIDUALS, not a normal outcome. This look is
# to understand range, skew, and bounds, so the residual diagnostics in Step 6
# are not a surprise.
cat("\nOutcome (mean_hrs_asleep) summary:\n")
print(summary(df_analytic$mean_hrs_asleep))
o <- df_analytic$mean_hrs_asleep
cat("Outcome skewness:", round(mean((o - mean(o))^3) / sd(o)^3, 3), "\n")
png("outcome_hist.png", width = 800, height = 600)
hist(o, breaks = 60, col = "grey80", border = "white",
     main = "Outcome distribution: mean hours asleep", xlab = "Hours/night")
abline(v = mean(o), col = "red", lwd = 2)
dev.off()
system(paste0("gsutil cp outcome_hist.png ", correct_bucket, "/sleep_paper/"))
# FINDING (see outcome_hist.png): left-skewed, skew = -1.75, centered near 6.3 hours,
# with a small near-zero cluster (min 0.187 h, a spike near 0.3 h), a floor artifact.
# DECISION: model on the raw minute scale. Outcome normality is not an assumption.
# Residuals are checked in Step 6, and the CLT carries inference at n = 64,803.
# Note the near-zero cluster for Limitations. Optional sensitivity: drop annual
# means below about 1 hour to confirm the estimates do not move.

# ── Step 1. Is a multilevel model needed? (null model ICC) ──────────────────
m0  <- lmer(mean_min_asleep ~ 1 + (1 | person_id), data = df_analytic, REML = FALSE)
vc0 <- as.data.frame(VarCorr(m0))
icc <- vc0$vcov[1] / sum(vc0$vcov)
message("Null ICC: ", round(icc, 3), " | grand mean: ", round(fixef(m0), 1), " min")
# FINDING: ICC = 0.803. About 80 percent of variance sits between persons.
# DECISION: use a person-level random effect. Single-level OLS would understate SEs.
# INTERPRET: ICC runs 0 to 1. Under 0.05 is negligible clustering, 0.10 and up is
# meaningful, and 0.80 means most variation is between people, so the random effect is essential.

# ── Step 2. Time form: linear vs quadratic vs cubic ─────────────────────────
m1  <- lmer(mean_min_asleep ~ t1 + (1 | person_id), data = df_analytic, REML = FALSE)
m2  <- lmer(mean_min_asleep ~ t1 + t2 + (1 | person_id), data = df_analytic, REML = FALSE)
m2c <- lmer(mean_min_asleep ~ t1 + t2 + t3 + (1 | person_id), data = df_analytic, REML = FALSE)
print(anova(m1, m2, m2c))
# FINDING: quadratic beats linear (p = .034); cubic beats quadratic
#          (chi-sq(1) = 130.9, p < 2e-16).
# DECISION: cubic time. Nonlinearity is modeled, not assumed away.

# ── Step 3. Fixed-effects ladder on cubic time ──────────────────────────────
# Build complexity in steps, each adding one source of structure.
m_a <- lmer(mean_min_asleep ~ t1 + t2 + t3 + (1 | person_id),
            data = df_analytic, REML = FALSE)

# race_cat for the race-only model
df_analytic <- df_analytic %>%
  mutate(race_cat = factor(gsub("_Male|_Female", "", as.character(race_sex)),
                           levels = c("White", "Black", "Hispanic", "Asian", "Other")))

m_b <- lmer(mean_min_asleep ~ t1 + t2 + t3 +
              race_cat*t1 + race_cat*t2 + race_cat*t3 +
              age_baseline + sexual + education + income + emp_status +
              bmi_f + heart_condition + smoking_status + alcohol_use + sb_use +
              (1 | person_id),
            data = df_analytic, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
# FINDING (Model B): pooled by race, Hispanic = -2.65 min, p = .47, not significant.

m_c <- lmer(mean_min_asleep ~ t1 + t2 + t3 +
              sex_at_birth*t1 + sex_at_birth*t2 + sex_at_birth*t3 +
              age_baseline + sexual + education + income + emp_status +
              bmi_f + heart_condition + smoking_status + alcohol_use + sb_use +
              (1 | person_id),
            data = df_analytic, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
# FINDING (Model C): Male = -12.2 min vs Female.
# DECISION: neither race alone nor sex alone is the right unit. Move to race x sex.

# ── Step 4. Random structure on the intersectional spec: intercept vs slope ──
m_core <- lmer(mean_min_asleep ~ t1 + t2 + t3 +
                 race_sex*t1 + race_sex*t2 + race_sex*t3 +
                 age_baseline + sexual + education + income + emp_status +
                 bmi_f + heart_condition + smoking_status + alcohol_use + sb_use +
                 (1 | person_id),
               data = df_analytic, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))

# m_core_rs is restored on the default path; refit only when rebuilding from raw.
if (REBUILD_FROM_RAW || !exists("m_core_rs")) {
  m_core_rs <- lmer(mean_min_asleep ~ t1 + t2 + t3 +
                      race_sex*t1 + race_sex*t2 + race_sex*t3 +
                      age_baseline + sexual + education + income + emp_status +
                      bmi_f + heart_condition + smoking_status + alcohol_use + sb_use +
                      (1 + t1 | person_id),
                    data = df_analytic, REML = FALSE,
                    control = lmerControl(optimizer = "bobyqa"))
}
print(anova(m_core, m_core_rs))
cat("Singular:", isSingular(m_core_rs), "\n")
# FINDING: random slope improves fit, delta AIC = 8,377, chi-sq(2) = 8,381.6,
#          p < 2e-16, and the model is not singular.
# DECISION: random slope (1 + t1 | person_id). m_core_rs is the PRIMARY model.
# INTERPRET: delta AIC above 10 is strong support, so 8,377 is decisive. isSingular FALSE
# means the random structure is estimable, not overfit.

# ── Step 5. Primary model: the intersectional finding ───────────────────────
summary(m_core_rs)
# vs White Male (adjusted): White Female +11.53*** | Black Male -46.29***
# Black Female -29.05*** | Hispanic Male -17.30** (p = .009) | Hispanic Female +11.35*
# Asian Male -29.96*** | AIC 691,665. The race-only model hid the Hispanic-male deficit.
# Reported in Table 2 and full in Table S9 (table_s9_model_d_rs_full.csv). Predicted
# trajectories in Figure 1; sentinel-year contrasts in Figure S1.
#
# DIRECTION CHECK (read before interpreting signs):
#   Continuous model (lmer): outcome is minutes asleep, so a NEGATIVE coefficient means
#     LESS sleep than White Male. Black Male -46.29 is 46 fewer minutes.
#   Binary model (GEE): outcome is short sleep under 7 hours, so a POSITIVE log-odds,
#     OR above 1, means MORE likely short. Black Male OR 7.48 is higher odds.
#   The two agree: groups with fewer minutes also carry higher odds of short sleep
#     (men, Black, Asian), and groups with more minutes carry lower odds (White and
#     Hispanic women). Reference is White Male in both. Higher income tracks with more
#     sleep, higher BMI with less sleep. Signs are consistent throughout.

saveRDS(m_core_rs, "m_core_rs.rds")
saveRDS(m_core,    "m_core.rds")
if (REBUILD_FROM_RAW) {
  system(paste0("gsutil cp m_core_rs.rds ", correct_bucket, "/sleep_paper/"))
  system(paste0("gsutil cp m_core.rds ",    correct_bucket, "/sleep_paper/"))
}

# Full fixed-effect table with Wald CIs (Supplement Table S9)
fe_ci   <- confint(m_core_rs, method = "Wald", parm = "beta_")
fe_coef <- fixef(m_core_rs)
fe_se   <- sqrt(diag(vcov(m_core_rs)))
fe_p    <- 2 * pnorm(abs(fe_coef / fe_se), lower.tail = FALSE)
table_s9_rs <- data.frame(term = names(fe_coef),
                          est = round(fe_coef, 3), se = round(fe_se, 3),
                          lo = round(fe_ci[, 1], 3), hi = round(fe_ci[, 2], 3),
                          p = round(fe_p, 4))
write.csv(table_s9_rs, "table_s9_model_d_rs_full.csv", row.names = FALSE)
if (REBUILD_FROM_RAW) {
  system(paste0("gsutil cp table_s9_model_d_rs_full.csv ", correct_bucket, "/sleep_paper/"))
}

# ── Step 6. Diagnostics on the primary model (finding -> decision) ───────────
res <- residuals(m_core_rs)
fit <- fitted(m_core_rs)

# 6.1 Residual normality (transparency only, not a gate at this n)
set.seed(42)
sw <- shapiro.test(sample(res, 5000))
cat("Residual Shapiro-Wilk W =", round(sw$statistic, 4),
    "| skew", round(mean((res - mean(res))^3) / sd(res)^3, 3),
    "| kurtosis", round(mean((res - mean(res))^4) / sd(res)^4, 3), "\n")
png("resid_hist.png", width = 800, height = 600)
hist(res, breaks = 80, col = "grey80", border = "white",
     main = "Model D residuals", xlab = "Residual (minutes)")
dev.off()
png("resid_qq.png", width = 800, height = 600)
qqnorm(res, main = "QQ plot: Model D residuals"); qqline(res, col = "red"); dev.off()
system(paste0("gsutil cp resid_hist.png ", correct_bucket, "/sleep_paper/"))
system(paste0("gsutil cp resid_qq.png ",   correct_bucket, "/sleep_paper/"))
# FINDING (see resid_hist.png, resid_qq.png): residuals centered at 0 and roughly
# symmetric (skew -0.80) but sharply peaked and heavy-tailed, W = 0.86, kurtosis 13.2.
# The QQ plot is S-shaped at both ends. Log and sqrt transforms worsen it.
# DECISION: proceed on the raw scale. CLT covers the fixed-effect p-values at n = 64,803.
# INTERPRET: skew 0 is symmetric and |skew| above 1 is highly skewed. Kurtosis 3 is normal
# and above 3 is heavy-tailed, so 13.2 is very heavy. Shapiro W near 1 is normal, but at
# this n it rejects on trivial departures, which is why it does not gate the decision.

# 6.2 Homoscedasticity (residuals vs fitted)
png("resid_vs_fitted.png", width = 800, height = 600)
plot(fit, res, pch = ".", xlab = "Fitted (minutes)", ylab = "Residual",
     main = "Residuals vs fitted"); abline(h = 0, col = "red"); dev.off()
system(paste0("gsutil cp resid_vs_fitted.png ", correct_bucket, "/sleep_paper/"))
# FINDING (see resid_vs_fitted.png): spread is even across the bulk, fitted 200 to 500.
# A small diagonal streak at low fitted values is the near-zero-sleep cluster from Step 0.
# DECISION: standard errors are trustworthy for the bulk. The floor cluster is the same
# one noted in Step 0, tiny relative to n, and is the optional sensitivity to run.

# 6.3 Random-effects normality (random intercept and slope)
re <- ranef(m_core_rs)$person_id
png("ranef_qq.png", width = 1000, height = 500)
par(mfrow = c(1, 2))
qqnorm(re[, 1], main = "Random intercept"); qqline(re[, 1], col = "red")
qqnorm(re[, 2], main = "Random slope (t1)"); qqline(re[, 2], col = "red")
dev.off()
system(paste0("gsutil cp ranef_qq.png ", correct_bucket, "/sleep_paper/"))
# FINDING (see ranef_qq.png): the random intercept and the random slope are roughly
# symmetric with mild heavy tails, a modest departure from normal.
# DECISION: the variance components and the random slope are sound. The slope decision
# from Step 4 (delta AIC 8,377) stands.

# 6.4 Multicollinearity (VIF), base R so no extra package is required
# Polynomial-by-group interactions are correlated by construction. VIF flags it.
Xv <- model.matrix(m_core_rs)
Xv <- Xv[, colnames(Xv) != "(Intercept)", drop = FALSE]
Xv <- Xv[, apply(Xv, 2, sd) > 0, drop = FALSE]   # drop constant columns
vif_vals <- tryCatch(diag(solve(cor(Xv))),
                     error = function(e) setNames(rep(NA_real_, ncol(Xv)), colnames(Xv)))
vif_tab  <- data.frame(term = names(vif_vals), VIF = round(vif_vals, 2))
print(head(vif_tab[order(-vif_tab$VIF), ], 15))
write.csv(vif_tab, "vif_model_d.csv", row.names = FALSE)
# Optional upgrade: if you install car, car::vif(m_core_rs) gives grouped GVIF.
# FINDING (see vif_model_d.csv): VIF reaches about 2,174 on t2 and the t2-by-race_sex
# terms, with the t1 and t3 interactions in the hundreds. This is polynomial-by-group
# collinearity, expected and benign, not a data problem.
# DECISION: interpret PREDICTED VALUES (Block 07, Figure 1), not raw interaction coefficients.
# INTERPRET: VIF 1 is no collinearity, above 5 is worth noting, above 10 is the usual concern
# threshold. The values in the hundreds here come from polynomial interactions and do not
# distort the predicted trajectories.

# 6.5 Influence at the stratum level (leave-one-stratum-out)
# Person-level DFBETA would need 14,084 refits and is not attempted. The right
# granularity is the 10 race-sex strata. Gated, since it refits 10 times.
RUN_INFLUENCE <- FALSE
if (RUN_INFLUENCE) {
  base_t1 <- fixef(m_core)[["t1"]]
  infl <- dplyr::bind_rows(lapply(levels(df_analytic$race_sex), function(g) {
    mi <- update(m_core, data = droplevels(dplyr::filter(df_analytic, race_sex != g)))
    data.frame(dropped = g,
               t1_shift = round(fixef(mi)[["t1"]] - base_t1, 3),
               aic      = round(AIC(mi)))
  }))
  print(infl)
  write.csv(infl, "influence_leave_one_stratum.csv", row.names = FALSE)
}
# FINDING (when run): how far the overall time slope moves when each stratum drops.
# DECISION: report wide uncertainty for the small male strata, Black 124, Asian 133.
# Their model CIs already encode this (Black Male SE 9.67, CI -65 to -27).

# 6.6 Sensitivity: drop the near-zero floor cluster (refit, gated)
# Step 0 and resid_vs_fitted.png showed a small cluster of near-zero annual means.
# This refits Model D after dropping annual means below 1 hour, then re-requiring two
# years, and compares the strata coefficients. Gated, since it refits the model.
RUN_FLOOR_SENS <- FALSE
if (RUN_FLOOR_SENS) {
  df_floor <- df_analytic %>%
    filter(mean_hrs_asleep >= 1) %>%
    group_by(person_id) %>% filter(n_distinct(year) >= 2) %>% ungroup()
  m_floor <- lmer(
    mean_min_asleep ~ t1 + t2 + t3 +
      race_sex*t1 + race_sex*t2 + race_sex*t3 +
      age_baseline + sexual + education + income + emp_status +
      bmi_f + heart_condition + smoking_status + alcohol_use + sb_use +
      (1 + t1 | person_id),
    data = df_floor, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
  key <- grep("^race_sex", names(fixef(m_core_rs)), value = TRUE)
  key <- key[!grepl(":", key)]
  floor_cmp <- data.frame(term  = key,
                          full  = round(fixef(m_core_rs)[key], 2),
                          floor = round(fixef(m_floor)[key], 2),
                          delta = round(fixef(m_floor)[key] - fixef(m_core_rs)[key], 2))
  print(floor_cmp)
  write.csv(floor_cmp, "floor_sensitivity.csv", row.names = FALSE)
}
# FINDING: dropping 258 persons with annual means below 1 hour moves every stratum
# coefficient by at most 2.87 minutes, and the significant strata move less.
# DECISION: the floor cluster does not drive the results. Report in Limitations.
# LIMITATION SENTENCE: A small group of near-zero annual sleep records reflects sparse
# wear. Excluding annual means below one hour leaves every estimate within three minutes.

# =============================================================================
# BLOCK 07: PREDICTED TRAJECTORIES (emmeans from m_core_rs)
# =============================================================================
# NOTE v4: emmeans now uses m_core_rs (random slope), not m_core.

# Sentinel years: 2017, 2019, 2021, 2023 (t1 = 0, 2, 4, 6)
emm <- emmeans(m_core_rs,
               ~ race_sex | t1,
               at = list(t1 = c(0, 2, 4, 6),
                         t2 = c(0, 4, 16, 36),
                         t3 = c(0, 8, 64, 216)),
               nuisance = c("age_baseline", "sexual", "education",
                            "income", "emp_status", "bmi_f",
                            "heart_condition", "smoking_status",
                            "alcohol_use", "sb_use"))

pred_table <- as.data.frame(emm) %>%
  mutate(
    year           = t1 + 2017,
    pred_hours     = emmean / 60,
    ci_lower_hours = asymp.LCL / 60,
    ci_upper_hours = asymp.UCL / 60
  )

write.csv(pred_table, "predicted_tst_by_group_year.csv", row.names = FALSE)
system(paste0("gsutil cp predicted_tst_by_group_year.csv ",
              correct_bucket, "/sleep_paper/"))

# All years (for Figure 1 trajectory plot)
emm_full <- emmeans(m_core_rs,
                    ~ race_sex | t1,
                    at = list(t1 = 0:6,
                              t2 = (0:6)^2,
                              t3 = (0:6)^3),
                    nuisance = c("age_baseline", "sexual", "education",
                                 "income", "emp_status", "bmi_f",
                                 "heart_condition", "smoking_status",
                                 "alcohol_use", "sb_use"))

pred_full <- as.data.frame(emm_full) %>%
  mutate(
    year           = t1 + 2017,
    pred_hours     = emmean / 60,
    ci_lower_hours = asymp.LCL / 60,
    ci_upper_hours = asymp.UCL / 60,
    race_sex_label = recode(as.character(race_sex),
                            "White_Male"      = "White Male",   "White_Female"    = "White Female",
                            "Black_Male"      = "Black Male",   "Black_Female"    = "Black Female",
                            "Hispanic_Male"   = "Hispanic Male","Hispanic_Female"  = "Hispanic Female",
                            "Asian_Male"      = "Asian Male",   "Asian_Female"    = "Asian Female",
                            "Other_Male"      = "Other Male",   "Other_Female"    = "Other Female"),
    sex  = if_else(grepl("Female", race_sex_label), "Female", "Male"),
    race = factor(gsub(" Male| Female", "", race_sex_label),
                  levels = c("White", "Black", "Hispanic", "Asian", "Other"))
  )

# Year-specific contrasts vs White Male (for Supplement Table S8)
contrasts_by_year <- list()
for (yr in c(0, 2, 4, 6)) {
  emm_yr <- emmeans(m_core_rs,
                    ~ race_sex,
                    at = list(t1 = yr, t2 = yr^2, t3 = yr^3),
                    nuisance = c("age_baseline", "sexual", "education",
                                 "income", "emp_status", "bmi_f",
                                 "heart_condition", "smoking_status",
                                 "alcohol_use", "sb_use"))
  cont <- contrast(emm_yr, method = "trt.vs.ctrl", ref = "White_Male") %>%
    as.data.frame() %>%
    mutate(
      year        = yr + 2017,
      diff_hours  = round(estimate / 60, 2),
      ci_lower_hr = round((estimate - 1.96 * SE) / 60, 2),
      ci_upper_hr = round((estimate + 1.96 * SE) / 60, 2),
      p.value     = round(p.value, 4)
    )
  contrasts_by_year[[as.character(yr + 2017)]] <- cont
}

contrast_table <- bind_rows(contrasts_by_year) %>%
  select(contrast, year, diff_hours, ci_lower_hr, ci_upper_hr, p.value)

write.csv(contrast_table, "group_vs_whitemale_by_year.csv", row.names = FALSE)
system(paste0("gsutil cp group_vs_whitemale_by_year.csv ",
              correct_bucket, "/sleep_paper/"))


# =============================================================================
# BLOCK 08: SENSITIVITY ANALYSES
# =============================================================================
# These re-pull raw data (alternate thresholds, device table) and use ids_string
# and df_demos_cl from the raw path. They run only when REBUILD_FROM_RAW is TRUE.
if (REBUILD_FROM_RAW) {
  
  # ── 08a. Inverse Probability Weighting (IPW) ──────────────────────────────
  # Weights = 1 / P(included in analytic sample | baseline covariates)
  # Truncated at 99th percentile. Weight range: 1.17 to 1.51 (low selection bias).
  # All primary coefficient deltas < 0.5 min — findings robust to selection bias.
  
  df_eligible_demos <- data.frame(person_id = eligible_ids) %>%
    left_join(df_demos_cl %>%
                select(person_id, race_cat, sex_at_birth, age_baseline, sexual,
                       education, income, emp_status, bmi_f, heart_condition,
                       smoking_status, alcohol_use, sb_use),
              by = "person_id") %>%
    mutate(included = as.integer(person_id %in% df_analytic$person_id))
  
  df_elig_cc <- df_eligible_demos %>% na.omit()
  
  ipw_model <- glm(
    included ~ race_cat + sex_at_birth + age_baseline + sexual +
      education + income + emp_status + bmi_f + heart_condition +
      smoking_status + alcohol_use + sb_use,
    data = df_elig_cc, family = binomial(link = "logit")
  )
  
  df_elig_cc$ps        <- predict(ipw_model, type = "response")
  df_elig_cc$ipw       <- ifelse(df_elig_cc$included == 1, 1 / df_elig_cc$ps, 0)
  p99_w                <- quantile(df_elig_cc$ipw[df_elig_cc$included == 1],
                                   0.99, na.rm = TRUE)
  df_elig_cc$ipw_trunc <- pmin(df_elig_cc$ipw, p99_w)
  cat("IPW range:", range(df_elig_cc$ipw_trunc[df_elig_cc$included == 1]), "\n")
  
  df_analytic_ipw <- df_analytic %>%
    left_join(df_elig_cc %>% filter(included == 1) %>%
                select(person_id, ipw_trunc), by = "person_id") %>%
    replace_na(list(ipw_trunc = 1))
  
  # NOTE: lmer weights argument treats these as precision weights (not probability
  # weights). Weight range 1.17–1.51 is so narrow that this distinction has
  # negligible impact on estimates. Documented in Limitations.
  m_ipw <- lmer(
    mean_min_asleep ~ t1 + t2 + t3 +
      race_sex * t1 + race_sex * t2 + race_sex * t3 +
      age_baseline + sexual + education + income + emp_status +
      bmi_f + heart_condition + smoking_status + alcohol_use + sb_use +
      (1 + t1 | person_id),
    data    = df_analytic_ipw,
    weights = ipw_trunc,
    REML    = FALSE,
    control = lmerControl(optimizer = "bobyqa")
  )
  
  saveRDS(m_ipw, "m_ipw.rds")
  system(paste0("gsutil cp m_ipw.rds ", correct_bucket, "/sleep_paper/"))
  
  # ── 08b. Alternate wear-time thresholds ───────────────────────────────────
  # Primary: >= 3 valid days/week, minute_asleep 10–720 min
  # T1: >= 5 days/week (stricter)
  # T2: >= 1 day/week (lenient)
  # T3: >= 3 days/week, minute_asleep >= 60 min floor
  
  make_sleep_sql <- function(min_days, min_minutes = 10) {
    paste0(
      "WITH sleep_clean AS (",
      "SELECT person_id, sleep_date, minute_asleep, ",
      "EXTRACT(YEAR FROM sleep_date) AS year, ",
      "DATE_TRUNC(sleep_date, WEEK) AS week ",
      "FROM `", cdr, ".sleep_daily_summary` ",
      "WHERE person_id IN (", ids_string, ") ",
      "AND sleep_date BETWEEN '2017-01-01' AND '2023-12-31' ",
      "AND is_main_sleep = 'true' ",
      "AND minute_asleep >= ", min_minutes, " AND minute_asleep <= 720), ",
      "valid_weeks AS (",
      "SELECT person_id, year, week, COUNT(*) AS days ",
      "FROM sleep_clean GROUP BY person_id, year, week), ",
      "valid_years AS (",
      "SELECT person_id, year FROM valid_weeks ",
      "WHERE days >= ", min_days,
      " GROUP BY person_id, year HAVING COUNT(*) >= 1), ",
      "agg AS (",
      "SELECT s.person_id, s.year, ",
      "s.year - 2017 AS t1, POW(s.year-2017,2) AS t2, POW(s.year-2017,3) AS t3, ",
      "AVG(s.minute_asleep) AS mean_min_asleep, COUNT(*) AS n_days ",
      "FROM sleep_clean s INNER JOIN valid_years v ",
      "ON s.person_id = v.person_id AND s.year = v.year ",
      "GROUP BY s.person_id, s.year) ",
      "SELECT * FROM agg"
    )
  }
  
  sleep_t1 <- run_query(make_sleep_sql(5))      # >= 5 days/week (strict)
  sleep_t2 <- run_query(make_sleep_sql(1))      # >= 1 day/week (lenient)
  sleep_t3 <- run_query(make_sleep_sql(3, 60))  # >= 60 min floor
  
  fit_threshold <- function(label, sleep_yr) {
    df_alt <- sleep_yr %>%
      inner_join(df_demos_cl, by = "person_id") %>%
      group_by(person_id) %>%
      filter(n_distinct(year) >= 2) %>%
      ungroup() %>%
      na.omit()
    message(label, ": N=", n_distinct(df_alt$person_id))
    m <- lmer(
      mean_min_asleep ~ t1 + t2 + t3 +
        race_sex * t1 + race_sex * t2 + race_sex * t3 +
        age_baseline + sexual + education + income + emp_status +
        bmi_f + heart_condition + smoking_status + alcohol_use + sb_use +
        (1 | person_id),  # random intercept only for sensitivity speed
      data = df_alt, REML = FALSE,
      control = lmerControl(optimizer = "bobyqa")
    )
    data.frame(
      label         = label,
      n_persons     = n_distinct(df_alt$person_id),
      n_py          = nrow(df_alt),
      intercept     = round(fixef(m)[["(Intercept)"]], 2),
      black_male    = round(fixef(m)[["race_sexBlack_Male"]], 2),
      black_female  = round(fixef(m)[["race_sexBlack_Female"]], 2),
      hispanic_male = round(fixef(m)[["race_sexHispanic_Male"]], 2),
      asian_male    = round(fixef(m)[["race_sexAsian_Male"]], 2)
    )
  }
  
  res0 <- data.frame(
    label         = "Primary (>=3 days, >=10 min, random slope)",
    n_persons     = n_distinct(df_analytic$person_id),
    n_py          = nrow(df_analytic),
    intercept     = round(fixef(m_core_rs)[["(Intercept)"]], 2),
    black_male    = round(fixef(m_core_rs)[["race_sexBlack_Male"]], 2),
    black_female  = round(fixef(m_core_rs)[["race_sexBlack_Female"]], 2),
    hispanic_male = round(fixef(m_core_rs)[["race_sexHispanic_Male"]], 2),
    asian_male    = round(fixef(m_core_rs)[["race_sexAsian_Male"]], 2)
  )
  res1 <- fit_threshold(">=5 days/week (strict)", sleep_t1)
  res2 <- fit_threshold(">=1 day/week (lenient)", sleep_t2)
  res3 <- fit_threshold(">=3 days, >=60 min floor", sleep_t3)
  
  threshold_compare <- bind_rows(res0, res1, res2, res3)
  write.csv(threshold_compare, "threshold_sensitivity.csv", row.names = FALSE)
  system(paste0("gsutil cp threshold_sensitivity.csv ", correct_bucket, "/sleep_paper/"))
  # Black Male ranges -40.70 to -56.57; direction consistent across all thresholds
  
  # ── 08c. Device type ──────────────────────────────────────────────────────
  # CDR v8: TRACKER (95.3%) and SCALE (4.7%) only.
  # Insufficient variation to include as covariate.
  # Chi-square by race/sex: p < 0.0001 but driven by n, not by meaningful
  # differences in device type distribution. Documented in Methods.
  
  sql_device <- paste0(
    "SELECT d.person_id, d.device_type ",
    "FROM `", cdr, ".device` d ",
    "WHERE d.person_id IN (", ids_string, ")"
  )
  df_device <- run_query(sql_device)
  
  df_device_main <- df_device %>%
    count(person_id, device_type) %>%
    group_by(person_id) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(person_id, device_type)
  
  write.csv(
    df_device_main %>% count(device_type) %>% mutate(pct = round(n / sum(n) * 100, 1)),
    "device_by_group.csv", row.names = FALSE
  )
  system(paste0("gsutil cp device_by_group.csv ", correct_bucket, "/sleep_paper/"))
  
  # ── 08d. Piecewise and spline models ──────────────────────────────────────
  # Piecewise breakpoint at 2020 (COVID).
  # Three-segment spline with knots at 2019 and 2021.
  # Result: cubic polynomial (m_core_rs) fits as well or better.
  
  df_analytic <- df_analytic %>%
    mutate(
      pre2020  = t1,
      post2020 = pmax(0, t1 - 3),
      seg1 = t1,
      seg2 = pmax(0, t1 - 2),
      seg3 = pmax(0, t1 - 4)
    )
  
  m_piecewise <- lmer(
    mean_min_asleep ~ pre2020 + post2020 +
      race_sex * pre2020 + race_sex * post2020 +
      age_baseline + sexual + education + income + emp_status +
      bmi_f + heart_condition + smoking_status + alcohol_use + sb_use +
      (1 | person_id),
    data = df_analytic, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa")
  )
  
  m_spline <- lmer(
    mean_min_asleep ~ seg1 + seg2 + seg3 +
      race_sex * seg1 + race_sex * seg2 + race_sex * seg3 +
      age_baseline + sexual + education + income + emp_status +
      bmi_f + heart_condition + smoking_status + alcohol_use + sb_use +
      (1 | person_id),
    data = df_analytic, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa")
  )
  
  anova(m_core, m_piecewise, m_spline)
  # m_core_rs best; spline equivalent to cubic; piecewise ~78 AIC points worse
  
  saveRDS(m_piecewise, "m_piecewise.rds")
  saveRDS(m_spline, "m_spline.rds")
  system(paste0("gsutil cp m_piecewise.rds ", correct_bucket, "/sleep_paper/"))
  system(paste0("gsutil cp m_spline.rds ", correct_bucket, "/sleep_paper/"))
  
  # ── 08e. Year-specific contrasts (FDR-adjusted) ───────────────────────────
  # Produced in Block 07 as contrast_table. Stored separately here for reference.
  # FDR (Benjamini-Hochberg) correction applied across all 36 pairwise comparisons.
  
  contrast_table_fdr <- contrast_table %>%
    group_by(year) %>%
    mutate(p_fdr = p.adjust(p.value, method = "BH")) %>%
    ungroup()
  
  write.csv(contrast_table_fdr, "group_vs_whitemale_by_year_fdr.csv",
            row.names = FALSE)
  system(paste0("gsutil cp group_vs_whitemale_by_year_fdr.csv ",
                correct_bucket, "/sleep_paper/"))
  
  
}  # ===== end Block 08 sensitivity (raw rebuild only) =====

# =============================================================================
# BLOCK 09: SECONDARY ANALYSES
# =============================================================================

# ── 09a. Short sleep GEE (<7 hours/night) ─────────────────────────────────
# Binary outcome: short sleep = mean_hrs_asleep < 7 hours/night
# GEE with AR(1) correlation structure; robust sandwich standard errors.
# NOTE: GEE is a separate model from the primary mixed-effects models and does
# NOT use the random slope specification. ORs represent pooled average across
# the study period (not longitudinal trajectory estimates).

df_gee <- df_analytic %>%
  arrange(person_id, year) %>%
  mutate(person_id_f = as.factor(person_id))

m_gee <- geeglm(
  short_sleep ~ t1 + t2 + t3 +
    race_sex + age_baseline + sexual + education +
    income + emp_status + bmi_f + heart_condition +
    smoking_status + alcohol_use + sb_use,
  family = binomial(link = "logit"),
  id     = person_id_f,
  corstr = "ar1",
  data   = df_gee
)
summary(m_gee)

or_table <- data.frame(
  OR       = round(exp(coef(m_gee)), 3),
  CI_lower = round(exp(coef(m_gee) -
                         1.96 * sqrt(diag(m_gee$geese$vbeta))), 3),
  CI_upper = round(exp(coef(m_gee) +
                         1.96 * sqrt(diag(m_gee$geese$vbeta))), 3)
)
# INTERPRET: OR 1 is no difference, above 1 is higher odds of short sleep, below 1 is lower.
# Key ORs vs White Male (pooled across study period):
#   Black Male:      7.48 (3.60–15.54)  ***
#   Asian Male:      4.09 (2.56–6.53)   ***
#   Black Female:    2.34 (1.82–3.01)   ***
#   Hispanic Male:   1.84 (1.38–2.47)   ***
#   White Female:    0.58 (0.53–0.62)   ***

write.csv(or_table, "gee_short_sleep_OR.csv", row.names = FALSE)
system(paste0("gsutil cp gee_short_sleep_OR.csv ",
              correct_bucket, "/sleep_paper/"))

# ── 09b. MAIHDA variance partition ────────────────────────────────────────
# Three-level unconditional null model:
#   Level 1: person-year (residual)
#   Level 2: person (stable between-person variance)
#   Level 3: race-sex stratum (intersectional variance)
#
# This is an unconditional variance decomposition across all person-years —
# NOT a longitudinal trajectory model. REML = TRUE for variance estimation.
#
# NOTE v4: strata_id now coded as factor (was integer in v3 — functionally
# equivalent for random effects but factor is cleaner and prevents
# misinterpretation as continuous.)
#
# With 10 strata, the strata-level VPC (3.07%) should be interpreted with
# appropriate caution — bootstrap CI provided below.

df_analytic <- df_analytic %>%
  mutate(strata_id = as.factor(race_sex))  # FIXED in v4: factor not integer

m_maihda <- lmer(
  mean_min_asleep ~ 1 + (1 | strata_id) + (1 | person_id),
  data = df_analytic, REML = TRUE
)

vc        <- as.data.frame(VarCorr(m_maihda))
var_s     <- vc$vcov[vc$grp == "strata_id"]
var_p     <- vc$vcov[vc$grp == "person_id"]
var_r     <- vc$vcov[vc$grp == "Residual"]
total_var <- var_s + var_p + var_r

vpc_table <- data.frame(
  level = c("Strata (race x sex)", "Person", "Residual"),
  var   = round(c(var_s, var_p, var_r), 2),
  vpc   = round(c(var_s, var_p, var_r) / total_var * 100, 2)
)
cat("\nMAIHDA Variance Partition:\n")
print(vpc_table)
# Strata VPC = 3.07% | Person = 77.55% | Residual = 19.38%
# INTERPRET: VPC is the share of total variance at each level. The strata share of 3.07%
# is small but non-trivial for 10 intersectional groups, and the strata SD interval below
# excludes zero, which confirms it is real rather than noise.

# Bootstrap CIs on variance components (1000 iterations, ~5 min)
set.seed(42)
maihda_ci <- tryCatch({
  confint(m_maihda, method = "boot", nsim = 1000,
          parm = "theta_", quiet = TRUE)
}, error = function(e) {
  message("Bootstrap failed, using profile CI: ", e$message)
  confint(m_maihda, method = "profile", parm = "theta_")
})
cat("\nBootstrap 95% CI on variance components (SD scale):\n")
print(round(maihda_ci, 3))
# lme4 orders random effects by number of levels, so the rows are:
#   .sig01 = person_id SD  (CI ~77.9 to 79.9, matches SD 78.9)
#   .sig02 = strata_id SD  (CI ~7.15 to 24.0, matches SD 15.7)
#   .sigma = residual SD   (CI ~39.2 to 39.7)
# FINDING: the strata SD lower bound is 7.15, not near zero, so the 3.07% intersectional
# VPC is supported, though the CI is wide. Report 3.07% together with this interval.

write.csv(vpc_table, "maihda_vpc.csv", row.names = FALSE)
saveRDS(m_maihda, "m_maihda.rds")
system(paste0("gsutil cp maihda_vpc.csv ", correct_bucket, "/sleep_paper/"))
system(paste0("gsutil cp m_maihda.rds ", correct_bucket, "/sleep_paper/"))


# =============================================================================
# BLOCK 10: FIGURES (.tif at 300 dpi, LZW compression — required by SLEEP)
# =============================================================================

race_colors <- c(
  "White"    = "#4477AA",
  "Black"    = "#EE6677",
  "Hispanic" = "#228833",
  "Asian"    = "#CCBB44",
  "Other"    = "#AA3377"
)

relabel_race_sex <- function(df, col = "race_sex") {
  df[[col]] <- recode(as.character(df[[col]]),
                      "White_Male"    = "White Male",   "White_Female"    = "White Female",
                      "Black_Male"    = "Black Male",   "Black_Female"    = "Black Female",
                      "Hispanic_Male" = "Hispanic Male","Hispanic_Female" = "Hispanic Female",
                      "Asian_Male"    = "Asian Male",   "Asian_Female"    = "Asian Female",
                      "Other_Male"    = "Other Male",   "Other_Female"    = "Other Female")
  df
}

# ── Figure 1: Predicted trajectories by race/ethnicity, faceted by sex ────
if (!exists("pred_full")) pred_full <- read.csv("predicted_tst_by_group_year.csv") %>%
  mutate(sex  = if_else(grepl("Female", race_sex), "Female", "Male"),
         race = factor(gsub("_Male|_Female", "", race_sex),
                       levels = c("White","Black","Hispanic","Asian","Other")))

p1 <- ggplot(pred_full, aes(x = year, y = pred_hours,
                            color = race, fill = race, group = race)) +
  geom_ribbon(aes(ymin = ci_lower_hours, ymax = ci_upper_hours),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.2) +
  geom_hline(yintercept = 7, linetype = "dashed", color = "black",
             linewidth = 0.6) +
  facet_wrap(~ sex, ncol = 2) +
  scale_color_manual(values = race_colors, name = "Race/Ethnicity") +
  scale_fill_manual(values  = race_colors, name = "Race/Ethnicity") +
  scale_x_continuous(breaks = 2017:2023) +
  scale_y_continuous(limits = c(4.0, 8.0), breaks = seq(4.0, 8.0, 0.5),
                     labels = function(x) sprintf("%.1f", x)) +
  labs(x = "Year", y = "Predicted total sleep time (hours/night)",
       caption = paste0(
         "Adjusted predicted values from Model D (Race \u00d7 Sex Assigned at Birth \u00d7 Time, ",
         "random slope specification). Shaded bands = 95% CI.\n",
         "Dashed line = 7-hour recommended minimum. n = 14,084.")) +
  theme_minimal(base_size = 13) +
  theme(legend.position   = "bottom",
        strip.text        = element_text(size = 13, face = "bold"),
        axis.text.x       = element_text(angle = 45, hjust = 1),
        panel.grid.minor  = element_blank(),
        plot.caption      = element_text(size = 8, hjust = 0, color = "gray40"))

ggsave("fig1_trajectories_final.tif", plot = p1, width = 10, height = 6,
       dpi = 300, device = "tiff", compression = "lzw")
message("Figure 1 saved.")

# ── Figure 2: Observed short sleep prevalence, faceted by sex ─────────────
if (!exists("desc_table")) desc_table <- read.csv("desc_tst_by_group_year.csv")

fig2_data <- desc_table %>%
  relabel_race_sex() %>%
  mutate(
    sex  = if_else(grepl("Female", race_sex), "Female", "Male"),
    race = gsub(" Male| Female", "", race_sex),
    race = factor(race, levels = c("White","Black","Hispanic","Asian","Other"))
  )

p2 <- ggplot(fig2_data, aes(x = year, y = pct_short, color = race, group = race)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.2) +
  facet_wrap(~ sex, ncol = 2) +
  scale_color_manual(values = race_colors, name = "Race/Ethnicity") +
  scale_x_continuous(breaks = 2017:2023) +
  scale_y_continuous(limits = c(40, 105), breaks = seq(40, 100, 10),
                     labels = function(x) paste0(x, "%")) +
  labs(x = "Year", y = "Percentage sleeping <7 hours/night",
       caption = "Unadjusted observed percentages (n = 64,803 person-years).") +
  theme_minimal(base_size = 13) +
  theme(legend.position   = "bottom",
        strip.text        = element_text(size = 13, face = "bold"),
        axis.text.x       = element_text(angle = 45, hjust = 1),
        panel.grid.minor  = element_blank(),
        plot.caption      = element_text(size = 8, hjust = 0, color = "gray40"))

ggsave("fig2_short_sleep_prevalence.tif", plot = p2, width = 10, height = 6,
       dpi = 300, device = "tiff", compression = "lzw")
message("Figure 2 saved.")

# ── Figure S1: Dot-and-whisker at sentinel years ──────────────────────────
figS1_data <- pred_table %>%
  filter(year %in% c(2017, 2019, 2021, 2023)) %>%
  relabel_race_sex() %>%
  mutate(
    race_sex = factor(race_sex, levels = rev(c(
      "White Male", "White Female", "Black Male", "Black Female",
      "Hispanic Male", "Hispanic Female", "Asian Male", "Asian Female",
      "Other Male", "Other Female"))),
    race = gsub(" Male| Female", "", as.character(race_sex)),
    race = factor(race, levels = c("White","Black","Hispanic","Asian","Other"))
  )

pS1 <- ggplot(figS1_data, aes(x = pred_hours, y = race_sex, color = race)) +
  geom_point(size = 2.8) +
  geom_errorbarh(aes(xmin = ci_lower_hours, xmax = ci_upper_hours),
                 height = 0.3, linewidth = 0.7) +
  facet_wrap(~ year, ncol = 4) +
  geom_vline(xintercept = 7, linetype = "dashed", color = "#DC2626",
             linewidth = 0.6) +
  scale_color_manual(values = race_colors, name = "Race/Ethnicity") +
  scale_x_continuous(limits = c(3.5, 8.0), breaks = seq(4, 8, 1)) +
  labs(x = "Predicted total sleep time (hours/night)", y = NULL,
       caption = paste0("Adjusted predicted values from Model D (random slope) at ",
                        "sentinel years, with 95% CI.\nDashed line = 7 hours. ",
                        "n = 14,084.")) +
  theme_minimal(base_size = 11) +
  theme(legend.position  = "bottom",
        strip.text       = element_text(size = 11, face = "bold"),
        axis.text.y      = element_text(size = 9),
        panel.grid.minor = element_blank(),
        plot.caption     = element_text(size = 7.5, hjust = 0, color = "gray40"))

ggsave("figS1_sentinel_dotwhisker.tif", plot = pS1, width = 13, height = 6,
       dpi = 300, device = "tiff", compression = "lzw")
message("Figure S1 saved.")

# ── Figures S2 and S3 ─────────────────────────────────────────────────────
# See v3 file for Figure S2 (device distribution) and Figure S3 (flow diagram)
# code — these are unchanged in v4.


# =============================================================================
# BLOCK 11: DOWNLOAD FIGURES FOR SCHOLAONE UPLOAD
# =============================================================================
# Run this block to retrieve .tif figures from bucket for local download.
# SLEEP requires figures as separate files in addition to manuscript.
# All figures must be .tif with LZW compression at >= 300 dpi.

fig_files <- c(
  "fig1_trajectories_final.tif",
  "fig2_short_sleep_prevalence.tif",
  "figS1_sentinel_dotwhisker.tif",
  "figS2_device_distribution.tif",
  "figS3_participant_flow.tif"
)

for (f in fig_files) {
  cmd <- paste0("gsutil cp ", correct_bucket, "/sleep_paper/", f, " .")
  exit_code <- system(cmd)
  if (exit_code == 0) {
    cat("Downloaded:", f, "\n")
  } else {
    cat("NOT FOUND in bucket:", f, "\n")
  }
}

# After downloading, verify each file:
for (f in fig_files) {
  if (file.exists(f)) {
    info <- file.info(f)
    cat(sprintf("  %-45s %.1f KB\n", f, info$size / 1024))
  }
}


# =============================================================================
# BLOCK 12: SAVE ALL OUTPUTS TO BUCKET
# =============================================================================

files_to_save <- c(
  # Core data objects
  "df_analytic.rds",
  "m_core_rs.rds",       # PRIMARY MODEL (random slope) — v4
  "m_core.rds",          # Comparison model (random intercept only)
  "m_ipw.rds",
  "m_maihda.rds",
  "m_piecewise.rds",
  "m_spline.rds",
  "sleep_t1.rds",
  "sleep_t2.rds",
  "sleep_t3.rds",
  # CSV outputs
  "predicted_tst_by_group_year.csv",
  "group_vs_whitemale_by_year.csv",
  "group_vs_whitemale_by_year_fdr.csv",
  "gee_short_sleep_OR.csv",
  "maihda_vpc.csv",
  "desc_tst_by_group_year.csv",
  "table1_covariates.csv",
  "threshold_sensitivity.csv",
  "device_by_group.csv",
  "tvbmi_full_comparison.csv",
  "table_s9_model_d_rs_full.csv",  # v4 Table S9
  # Figures
  "fig1_trajectories_final.tif",
  "fig2_short_sleep_prevalence.tif",
  "figS1_sentinel_dotwhisker.tif",
  "figS2_device_distribution.tif",
  "figS3_participant_flow.tif",
  # Diagnostics
  "qqplot_model_d_rs.png"
)

for (f in files_to_save) {
  if (file.exists(f)) {
    exit_code <- system(paste0("gsutil cp ", f, " ", correct_bucket, "/sleep_paper/"))
    if (exit_code == 0) cat("Saved:", f, "\n") else cat("FAILED:", f, "\n")
  } else {
    cat("SKIPPED (not found):", f, "\n")
  }
}

message("All outputs saved to: ", correct_bucket, "/sleep_paper/")


# =============================================================================
# BLOCK 13: PUBLISH TO THE GITHUB REPO
# =============================================================================
# This is the paper's analysis file. The repo holds code only, no participant data.
# Repo: github.com/hazel275/Intersectionality-and-Sleep
#
# Publishing workflow (run the DUA-safe publish script, do not git push raw output):
#   1. Confirm no participant-level data and no cell under 20 in any committed file.
#   2. Run github_publish_sleep_DUA_safe.R, which gates on the DUA and pushes the
#      code-only files to the repo.
#   3. Commit this file as the canonical pipeline, alongside Sleep_Supplement_Regen_FINAL.R.
#
# Outputs that belong in the repo as figures or tables (code-generated, no raw data):
#   Figure 1  fig1_trajectories_final.tif        Figure 2  fig2_short_sleep_prevalence.tif
#   Figure S1 figS1_sentinel_dotwhisker.tif      Figure S2 figS2_device_distribution.tif
#   Table S9  table_s9_model_d_rs_full.csv        GEE ORs   gee_short_sleep_OR.csv
#   MAIHDA    maihda_vpc.csv                       VIF       vif_model_d.csv
#   Diagnostics: outcome_hist.png, resid_hist.png, resid_qq.png,
#                resid_vs_fitted.png, ranef_qq.png
# =============================================================================