# Intersectionality and Sleep

Analysis code for a study of wearable-measured sleep disparities by race and sex in the All of Us Research Program. The data come from bring-your-own-device Fitbit records in CDR v8, covering 2017 to 2023, for 14,084 participants and 64,803 person-years. The work shows that race-by-sex strata carry sleep deficits that race-only models hide. Hispanic men, for example, show no significant deficit when race is pooled across sex, and a 17-minute deficit once the model separates race and sex.

## What this code does

The pipeline runs in a single annotated script. It pulls and cleans the data, builds the analytic sample, produces Table 1, and then builds the models as an ordered set of decisions. Each modeling step runs a test, reports the finding, and records the decision: the null model gives the intraclass correlation, a likelihood ratio test selects the cubic time form, and an AIC comparison selects the random-slope structure. The script then fits the primary intersectional model, runs the assumption diagnostics, and produces the secondary analyses and the figures.

## Data statement

This repository holds analysis code only. It contains no All of Us participant-level data. Any aggregate result follows the All of Us Data and Statistics Dissemination Policy, so no reported cell corresponds to fewer than 20 participants. Reproduce all results inside the All of Us Researcher Workbench under Controlled Tier access and your institution's data use agreement.

## Requirements

- R, run in RStudio rather than the Python kernel, because the mixed models use lme4.
- Packages: tidyverse, lme4, lmerTest, emmeans, geepack, broom.mixed, bigrquery.
- An All of Us account with Controlled Tier access.

## How to run

The main script has two paths, set by REBUILD_FROM_RAW near the top.

- Default, FALSE. The script restores the fitted primary model and reads the cleaned analytic data from its model frame. It runs end to end, reproduces the published results, and needs no raw data pull.
- Rebuild, TRUE. The script rebuilds the analytic sample from the CDR. This path needs the v3 demographics recoding pasted into Block 02, which is provided as a labeled scaffold.

Run the blocks in order, 00 through 13.

## Files

- Sleep_Analysis_CDRv8_Annotated.R. The complete, annotated pipeline. This is the canonical analysis file.

## Models

- Models A to D. Linear mixed growth models on total sleep time in minutes. Model D, the primary model, adds the race-by-sex strata interacted with time and a random slope.
- GEE. A population-averaged logistic model for short sleep, under 7 hours, with an AR(1) working correlation and robust standard errors.
- MAIHDA. A variance partition that places the race-by-sex strata as a level above the person.

## Citation

Cook SH, Holm J, Karaj A, Rodrigues M, Wood EP, Patippe C, Tian D. Wearable-measured sleep disparities by race and sex in the All of Us Research Program. Target journal: SLEEP.
