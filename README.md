Two datasets were used, linked by participant ID.

CHILLER time-series data (data.xlsx, sheet: "Cleaned"), wide format with one row per participant per measure:

participant_id: numeric session identifier (session start time, YYYYMMDDHHMM format); linking key to survey 

datasession_id: session identifier

video_order: counterbalancing condition ("Version A (ectoparasite -> wash-out -> pathogen)" or "Version B (pathogen -> wash-out -> ectoparasite)")

measure: one of three row types per participant: "Time (seconds)", "% change in amplitude", or "spectral units", or "presence of goosebumps?"

timepoint columns: one column per snapshot frame (Baseline, then sequential columns), each holding the value of the relevant measure at that timepoint



After reshaping to long format, the analysis dataset contained:

time_seconds: elapsed time in seconds since session start

pct_change_amplitude: goosebump amplitude as percentage change from individual baseline; derived from raw MATLAB output

spectral_units: raw goosebump spectral amplitude in arbitrary FFT power units, as output by the CHILLER device

goosebump_present: binary indicator of goosebump detection (1 = present, 0 = absent), as coded by the CHILLER algorithm (threshold: 30% change from baseline)

condition: experimental block label assigned based on time_seconds and video_order: "control/baseline" (0–90 s), "ectoparasite" or "pathogen" (Block 1: 90–240 s; Block 2: 540–690 s, counterbalanced by video order), or "washout" (240–540 s; excluded from analysis)


Qualtrics survey data (chiller qualtrics.xlsx, sheet: "Sheet1"):

gender: self-reported gender; factor with levels: Male, Female, Non-binary, Other (Q1)

birth_year: extracted from free-text birth year response (Q2_1); used to compute age

age: calculated as 2026 minus birth_year


attention_check_passed: logical; TRUE if participant answered "Blue" to the attention check item (Q12: "look up on a clear day, what color do you see?")

seriousness: self-reported survey engagement; ordered factor: "Not at all seriously" < "Somewhat seriously" < "Quite seriously" (Q11)

hirsuteness_self: self-rated forearm hair density; ordered factor: "Sparse" < "Moderate" < "Dense" (Q7)outdoors_activities: self-reported frequency of nature-based outdoor activities; ordered factor with five levels from "I do not participate" to "I frequently participate" (Q8)

natural_area_access: self-reported childhood access to undeveloped natural areas; ordered factor with four levels from urban/no access to sparsely populated/wilderness (Q9)

pathogen_disgust: sum of 7 pathogen-subscale items from the Three Domains of Disgust scale (Q4_2, Q4_4, Q5_1, Q5_3, Q6_1, Q6_3, Q6_5); each item scored 0–6; range 0–42

moral_disgust: sum of 7 moral-subscale items from the Three Domains of Disgust scale (Q4_1, Q4_3, Q4_5, Q5_2, Q5_4, Q6_2, Q6_4); each item scored 0–6; range 0–42


Analyses were conducted using R.

The R script read and reshaped the CHILLER time-series data from wide to long format; joined time-in-seconds from a separate measure row; assigned condition labels based on block timing and video order; collapsed to one row per participant per condition by taking the peak (maximum) value of pct_change_amplitude and goosebump_present within each block; joined survey covariates; and fit all statistical models (see below).

packages used:

readxl — Excel file import

dplyr, tidyr — data wrangling and reshaping

stringr, readr — string cleaning and CSV import

lme4 — linear mixed-effects models (lmer) and generalized linear mixed-effects models (glmer)

lmerTest — Satterthwaite-adjusted degrees of freedom and p-values for lmer models

blme — bias-reduced generalized linear mixed-effects models (bglmer) with weakly informative priors on fixed effects

emmeans — estimated marginal means and pairwise/Tukey-adjusted contrasts

broom.mixed — tidy model output tables

effsize — Cohen's d for gender comparisons on disgust subscales


MATLAB code used for real-time image acquisition and goosebump amplitude computation during data collection is also archived at the ripolleslab GitHub.

Raw CHILLER output was exported from MATLAB and cleaned in Excel, then imported into R. Qualtrics survey data were imported separately and cleaned: gender was recoded as a factor, age was computed as 2026 minus self-reported birth year, disgust subscale scores were computed by summing the relevant items (treating any skipped slider items as 0), and ordered factors were applied to hirsuteness, outdoor activity, and natural area access. Survey data were linked to CHILLER data via a participant crosswalk file matching ResponseId to participant_id.Of 200 initial participants, three were excluded prior to analysis: two (IDs 42, 78) failed the attention check and one (ID 169) had a blank survey response, yielding a final analytic sample of N = 167.

CHILLER data were reshaped from wide to long format (one row per participant per timepoint), with elapsed time in seconds recovered from a dedicated "Time (seconds)" measure row and joined to the remaining measures. Condition labels were assigned based on elapsed time and video order: the control period covered 0–90 s; Block 1 (90–240 s) and Block 2 (540–690 s) were labelled ectoparasite or pathogen according to the participant's counterbalancing version; the washout period (240–540 s) was excluded. Data were then collapsed to one observation per participant per condition by taking the peak (maximum) value of pct_change_amplitude and goosebump_present within each block, so that goosebump_present = 1 if a goosebump was detected at any point during the block.


Main analyses (H1, H2): Two dependent variables were modelled: peak percentage change in goosebump amplitude (continuous) and goosebump presence (binary). For each, a mixed-effects model was fit with condition (control, pathogen, ectoparasite) as the fixed effect and participant as a random intercept. Continuous models used lmer with Satterthwaite-adjusted degrees of freedom (lmerTest); binary models used bglmer (blme) with a weakly informative Student-t prior (df = 3) on fixed effects, chosen because the near-zero goosebump base rate in the control condition (1/168 observations, ~0.6%) produced quasi-complete separation that rendered plain glmer estimates unreliable. 

H1 (ectoparasite > pathogen) was tested via the ectoparasite − pathogen pairwise contrast from emmeans. 

H2 (both video conditions > control) was tested by relevelling the reference category to "control" and inspecting the pathogen and ectoparasite coefficients; Tukey-adjusted confidence intervals for all pairwise contrasts are reported in supplementary tables. Odds ratios for the binary model were obtained by back-transforming on the response scale (type = "response" in emmeans).


Exploratory analyses (E1–E6): Six exploratory models examined potential moderators of the piloerection response. 

E1 tested gender differences and the gender × condition interaction; E1b extended E1 by adding pathogen disgust × condition and moral disgust × condition interactions. 

E2 tested condition × disgust subscale interactions without gender. 

E3 and E4 tested condition × outdoor activity history and condition × childhood natural area access interactions, respectively, treating both variables as ordered factors. 

E5 tested whether self-rated forearm hirsuteness predicted goosebump amplitude and presence. 

E6 tested whether pathogen and moral disgust subscale scores predicted piloerection independent of condition. Gender comparisons on disgust subscale scores (reported as descriptive context) used Welch's two-sample t-tests (Male vs. Female only; the single Non-binary respondent was excluded from this comparison) with Cohen's d as the effect size estimate.
