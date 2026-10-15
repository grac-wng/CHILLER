library(readxl)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)   # add p-values to lmer via Satterthwaite
library(emmeans)     # pairwise/contrast tests (ecto vs patho for H1)
library(broom.mixed)

# ----------------------------------------------------------------------------

raw <- read_excel(
  "data.xlsx",
  sheet = "Cleaned",
  col_names = FALSE
)

# first row is the header; rename columns
names(raw) <- as.character(unlist(raw[1, ]))
raw <- raw[-1, ]

names(raw)[1:4] <- c("participant_id", "session_id", "video_order", "measure")

# timepoint columns to numeric (Baseline + Column5...Column407)
time_cols <- names(raw)[5:ncol(raw)]
raw <- raw %>%
  mutate(across(all_of(time_cols), as.numeric)) %>%
  filter(!is.na(participant_id)) %>%  # drop trailing blank rows
  mutate(participant_id = as.integer(as.numeric(participant_id)))  # IDs are numbers in the sheet; keep them numeric

# reshape (one row per participant x measure x timepoint)
long_raw <- raw %>%
  pivot_longer(
    cols = all_of(time_cols),
    names_to = "timepoint_col",
    values_to = "value"
  ) %>%
  # recover a chronological timepoint index using column POSITION
  # use match() against time_cols instead
  mutate(
    timepoint_index = match(timepoint_col, time_cols) - 1L  # Baseline -> 0, next column -> 1, ...
  )

# pull "Time (seconds)" into its own column so every measure has a real time value
time_lookup <- long_raw %>%
  filter(measure == "Time (seconds)") %>%
  select(participant_id, video_order, timepoint_index, time_seconds = value)

long_data <- long_raw %>%
  filter(measure != "Time (seconds)") %>%
  left_join(time_lookup, by = c("participant_id", "video_order", "timepoint_index")) %>%
  pivot_wider(
    id_cols = c(participant_id, session_id, video_order, timepoint_index, time_seconds),
    names_from = measure,
    values_from = value
  ) %>%
  rename(
    pct_change_amplitude = `% change in amplitude`,
    spectral_units        = `spectral units`,
    goosebump_present     = `presence of goosebumps?`
  ) %>%
  arrange(participant_id, timepoint_index)

# long_data now has one row per participant per timepoint, with:
#   pct_change_amplitude (continuous), spectral_units, goosebump_present (0/1)
#   video_order = which sequence (A: ecto->wash->patho, B: patho->wash->ecto)
#   time_seconds = elapsed time since session start

# attach condition labels based on block timing
#   baseline / control : 0s   - 90s
#   block 1             : 90s  - 240s   (4 min)
#   washout             : 240s - 540s   (9 min)
#   block 2             : 540s - 690s   (11 min 30 sec)
# Block 1 = ectoparasite for Version A, pathogen for Version B (and vice versa
# for block 2), since the video order swaps which condition comes first.
# ----------------------------------------------------------------------------

baseline_end <- 90
block1_end   <- 240   # 4 min
washout_end  <- 540   # 9 min
block2_end   <- 690   # 11 min 30 sec

long_data <- long_data %>%
  mutate(
    condition = case_when(
      time_seconds <= baseline_end ~ "control",
      time_seconds > baseline_end & time_seconds <= block1_end &
        video_order == "Version A (ectoparasite -> wash-out -> pathogen)" ~ "ectoparasite",
      time_seconds > baseline_end & time_seconds <= block1_end &
        video_order == "Version B (pathogen -> wash-out -> ectoparasite)" ~ "pathogen",
      time_seconds > washout_end & time_seconds <= block2_end &
        video_order == "Version A (ectoparasite -> wash-out -> pathogen)" ~ "pathogen",
      time_seconds > washout_end & time_seconds <= block2_end &
        video_order == "Version B (pathogen -> wash-out -> ectoparasite)" ~ "ectoparasite",
      TRUE ~ "washout"
    )
  ) %>%
  filter(condition != "washout")

# collapse to one value per participant x condition
# peak % change during the block

safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

analysis_data <- long_data %>%
  group_by(participant_id, condition) %>%
  summarise(
    amplitude_pct_change = safe_max(pct_change_amplitude),
    goosebump_present    = safe_max(goosebump_present),  # 1 if present at any point in block
    .groups = "drop"
  ) %>%
  mutate(
    condition = factor(condition, levels = c("control", "pathogen", "ectoparasite"))
  )
analysis_data <- analysis_data %>%
  filter(!participant_id %in% c(42L, 78L, 169L))  # failed attention check (42, 78) / blank survey (169)

# how many participant x condition cells have no usable data 
analysis_data %>%
  summarise(
    missing_amplitude = sum(is.na(amplitude_pct_change)),
    missing_presence   = sum(is.na(goosebump_present))
  )


# H1: ectoparasite cues elicit greater piloerection than pathogen cues
# --- continuous DV: % change in goosebump amplitude ---
h1_amplitude_model <- lmer(
  amplitude_pct_change ~ condition + (1 | participant_id),
  data = analysis_data
)
summary(h1_amplitude_model)

# direct ectoparasite vs pathogen contrast 
h1_emm <- emmeans(h1_amplitude_model, ~ condition)
pairs(h1_emm)  # look at the ectoparasite - pathogen row

# --- binary DV: presence of goosebumps ---
h1_presence_model <- glmer(
  goosebump_present ~ condition + (1 | participant_id),
  data = analysis_data,
  family = binomial
)
summary(h1_presence_model)
pairs(emmeans(h1_presence_model, ~ condition))

# --- bias-reduced alternative (recommended) ---
# confirmed quasi-complete separation: control = 0.6% presence (1/168) vs.
# ~3-4% in the video conditions. plain glmer's fixed-effect estimates and
# p-values for a boundary this close to 0 aren't reliable 
# blme puts a weakly-informative t-prior on the fixed effects, which
# regularizes the near-boundary estimate without changing the model
# structure. compare these p-values/estimates to the plain glmer ones above.
library(blme)

h1_presence_model_br <- bglmer(
  goosebump_present ~ condition + (1 | participant_id),
  data = analysis_data,
  family = binomial,
  fixef.prior = t(df = 3)
)
summary(h1_presence_model_br)
pairs(emmeans(h1_presence_model_br, ~ condition))

# H2: both ectoparasite and pathogen cues elicit greater piloerection than
#     control

analysis_data <- analysis_data %>%
  mutate(condition = relevel(condition, ref = "control"))

# --- continuous DV ---
h2_amplitude_model <- lmer(
  amplitude_pct_change ~ condition + (1 | participant_id),
  data = analysis_data
)
summary(h2_amplitude_model)

# --- binary DV ---
h2_presence_model <- glmer(
  goosebump_present ~ condition + (1 | participant_id),
  data = analysis_data,
  family = binomial
)
summary(h2_presence_model)

# exploratory stuff 
#
survey_data <- readRDS("survey_analysis.rds") %>%
  mutate(participant_id = as.integer(as.numeric(participant_id)),
         hirsuteness = as.numeric(hirsuteness_self)) %>%   # RA rating isn't in the survey file
  select(participant_id, gender, pathogen_disgust, moral_disgust,
         outdoors_activities, natural_area_access, hirsuteness)

analysis_data <- analysis_data %>%
  left_join(survey_data, by = "participant_id")

# E1: gender differences in piloerection to ectoparasite/pathogen cues,
#     possibly moderated by disgust sensitivity

# gender main effect / interaction with condition
e1_amplitude_model <- lmer(
  amplitude_pct_change ~ condition * gender + (1 | participant_id),
  data = analysis_data
)
summary(e1_amplitude_model)

e1_presence_model <- bglmer(
  goosebump_present ~ condition * gender + (1 | participant_id),
  data = analysis_data,
  family = binomial,
  fixef.prior = t(df = 3)
)
summary(e1_presence_model)

# extended model adding disgust sensitivity as a moderator
e1b_amplitude_model <- lmer(
  amplitude_pct_change ~ condition * gender + condition * pathogen_disgust +
    condition * moral_disgust + (1 | participant_id),
  data = analysis_data
)
summary(e1b_amplitude_model)

e1b_presence_model <- bglmer(
  goosebump_present ~ condition * gender + condition * pathogen_disgust +
    condition * moral_disgust + (1 | participant_id),
  data = analysis_data,
  family = binomial,
  fixef.prior = t(df = 3)
)
summary(e1b_presence_model)


# E2: pathogen/moral disgust subscale scores predict piloerection

e2_amplitude_model <- lmer(
  amplitude_pct_change ~ condition * pathogen_disgust + condition * moral_disgust +
    (1 | participant_id),
  data = analysis_data
)
summary(e2_amplitude_model)

e2_presence_model <- bglmer(
  goosebump_present ~ condition * pathogen_disgust + condition * moral_disgust +
    (1 | participant_id),
  data = analysis_data,
  family = binomial,
  fixef.prior = t(df = 3)
)
summary(e2_presence_model)

# E3: prior outdoor-activity experience moderates sensitivity to the stimuli

analysis_data <- analysis_data %>%
  mutate(outdoors_activities = factor(outdoors_activities, ordered = TRUE))

e3_amplitude_model <- lmer(
  amplitude_pct_change ~ condition * outdoors_activities + (1 | participant_id),
  data = analysis_data
)
summary(e3_amplitude_model)

e3_presence_model <- bglmer(
  goosebump_present ~ condition * outdoors_activities + (1 | participant_id),
  data = analysis_data,
  family = binomial,
  fixef.prior = t(df = 3)
)
summary(e3_presence_model)

# E4: Childhood access to undeveloped natural areas moderates sensitivity

analysis_data <- analysis_data %>%
  mutate(natural_area_access = factor(natural_area_access, ordered = TRUE))

# cor(as.numeric(analysis_data$outdoors_activities),
#     as.numeric(analysis_data$natural_area_access), use = "complete.obs")

e4_amplitude_model <- lmer(
  amplitude_pct_change ~ condition * natural_area_access + (1 | participant_id),
  data = analysis_data
)
summary(e4_amplitude_model)

e4_presence_model <- bglmer(
  goosebump_present ~ condition * natural_area_access + (1 | participant_id),
  data = analysis_data,
  family = binomial,
  fixef.prior = t(df = 3)
)
summary(e4_presence_model)

# E5: Forearm hirsuteness affects the measured goosebump amplitude

e5_amplitude_model <- lmer(
  amplitude_pct_change ~ hirsuteness + (1 | participant_id),
  data = analysis_data
)
summary(e5_amplitude_model)

e5_presence_model <- bglmer(
  goosebump_present ~ hirsuteness + (1 | participant_id),
  data = analysis_data,
  family = binomial,
  fixef.prior = t(df = 3)
)
summary(e5_presence_model)


# E6: individual variation in disgust response 

e6_amplitude_model <- lmer(
  amplitude_pct_change ~ pathogen_disgust + moral_disgust + (1 | participant_id),
  data = analysis_data
)
summary(e6_amplitude_model)

e6_presence_model <- bglmer(
  goosebump_present ~ pathogen_disgust + moral_disgust + (1 | participant_id),
  data = analysis_data,
  family = binomial,
  fixef.prior = t(df = 3)
)
summary(e6_presence_model)

##############################################################################
# diagnostics (for every model above before trusting the p-values)
##############################################################################

# check_model_fit <- function(model) {
#   print(performance::check_singularity(model))
#   plot(residuals(model))
# }
# check_model_fit(h1_amplitude_model)
# check_model_fit(h2_amplitude_model)

##############################################################################
# tidy output tables (for the appendix / supplementary materials)
##############################################################################

broom.mixed::tidy(h1_amplitude_model, effects = "fixed", conf.int = TRUE)

# Use the bias-reduced model here, not h1_presence_model — plain glmer's
# estimates/CIs are unusable given the near-zero control base rate (see the
# ~300x-1200x odds ratios with CIs into the hundreds of thousands if you tidy
# h1_presence_model directly).
broom.mixed::tidy(h1_presence_model_br, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)


# tukey-adjusted CIs for the H2 tables

# --- table 2: amplitude vs. control (Tukey-adjusted) ---
h2_amplitude_contrasts <- confint(pairs(emmeans(h1_amplitude_model, ~ condition)))
print(h2_amplitude_contrasts)
#  pull the "control - pathogen" and "control - ectoparasite" rows
# "ectoparasite vs. control", to match the direction already used in the draft)

# --- table 3: presence vs. control (Tukey-adjusted, odds-ratio scale) ---
h2_presence_contrasts <- confint(pairs(emmeans(h1_presence_model_br, ~ condition)), type = "response")
print(h2_presence_contrasts)
# type = "response" back-transforms the log-odds CI to the odds-ratio scale
# directly, don't have to exponentiate by hand

# demographic stuff

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)

# read
raw <- read_excel(
  "chiller qualtrics.xlsx",
  sheet = "Sheet1",
  col_names = FALSE
)

var_names <- as.character(unlist(raw[1, ]))
question_text <- as.character(unlist(raw[2, ]))  # kept for reference, not used below
survey <- raw[-c(1, 2), ]
names(survey) <- var_names

# clean, recode

survey_clean <- survey %>%
  mutate(
    # --- gender ---
    gender = factor(Q1, levels = c("Male", "Female", "Non-binary", "Other")),
    
    # --- age from birth year ---
    # inconsistent formatting
    birth_year = as.integer(str_extract(Q2_1, "(19|20)\\d{2}")),
    age = 2026 - birth_year,  # survey collected Jan 2026 per StartDate
    
    # --- attention / data-quality checks ---
    # Q12 is an attention check ("look up on a clear day, what color do you
    # see?") with an expected answer of "Blue" — not a real demographic item.
    attention_check_passed = Q12 == "Blue",
    seriousness = factor(
      Q11,
      levels = c("Not at all seriously", "Somewhat seriously", "Quite seriously"),
      ordered = TRUE
    ),
    
    # --- hirsuteness (self-rating) ---
    hirsuteness_self = factor(
      Q7,
      levels = c("Sparse", "Moderate", "Dense"),
      ordered = TRUE
    ),
    # mean with RA rating? ...
    
    # --- outdoor activity history ---
    # trim trailing whitespace on read
    # normalize both sides before the factor match
    outdoors_activities = factor(
      str_squish(Q8),
      levels = c(
        "I do not participate in nature-based outdoor activities",
        "I rarely participate in nature-based outdoor activities",
        "I occasionally participate in nature-based outdoor activities",
        "I regularly participate in nature-based outdoor activities",
        "I frequently participate in nature-based outdoor activities"
      ),
      ordered = TRUE
    ),
    
    # --- childhood access to undeveloped natural areas ---
    natural_area_access = factor(
      str_squish(Q9),
      levels = c(
        "I grew up in an urban area with no undeveloped natural area in the immediate surroundings.",
        "I grew up in the suburbs of a city or a large town.",
        "I grew up in a rural area or on a farm, surrounded by croplands and livestock, with some undeveloped natural area available.",
        "I grew up in a sparsely populated area, surrounded by wilderness and undeveloped natural areas."
      ),
      ordered = TRUE
    )
  ) %>%
  mutate(across(
    c(Q4_1, Q4_2, Q4_3, Q4_4, Q4_5, Q5_1, Q5_2, Q5_3, Q5_4,
      Q6_1, Q6_2, Q6_3, Q6_4, Q6_5),
    ~ replace_na(as.numeric(.x), 0)
  ))

# 3 domains of disgust subscales (pathogen + moral only).
# item-to-subscale mapping identified from the question text:
#   pathogen: Q4_2 (dog poop), Q4_4 (red sores), Q5_1 (sweaty palms),
#             Q5_3 (mold), Q6_1 (body odor), Q6_3 (cockroach), Q6_5 (bloody cut)
#   moral:    Q4_1 (shoplifting), Q4_3 (stealing), Q4_5 (cheating),
#             Q5_2 (deceiving), Q5_4 (forging), Q6_2 (cutting in line),
#             Q6_4 (lying in business)
# each item is on a 0-6 scale

pathogen_items <- c("Q4_2", "Q4_4", "Q5_1", "Q5_3", "Q6_1", "Q6_3", "Q6_5")
moral_items    <- c("Q4_1", "Q4_3", "Q4_5", "Q5_2", "Q5_4", "Q6_2", "Q6_4")

survey_clean <- survey_clean %>%
  rowwise() %>%
  mutate(
    pathogen_disgust = sum(c_across(all_of(pathogen_items)), na.rm = FALSE),
    moral_disgust    = sum(c_across(all_of(moral_items)), na.rm = FALSE)
  ) %>%
  ungroup()

# data-quality summary

cat("N responses in file:", nrow(survey_clean), "\n")
cat("N passing attention check (Q12 = 'Blue'):",
    sum(survey_clean$attention_check_passed, na.rm = TRUE), "/",
    sum(!is.na(survey_clean$attention_check_passed)), "\n")

table(survey_clean$seriousness, useNA = "ifany")

# filter to attention-check-passed + "Quite/Somewhat seriously"
# survey_analysis <- survey_clean %>%
# --- exclusions: incomplete survey, failed attention checks, incomplete session ---
survey_items <- c(pathogen_items, moral_items)

crosswalk <- read_csv("participant_survey_crosswalk.csv") %>%
  select(participant_id, ResponseId)

survey_analysis <- survey_clean %>%
  inner_join(crosswalk, by = "ResponseId") %>%
  filter(attention_check_passed %in% TRUE)     # drops 42, 78, 169 -> N = 167

saveRDS(survey_analysis, "survey_analysis.rds")  # for the goosebump script

cat("Final N:", nrow(survey_analysis), "(expected 170)\n")
table(survey_analysis$gender, useNA = "ifany")
survey_analysis %>%
  summarise(mean_age = mean(age), sd_age = sd(age), min_age = min(age), max_age = max(age))

# check against the reported final sample (n = 167, women = 123, M = 20.01, SD = 1.57)
cat("Final N:", nrow(survey_analysis), "\n")
table(survey_analysis$gender, useNA = "ifany")
survey_analysis %>%
  summarise(mean_age = mean(age), sd_age = sd(age), min_age = min(age), max_age = max(age))

# demographic breakdown

# --- gender ---
survey_analysis %>%
  count(gender) %>%
  mutate(pct = round(100 * n / sum(n), 1))

# --- age ---
survey_analysis %>%
  summarise(
    n_valid   = sum(!is.na(age)),
    mean_age  = mean(age, na.rm = TRUE),
    sd_age    = sd(age, na.rm = TRUE),
    min_age   = min(age, na.rm = TRUE),
    max_age   = max(age, na.rm = TRUE)
  )

# --- hirsuteness (self-rating) ---
survey_analysis %>%
  count(hirsuteness_self) %>%
  mutate(pct = round(100 * n / sum(n), 1))

# --- outdoor activity history ---
survey_analysis %>%
  count(outdoors_activities) %>%
  mutate(pct = round(100 * n / sum(n), 1))

# --- childhood natural area access ---
survey_analysis %>%
  count(natural_area_access) %>%
  mutate(pct = round(100 * n / sum(n), 1))

# disgust subscale descriptives

survey_analysis %>%
  summarise(
    n                     = sum(!is.na(pathogen_disgust)),
    pathogen_mean         = mean(pathogen_disgust, na.rm = TRUE),
    pathogen_sd           = sd(pathogen_disgust, na.rm = TRUE),
    moral_mean            = mean(moral_disgust, na.rm = TRUE),
    moral_sd              = sd(moral_disgust, na.rm = TRUE)
  )

# by gender (E1 - gender differences in disgust sensitivity)
survey_analysis %>%
  filter(!is.na(gender)) %>%
  group_by(gender) %>%
  summarise(
    n             = n(),
    pathogen_mean = mean(pathogen_disgust, na.rm = TRUE),
    pathogen_sd   = sd(pathogen_disgust, na.rm = TRUE),
    moral_mean    = mean(moral_disgust, na.rm = TRUE),
    moral_sd      = sd(moral_disgust, na.rm = TRUE)
  )

# correlation between pathogen and moral disgust, and between outdoor
# activity and natural-area access (appendix E4 flags these two might need
# to be collapsed into one variable if highly correlated)
cor(survey_analysis$pathogen_disgust, survey_analysis$moral_disgust, use = "complete.obs")
cor(as.numeric(survey_analysis$outdoors_activities),
    as.numeric(survey_analysis$natural_area_access), use = "complete.obs")

# gender differences in disgust sensitivity
# non-binary (n = 1) too small for reliable group comparison, formal
# test just male vs. female

gender_comparison_data <- survey_analysis %>%
  filter(gender %in% c("Male", "Female")) %>%
  mutate(gender = droplevels(gender))

# welch's t-test (does not assume equal variance/n between groups, 
# unbalanced 43/123 split)
t.test(pathogen_disgust ~ gender, data = gender_comparison_data, var.equal = FALSE)
t.test(moral_disgust ~ gender, data = gender_comparison_data, var.equal = FALSE)

# effect size/cohen's d
library(effsize)
cohen.d(pathogen_disgust ~ gender, data = gender_comparison_data)
cohen.d(moral_disgust ~ gender, data = gender_comparison_data)
