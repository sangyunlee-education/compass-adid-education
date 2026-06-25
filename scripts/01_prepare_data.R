# ==============================================================================
# 01_prepare_data.R
# Construct the analytic sample for the empirical aDID example
# ==============================================================================

library(haven)
library(dplyr)
library(tidyr)

source("R/functions.R")

# ------------------------------------------------------------------------------
# 0. User settings
# ------------------------------------------------------------------------------

data_dir <- "data"
output_dir <- "output"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Y8 civic-consciousness items corresponding to the Y10 posttest outcome
pre_civic_vars <- c(
  "Y8H_ST19_1", "Y8H_ST19_2", "Y8H_ST19_3", "Y8H_ST19_5",
  "Y8H_ST19_6", "Y8H_ST19_7", "Y8H_ST19_8", "Y8H_ST19_10",
  "Y8H_ST19_11", "Y8H_ST19_15", "Y8H_ST19_16", "Y8H_ST19_17"
)

# Reverse-coded Y8 civic-consciousness items
pre_civic_reverse <- c(
  "Y8H_ST19_8", "Y8H_ST19_15", "Y8H_ST19_16", "Y8H_ST19_17"
)

# Compass-variable candidates: learning motivation
detector1_vars <- paste0("Y8H_ST24_", 9:12)   # identified regulation
detector2_vars <- paste0("Y8H_ST24_", 13:16)  # intrinsic motivation
detector3_vars <- paste0("Y8H_ST24_", 17:20)  # amotivation, later reverse-coded

# ------------------------------------------------------------------------------
# 1. Load raw data
# ------------------------------------------------------------------------------

d_y8 <- read_sav(file.path(data_dir, "y8STU.sav")) %>%
  rename_with(toupper)

d_y9 <- read_sav(file.path(data_dir, "y9STU.sav")) %>%
  rename_with(toupper)

d_y10 <- read_sav(file.path(data_dir, "Y10STU_학술대회.sav")) %>%
  rename_with(toupper)

# ------------------------------------------------------------------------------
# 2. Merge waves and define treatment
# ------------------------------------------------------------------------------

merged_base <- d_y9 %>%
  select(STUID, Y9H_ST5) %>%
  inner_join(
    d_y8 %>%
      select(
        STUID,
        Y8H_ST5,
        Y8H_ST19_1:Y8H_ST19_23,
        Y8H_ST24_9:Y8H_ST24_20
      ),
    by = "STUID"
  ) %>%
  inner_join(
    d_y10 %>%
      select(STUID, A3_1:A3_10),
    by = "STUID"
  ) %>%
  mutate(
    Y8H_ST5_NUM = as.numeric(Y8H_ST5),
    Y9H_ST5_NUM = as.numeric(Y9H_ST5),

    # Exclude students with part-time job experience at the pretest wave.
    pre_treated = if_else(Y8H_ST5_NUM == 1, 1, 0, missing = NA_real_),

    # Treatment: newly started part-time work by wave 9.
    # T = 1: no part-time work in wave 8, started by wave 9.
    # T = 0: no part-time work in both waves 8 and 9.
    T = case_when(
      Y9H_ST5_NUM == 1 ~ 1,
      Y9H_ST5_NUM == 2 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(pre_treated == 0)

# ------------------------------------------------------------------------------
# 3. Construct raw scores
# ------------------------------------------------------------------------------

analysis_raw <- merged_base %>%
  mutate(
    # Y10 posttest civic-consciousness items
    # Original scale: 1 = very important, 5 = not important at all.
    # Transformed scale: higher values indicate higher civic consciousness.
    across(A3_1:A3_10, reverse_likert_1_5),

    # Reverse-coded Y8 pretest civic-consciousness items.
    across(all_of(pre_civic_reverse), reverse_likert_1_5)
  ) %>%
  mutate(
    P_raw = row_mean_min(
      data = .,
      vars = pre_civic_vars,
      min_items = 8
    ),

    Y_raw = row_mean_min(
      data = .,
      vars = paste0("A3_", 1:10),
      min_items = 7
    ),

    D1_raw = row_mean_min(
      data = .,
      vars = detector1_vars,
      min_items = 3
    ),

    D2_raw = row_mean_min(
      data = .,
      vars = detector2_vars,
      min_items = 3
    ),

    # Reverse amotivation so that higher values indicate lower amotivation.
    D3_raw = 6 - row_mean_min(
      data = .,
      vars = detector3_vars,
      min_items = 3
    )
  )

# ------------------------------------------------------------------------------
# 4. Missing-data summary before complete-case restriction
# ------------------------------------------------------------------------------

missing_summary <- analysis_raw %>%
  summarise(
    n_before_complete_case = n(),
    missing_T = sum(is.na(T)),
    missing_P_raw = sum(is.na(P_raw)),
    missing_Y_raw = sum(is.na(Y_raw)),
    missing_D1_raw = sum(is.na(D1_raw)),
    missing_D2_raw = sum(is.na(D2_raw)),
    missing_D3_raw = sum(is.na(D3_raw))
  )

print("=== Missing-data summary before complete-case restriction ===")
print(missing_summary)

# ------------------------------------------------------------------------------
# 5. Complete-case analytic sample and within-sample standardization
# ------------------------------------------------------------------------------

final_data <- analysis_raw %>%
  select(STUID, T, P_raw, Y_raw, D1_raw, D2_raw, D3_raw) %>%
  drop_na(T, P_raw, Y_raw, D1_raw, D2_raw, D3_raw) %>%
  mutate(
    P = standardize_in_sample(P_raw),
    Y = standardize_in_sample(Y_raw),
    D1 = standardize_in_sample(D1_raw),
    D2 = standardize_in_sample(D2_raw),
    D3 = standardize_in_sample(D3_raw),
    diff_Y = Y - P
  ) %>%
  select(STUID, T, P, Y, diff_Y, D1, D2, D3,
         P_raw, Y_raw, D1_raw, D2_raw, D3_raw)

# ------------------------------------------------------------------------------
# 6. Safety checks
# ------------------------------------------------------------------------------

print("=== Final analytic sample size ===")
print(nrow(final_data))

print("=== Treatment-group distribution ===")
print(table(final_data$T, useNA = "ifany"))

stopifnot(all(stats::na.omit(final_data$T) %in% c(0, 1)))

range_check <- final_data %>%
  summarise(
    P_raw_min = min(P_raw, na.rm = TRUE),
    P_raw_max = max(P_raw, na.rm = TRUE),
    Y_raw_min = min(Y_raw, na.rm = TRUE),
    Y_raw_max = max(Y_raw, na.rm = TRUE),
    D1_raw_min = min(D1_raw, na.rm = TRUE),
    D1_raw_max = max(D1_raw, na.rm = TRUE),
    D2_raw_min = min(D2_raw, na.rm = TRUE),
    D2_raw_max = max(D2_raw, na.rm = TRUE),
    D3_raw_min = min(D3_raw, na.rm = TRUE),
    D3_raw_max = max(D3_raw, na.rm = TRUE)
  )

print("=== Raw-score range check ===")
print(range_check)

print("=== Correlation matrix for standardized variables ===")
print(round(cor(final_data %>% select(P, Y, D1, D2, D3), use = "complete.obs"), 3))

# ------------------------------------------------------------------------------
# 7. Save analytic data and summaries
# ------------------------------------------------------------------------------

saveRDS(final_data, file.path(output_dir, "final_data.rds"))
write.csv(missing_summary, file.path(output_dir, "missing_summary.csv"), row.names = FALSE)
write.csv(range_check, file.path(output_dir, "range_check.csv"), row.names = FALSE)

print("Saved output/final_data.rds")
