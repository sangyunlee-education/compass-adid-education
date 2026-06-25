# ==============================================================================
# 02_run_analysis.R
# Run OLS, DID, compass-variable diagnostics, and adjusted DID
# ==============================================================================

library(dplyr)

source("R/functions.R")

output_dir <- "output"

final_data_path <- file.path(output_dir, "final_data.rds")
if (!file.exists(final_data_path)) {
  stop("Run scripts/01_prepare_data.R first. Missing file: ", final_data_path, call. = FALSE)
}

final_data <- readRDS(final_data_path)

detectors <- c("D1", "D2", "D3")
detector_labels <- c("확인된 조절", "내재동기", "낮은 무동기")

# ------------------------------------------------------------------------------
# 1. OLS and conventional DID
# ------------------------------------------------------------------------------

fit_ols_unadjusted <- lm(Y ~ T, data = final_data)
fit_ols_pretest <- lm(Y ~ T + P, data = final_data)
fit_did <- lm(diff_Y ~ T, data = final_data)

basic_model_results <- bind_rows(
  extract_lm_effect(
    fit = fit_ols_unadjusted,
    term = "T",
    model_name = "OLS: no pretest adjustment"
  ),
  extract_lm_effect(
    fit = fit_ols_pretest,
    term = "T",
    model_name = "OLS: pretest-adjusted"
  ),
  extract_lm_effect(
    fit = fit_did,
    term = "T",
    model_name = "DID"
  )
)

print("=== OLS and conventional DID estimates ===")
print(basic_model_results)

# ------------------------------------------------------------------------------
# 2. Relevance-condition diagnostics
# ------------------------------------------------------------------------------

relevance_results <- bind_rows(
  diagnose_relevance(data = final_data, detector = "D1", detector_label = "확인된 조절"),
  diagnose_relevance(data = final_data, detector = "D2", detector_label = "내재동기"),
  diagnose_relevance(data = final_data, detector = "D3", detector_label = "낮은 무동기")
)

print("=== Relevance-condition diagnostics: P ~ C + T ===")
print(relevance_results)

# ------------------------------------------------------------------------------
# 3. Conditional local independence diagnostics
# ------------------------------------------------------------------------------

tetrad_results <- run_pairwise_tetrad_tests(
  data = final_data,
  detectors = detectors,
  detector_labels = detector_labels,
  outcome = "Y",
  pre = "P",
  group = "T",
  alpha = 0.05,
  adjust_method = "bonferroni"
)

print("=== Tetrad diagnostics for conditional local independence ===")
print(tetrad_results)

# ------------------------------------------------------------------------------
# 4. Adjusted DID estimates
# ------------------------------------------------------------------------------

adid_results <- bind_rows(
  adjusted_did_analytic(
    data = final_data,
    treatment = "T",
    outcome = "Y",
    pre = "P",
    detector = "D1",
    detector_label = "확인된 조절",
    verbose = FALSE
  ),
  adjusted_did_analytic(
    data = final_data,
    treatment = "T",
    outcome = "Y",
    pre = "P",
    detector = "D2",
    detector_label = "내재동기",
    verbose = FALSE
  ),
  adjusted_did_analytic(
    data = final_data,
    treatment = "T",
    outcome = "Y",
    pre = "P",
    detector = "D3",
    detector_label = "낮은 무동기",
    verbose = FALSE
  )
)

print("=== Adjusted DID estimates ===")
print(adid_results)

# ------------------------------------------------------------------------------
# 5. Combined table for paper
# ------------------------------------------------------------------------------

paper_table <- bind_rows(
  basic_model_results %>%
    transmute(
      model = model,
      covariate_or_detector = if_else(model == "OLS: pretest-adjusted", "P", "-"),
      estimate = estimate,
      se = se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = p_value
    ),
  adid_results %>%
    transmute(
      model = "aDID",
      covariate_or_detector = detector,
      estimate = estimate,
      se = se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = p_value
    )
)

print("=== Combined paper table ===")
print(paper_table)

# ------------------------------------------------------------------------------
# 6. Save results
# ------------------------------------------------------------------------------

write.csv(basic_model_results, file.path(output_dir, "basic_model_results.csv"), row.names = FALSE)
write.csv(relevance_results, file.path(output_dir, "relevance_results.csv"), row.names = FALSE)
write.csv(tetrad_results, file.path(output_dir, "tetrad_results.csv"), row.names = FALSE)
write.csv(adid_results, file.path(output_dir, "adid_results.csv"), row.names = FALSE)
write.csv(paper_table, file.path(output_dir, "paper_table.csv"), row.names = FALSE)

print("Analysis completed. Results saved in output/.")
