# ==============================================================================
# Core functions for compass variable diagnostics and adjusted DID
# ==============================================================================

# This file contains:
# 1. data-processing helpers,
# 2. conditional local independence diagnostics using tetrad tests,
# 3. OLS, DID, relevance diagnostics, and adjusted DID estimation.
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Data-processing helper functions
# ------------------------------------------------------------------------------

clean_likert_1_5 <- function(x) {
  x <- as.numeric(x)
  ifelse(x %in% 1:5, x, NA_real_)
}

reverse_likert_1_5 <- function(x) {
  x <- clean_likert_1_5(x)
  ifelse(is.na(x), NA_real_, 6 - x)
}

row_mean_min <- function(data, vars, min_items = 1) {
  stopifnot(is.data.frame(data))

  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing variables: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  x <- data[, vars, drop = FALSE]
  x <- as.data.frame(lapply(x, clean_likert_1_5))

  valid_n <- rowSums(!is.na(x))
  out <- rowMeans(x, na.rm = TRUE)
  out[valid_n < min_items] <- NA_real_
  out
}

standardize_in_sample <- function(x) {
  as.vector(scale(x))
}


# ------------------------------------------------------------------------------
# 2. Utilities for model summaries
# ------------------------------------------------------------------------------

extract_lm_effect <- function(fit, term = "T", model_name, conf_level = 0.95) {
  if (!inherits(fit, "lm")) {
    stop("fit must be an lm object.", call. = FALSE)
  }

  coef_mat <- coef(summary(fit))
  if (!term %in% rownames(coef_mat)) {
    stop("The term '", term, "' was not found in the model.", call. = FALSE)
  }

  ci <- confint(fit, parm = term, level = conf_level)

  data.frame(
    model = model_name,
    term = term,
    estimate = unname(coef_mat[term, "Estimate"]),
    se = unname(coef_mat[term, "Std. Error"]),
    ci_lower = unname(ci[1]),
    ci_upper = unname(ci[2]),
    p_value = unname(coef_mat[term, "Pr(>|t|)"]),
    n = stats::nobs(fit),
    row.names = NULL
  )
}


# ------------------------------------------------------------------------------
# 3. Relevance-condition diagnostics
# ------------------------------------------------------------------------------

diagnose_relevance <- function(data, detector, detector_label = detector,
                               pre = "P", treatment = "T", conf_level = 0.95) {
  required <- c(pre, treatment, detector)
  missing_vars <- setdiff(required, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing variables: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  f <- stats::reformulate(termlabels = c(detector, treatment), response = pre)
  fit <- stats::lm(f, data = data)

  coef_mat <- coef(summary(fit))
  ci <- confint(fit, parm = detector, level = conf_level)

  t_value <- coef_mat[detector, "t value"]
  df_resid <- df.residual(fit)
  partial_r2 <- t_value^2 / (t_value^2 + df_resid)

  data.frame(
    detector = detector_label,
    detector_variable = detector,
    b_P = unname(coef_mat[detector, "Estimate"]),
    se = unname(coef_mat[detector, "Std. Error"]),
    ci_lower = unname(ci[1]),
    ci_upper = unname(ci[2]),
    p_value = unname(coef_mat[detector, "Pr(>|t|)"]),
    partial_r2 = unname(partial_r2),
    model_r2 = summary(fit)$r.squared,
    n = stats::nobs(fit),
    row.names = NULL
  )
}


# ------------------------------------------------------------------------------
# 4. Conditional local independence diagnostics: tetrad test
# ------------------------------------------------------------------------------

test_tetrad_assumption <- function(data, outcome, pre, detector_1, detector_2,
                                   group, alpha = 0.05,
                                   adjust_method = "bonferroni") {
  vars_needed <- c(outcome, pre, detector_1, detector_2, group)
  missing_vars <- setdiff(vars_needed, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing variables: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  numeric_vars <- c(outcome, pre, detector_1, detector_2)
  non_numeric <- numeric_vars[!vapply(data[numeric_vars], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("The following variables must be numeric: ",
         paste(non_numeric, collapse = ", "), call. = FALSE)
  }

  groups <- sort(unique(stats::na.omit(data[[group]])))
  result_list <- vector("list", length(groups))

  for (i in seq_along(groups)) {
    g <- groups[i]
    sub_data <- data[data[[group]] == g, numeric_vars, drop = FALSE]
    sub_data <- stats::na.omit(sub_data)
    n <- nrow(sub_data)

    if (n < 5) {
      result_list[[i]] <- data.frame(
        group = g, n = n, tetrad = NA_real_, se = NA_real_,
        z_value = NA_real_, p_value = NA_real_
      )
      next
    }

    s <- stats::cov(sub_data)

    idx_p  <- which(colnames(s) == pre)
    idx_y  <- which(colnames(s) == outcome)
    idx_c1 <- which(colnames(s) == detector_1)
    idx_c2 <- which(colnames(s) == detector_2)

    s_pc1 <- s[idx_p, idx_c1]
    s_yc2 <- s[idx_y, idx_c2]
    s_pc2 <- s[idx_p, idx_c2]
    s_yc1 <- s[idx_y, idx_c1]

    tetrad <- s_pc1 * s_yc2 - s_pc2 * s_yc1

    gradient <- c(s_yc2, s_pc1, -s_yc1, -s_pc2)
    pairs_idx <- list(
      c(idx_p, idx_c1),
      c(idx_y, idx_c2),
      c(idx_p, idx_c2),
      c(idx_y, idx_c1)
    )

    cov_of_covs <- matrix(0, 4, 4)

    for (r in 1:4) {
      for (c in 1:4) {
        a <- pairs_idx[[r]][1]
        b <- pairs_idx[[r]][2]
        cc <- pairs_idx[[c]][1]
        d <- pairs_idx[[c]][2]

        cov_of_covs[r, c] <- (s[a, cc] * s[b, d] + s[a, d] * s[b, cc]) / (n - 1)
      }
    }

    var_tetrad <- as.numeric(t(gradient) %*% cov_of_covs %*% gradient)
    se <- sqrt(var_tetrad)

    if (is.na(se) || se <= 0) {
      z_value <- NA_real_
      p_value <- NA_real_
    } else {
      z_value <- tetrad / se
      p_value <- 2 * (1 - stats::pnorm(abs(z_value)))
    }

    result_list[[i]] <- data.frame(
      group = g,
      n = n,
      tetrad = tetrad,
      se = se,
      z_value = z_value,
      p_value = p_value
    )
  }

  result_table <- do.call(rbind, result_list)
  result_table$adjusted_p <- NA_real_

  valid_p <- !is.na(result_table$p_value)
  if (any(valid_p)) {
    result_table$adjusted_p[valid_p] <- stats::p.adjust(
      result_table$p_value[valid_p],
      method = adjust_method
    )
    reject <- any(result_table$adjusted_p[valid_p] < alpha)
    min_adjusted_p <- min(result_table$adjusted_p[valid_p])
  } else {
    reject <- NA
    min_adjusted_p <- NA_real_
  }

  output <- list(
    input = list(
      outcome = outcome,
      pre = pre,
      detectors = c(detector_1, detector_2),
      group = group,
      alpha = alpha,
      adjust_method = adjust_method
    ),
    result_table = result_table,
    decision = list(
      reject = reject,
      min_adjusted_p = min_adjusted_p
    )
  )

  class(output) <- "tetrad_test"
  output
}

print.tetrad_test <- function(x, digits = 4, ...) {
  cat("\n=== Tetrad Test for Conditional Local Independence ===\n")
  cat("Detectors:", paste(x$input$detectors, collapse = " and "), "\n")
  cat("Adjustment method:", x$input$adjust_method, "\n")
  cat("Alpha:", x$input$alpha, "\n\n")

  out <- x$result_table
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], round, digits)

  print(out, row.names = FALSE)

  cat("\nInterpretation:\n")
  if (isTRUE(x$decision$reject)) {
    cat("At least one adjusted p-value is below alpha. The tetrad implication is not supported.\n")
  } else if (isFALSE(x$decision$reject)) {
    cat("No adjusted p-value is below alpha. This does not prove the condition; it only indicates that this diagnostic did not find clear evidence against the tetrad implication.\n")
  } else {
    cat("The test could not be evaluated because valid p-values were not available.\n")
  }

  invisible(x)
}


# ------------------------------------------------------------------------------
# 5. Adjusted DID estimation
# ------------------------------------------------------------------------------

adjusted_did_analytic <- function(data, treatment = "T", outcome = "Y",
                                  pre = "P", detector = "C",
                                  detector_label = detector,
                                  conf_level = 0.95, verbose = TRUE) {
  required <- c(treatment, outcome, pre, detector)
  missing_vars <- setdiff(required, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing variables: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  d <- data[, required, drop = FALSE]
  d <- stats::na.omit(d)

  if (!all(d[[treatment]] %in% c(0, 1))) {
    stop("The treatment variable must be coded as 0/1.", call. = FALSE)
  }

  f_pre <- stats::reformulate(termlabels = c(detector, treatment), response = pre)
  f_out <- stats::reformulate(termlabels = c(detector, treatment), response = outcome)

  fit_pre <- stats::lm(f_pre, data = d, x = TRUE)
  fit_out <- stats::lm(f_out, data = d, x = TRUE)

  b_pre <- coef(fit_pre)[detector]
  g_pre <- coef(fit_pre)[treatment]
  b_out <- coef(fit_out)[detector]
  g_out <- coef(fit_out)[treatment]

  if (is.na(b_pre) || abs(b_pre) < .Machine$double.eps) {
    stop("The detector coefficient in the pretest regression is too close to zero.", call. = FALSE)
  }

  delta_hat <- as.numeric(b_out / b_pre)
  tau_hat <- as.numeric(g_out - delta_hat * g_pre)

  x <- model.matrix(fit_out)
  xtx_inv <- solve(crossprod(x))
  k <- ncol(x)
  n <- nrow(d)

  e_pre <- resid(fit_pre)
  e_out <- resid(fit_out)

  sigma_pp <- sum(e_pre^2) / (n - k)
  sigma_yy <- sum(e_out^2) / (n - k)
  sigma_py <- sum(e_pre * e_out) / (n - k)

  idx_c <- which(colnames(x) == detector)
  idx_t <- which(colnames(x) == treatment)

  inv_cc <- xtx_inv[idx_c, idx_c]
  inv_tt <- xtx_inv[idx_t, idx_t]
  inv_ct <- xtx_inv[idx_c, idx_t]

  # Parameter order: b_pre, b_out, g_pre, g_out
  v <- matrix(0, 4, 4)
  v[1, 1] <- sigma_pp * inv_cc
  v[2, 2] <- sigma_yy * inv_cc
  v[3, 3] <- sigma_pp * inv_tt
  v[4, 4] <- sigma_yy * inv_tt

  v[1, 2] <- v[2, 1] <- sigma_py * inv_cc
  v[3, 4] <- v[4, 3] <- sigma_py * inv_tt
  v[1, 3] <- v[3, 1] <- sigma_pp * inv_ct
  v[2, 4] <- v[4, 2] <- sigma_yy * inv_ct
  v[1, 4] <- v[4, 1] <- sigma_py * inv_ct
  v[2, 3] <- v[3, 2] <- sigma_py * inv_ct

  # tau = g_out - (b_out / b_pre) * g_pre
  grad <- c(
    g_pre * b_out / b_pre^2,
    -g_pre / b_pre,
    -b_out / b_pre,
    1
  )

  var_tau <- as.numeric(t(grad) %*% v %*% grad)
  se_tau <- sqrt(var_tau)

  z_stat <- tau_hat / se_tau
  p_value <- 2 * (1 - stats::pnorm(abs(z_stat)))
  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  ci <- tau_hat + c(-z_crit, z_crit) * se_tau

  result <- data.frame(
    model = "aDID",
    detector = detector_label,
    detector_variable = detector,
    estimate = tau_hat,
    se = se_tau,
    ci_lower = ci[1],
    ci_upper = ci[2],
    p_value = p_value,
    delta = delta_hat,
    b_P = as.numeric(b_pre),
    b_Y = as.numeric(b_out),
    g_P = as.numeric(g_pre),
    g_Y = as.numeric(g_out),
    n = n,
    row.names = NULL
  )

  if (isTRUE(verbose)) {
    print(result)
  }

  invisible(result)
}


# ------------------------------------------------------------------------------
# 6. Pairwise tetrad diagnostics wrapper
# ------------------------------------------------------------------------------

run_pairwise_tetrad_tests <- function(data, detectors, detector_labels = detectors,
                                      outcome = "Y", pre = "P", group = "T",
                                      alpha = 0.05, adjust_method = "bonferroni") {
  if (length(detectors) < 2) {
    stop("At least two detectors are required.", call. = FALSE)
  }

  pairs <- utils::combn(detectors, 2, simplify = FALSE)

  out <- lapply(pairs, function(pair) {
    test <- test_tetrad_assumption(
      data = data,
      outcome = outcome,
      pre = pre,
      detector_1 = pair[1],
      detector_2 = pair[2],
      group = group,
      alpha = alpha,
      adjust_method = adjust_method
    )

    label_1 <- detector_labels[match(pair[1], detectors)]
    label_2 <- detector_labels[match(pair[2], detectors)]

    cbind(
      detector_pair = paste(label_1, label_2, sep = "-"),
      detector_1 = pair[1],
      detector_2 = pair[2],
      test$result_table,
      row.names = NULL
    )
  })

  do.call(rbind, out)
}
