#' Multi-Engine Consensus Differential Abundance Analysis
#'
#' @description
#' Implements a state-of-the-art consensus differential abundance testing framework
#' inspired by recent bioRxiv preprints (*ConsensusMetaDA*, *dar*, *LinDA*). Supports
#' full formula specifications (`formula = ~ treatment + age + batch`) to control
#' for clinical covariates, continuous biomarkers (e.g. `formula = ~ bmi`), and multi-group
#' experimental designs. Runs multiple complementary engines (LinDA compositional regression,
#' CLR linear models, two-part CAFT, and rank/Wilcoxon tests) and computes consensus
#' p-values (via the Cauchy combination test), agreement scores, and effect sizes.
#'
#' @param tb A `tidy_microbiome` object.
#' @param group Optional character string naming the primary predictor column in sample metadata.
#' @param formula Optional model formula specifying the target predictor and covariates (e.g. `~ treatment + age + batch`).
#' @param covariates Optional character vector of covariate names to adjust for if `formula` is omitted.
#' @param contrast Optional character string naming the specific predictor term of interest in `formula`.
#'   Defaults to the first term specified in `formula` or `group`.
#' @param methods Character vector of engines to run. Options:
#'   * `"consensus"`: Runs all available engines and aggregates them into a consensus verdict.
#'   * `"caft"`: Compositional Log-Linear Model with Zero Cells (bioRxiv Dec 2025),
#'     jointly modeling presence/absence probability and conditional log-abundance.
#'   * `"linda"`: Linear Models for Differential Abundance with compositional bias correction.
#'   * `"clr_linear"`: Linear regression on centered log-ratio abundances with covariate adjustment.
#'   * `"wilcoxon"`: Non-parametric two-sample Wilcoxon rank-sum or rank-regression test.
#' @param fdr_cutoff Numeric false discovery rate threshold (default: 0.05).
#' @param assay Character string naming the count assay to use. Defaults to `"counts"`.
#' @param pseudocount Numeric pseudocount for log transformations. Defaults to 0.5.
#'
#' @return A [tibble::tbl_df] with one row per taxon containing:
#'   * `taxon_id`: Taxon identifier.
#'   * `log2fc`: Average log2 fold change or regression slope across engines.
#'   * `p_consensus`: Combined p-value using Cauchy combination test.
#'   * `padj_consensus`: Benjamini-Hochberg adjusted consensus p-value.
#'   * `agreement_score`: Number of engines agreeing on significance.
#'   * `is_significant`: Logical flag whether taxon meets `fdr_cutoff`.
#'   * Per-method statistics and taxonomic annotations.
#' @export
calc_differential_abundance <- function(tb,
                                        group = NULL,
                                        formula = NULL,
                                        covariates = NULL,
                                        contrast = NULL,
                                        methods = c("consensus", "caft", "linda", "clr_linear", "wilcoxon"),
                                        fdr_cutoff = 0.05,
                                        assay = "counts",
                                        pseudocount = 0.5) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  sample_df <- as.data.frame(tibble::as_tibble(tb))

  # Resolve formula / group / covariates
  if (is.null(formula)) {
    if (is.null(group)) {
      stop("Either `group` or `formula` must be specified.", call. = FALSE)
    }
    if (!group %in% colnames(sample_df)) {
      stop(sprintf("Grouping variable '%s' not found in sample metadata.", group), call. = FALSE)
    }
    cov_terms <- if (!is.null(covariates)) paste(covariates, collapse = " + ") else ""
    f_str <- if (nchar(cov_terms) > 0) paste("~", group, "+", cov_terms) else paste("~", group)
    formula <- stats::as.formula(f_str)
    target_var <- group
  } else {
    all_vars <- all.vars(formula)
    if (length(all_vars) == 0) {
      stop("`formula` must specify at least one predictor variable.", call. = FALSE)
    }
    missing_vars <- setdiff(all_vars, colnames(sample_df))
    if (length(missing_vars) > 0) {
      stop(sprintf("Variables not found in sample metadata: %s", paste(missing_vars, collapse = ", ")), call. = FALSE)
    }
    target_var <- if (!is.null(contrast)) contrast else all_vars[1]
    if (!target_var %in% all_vars) {
      stop(sprintf("Target contrast '%s' is not among the variables in formula.", target_var), call. = FALSE)
    }
  }

  covariate_vars <- setdiff(all.vars(formula), target_var)

  target_vals <- sample_df[[target_var]]
  is_continuous <- is.numeric(target_vals)

  if (!is_continuous) {
    sample_df[[target_var]] <- as.factor(target_vals)
  }

  assays <- attr(tb, "assays")
  if (!assay %in% names(assays)) {
    stop(sprintf("Assay '%s' not found in tidy_microbiome.", assay), call. = FALSE)
  }

  counts <- assays[[assay]]
  taxa_names <- rownames(counts)
  n_taxa <- length(taxa_names)

  # Check which engines to execute
  run_all <- "consensus" %in% methods
  do_caft  <- run_all || "caft" %in% methods
  do_linda <- run_all || "linda" %in% methods
  do_clr   <- run_all || "clr_linear" %in% methods
  do_wilc  <- run_all || "wilcoxon" %in% methods

  results_df <- tibble::tibble(taxon_id = taxa_names)
  p_matrix <- matrix(NA, nrow = n_taxa, ncol = 0)
  log2fc_matrix <- matrix(NA, nrow = n_taxa, ncol = 0)

  pred_df <- sample_df[, c(target_var, covariate_vars), drop = FALSE]

  # 1. CAFT Engine (Zero-cell compositional hurdle model)
  if (do_caft) {
    caft_res <- run_caft_engine(counts, pred_df, target_var, covariate_vars, is_continuous)
    results_df$log2fc_caft <- caft_res$log2fc
    results_df$p_caft <- caft_res$pvalue
    results_df$padj_caft <- stats::p.adjust(caft_res$pvalue, method = "BH")
    p_matrix <- cbind(p_matrix, caft = caft_res$pvalue)
    log2fc_matrix <- cbind(log2fc_matrix, caft = caft_res$log2fc)
  }

  # 2. LinDA Engine
  if (do_linda) {
    linda_res <- run_linda_engine(counts, pred_df, target_var, covariate_vars, is_continuous, pseudocount = pseudocount)
    results_df$log2fc_linda <- linda_res$log2fc
    results_df$p_linda <- linda_res$pvalue
    results_df$padj_linda <- stats::p.adjust(linda_res$pvalue, method = "BH")
    p_matrix <- cbind(p_matrix, linda = linda_res$pvalue)
    log2fc_matrix <- cbind(log2fc_matrix, linda = linda_res$log2fc)
  }

  # 3. CLR Linear Engine
  if (do_clr) {
    clr_res <- run_clr_linear_engine(counts, pred_df, target_var, covariate_vars, is_continuous, pseudocount = pseudocount)
    results_df$log2fc_clr <- clr_res$log2fc
    results_df$p_clr <- clr_res$pvalue
    results_df$padj_clr <- stats::p.adjust(clr_res$pvalue, method = "BH")
    p_matrix <- cbind(p_matrix, clr = clr_res$pvalue)
    log2fc_matrix <- cbind(log2fc_matrix, clr = clr_res$log2fc)
  }

  # 4. Wilcoxon / Rank Regression Engine
  if (do_wilc) {
    wilc_res <- run_wilcoxon_engine(counts, pred_df, target_var, covariate_vars, is_continuous)
    results_df$log2fc_wilcoxon <- wilc_res$log2fc
    results_df$p_wilcoxon <- wilc_res$pvalue
    results_df$padj_wilcoxon <- stats::p.adjust(wilc_res$pvalue, method = "BH")
    p_matrix <- cbind(p_matrix, wilcoxon = wilc_res$pvalue)
    log2fc_matrix <- cbind(log2fc_matrix, wilcoxon = wilc_res$log2fc)
  }

  # Consensus Combination
  k_engines <- ncol(p_matrix)
  if (k_engines > 1) {
    comb_p <- apply(p_matrix, 1, function(p_row) {
      valid_p <- p_row[!is.na(p_row)]
      if (length(valid_p) == 0) return(NA_real_)
      valid_p <- pmax(1e-15, pmin(1 - 1e-15, valid_p))
      t_stat <- mean(tan((0.5 - valid_p) * pi))
      p_cauchy <- 0.5 - (atan(t_stat) / pi)
      pmax(0, pmin(1, p_cauchy))
    })

    results_df$p_consensus <- comb_p
    results_df$padj_consensus <- stats::p.adjust(comb_p, method = "BH")
    results_df$log2fc <- rowMeans(log2fc_matrix, na.rm = TRUE)

    # Agreement score
    results_df$agreement_score <- apply(p_matrix, 1, function(p_row) {
      sum(p_row <= fdr_cutoff & !is.na(p_row))
    })
    results_df$is_significant <- results_df$padj_consensus <= fdr_cutoff & results_df$agreement_score >= 1
  } else {
    results_df$log2fc <- log2fc_matrix[, 1]
    results_df$p_consensus <- p_matrix[, 1]
    results_df$padj_consensus <- stats::p.adjust(p_matrix[, 1], method = "BH")
    results_df$agreement_score <- as.integer(p_matrix[, 1] <= fdr_cutoff)
    results_df$is_significant <- results_df$padj_consensus <= fdr_cutoff
  }

  # Merge with taxonomy
  tax_df <- attr(tb, "tax_table")
  if (!is.null(tax_df) && nrow(tax_df) > 0) {
    results_df <- merge(results_df, tax_df, by = "taxon_id", all.x = TRUE, sort = FALSE)
  }

  # Reorder by significance
  results_df <- results_df[order(results_df$padj_consensus, -abs(results_df$log2fc)), ]
  results_df <- tibble::as_tibble(results_df)

  attr(tb, "metadata")$da_results[[target_var]] <- results_df
  results_df
}

run_linda_engine <- function(counts, pred_df, target_var, covariate_vars, is_continuous, pseudocount = 0.5) {
  rel_mat <- calc_relabundance_matrix(counts)
  log_prop <- log2(rel_mat + pseudocount / colSums(counts))
  n_taxa <- nrow(counts)
  log2fc <- numeric(n_taxa)
  pvals  <- numeric(n_taxa)

  f_str <- paste("y ~", paste(c(target_var, covariate_vars), collapse = " + "))
  fit_formula <- stats::as.formula(f_str)
  fits <- vector("list", n_taxa)

  for (i in seq_len(n_taxa)) {
    pred_df$y <- log_prop[i, ]
    fit <- tryCatch(stats::lm(fit_formula, data = pred_df), error = function(e) NULL)
    fits[[i]] <- fit
    if (!is.null(fit)) {
      coefs <- summary(fit)$coefficients
      pattern <- paste0("^", target_var)
      target_idx <- grep(pattern, rownames(coefs))
      if (length(target_idx) >= 1) {
        log2fc[i] <- coefs[target_idx[1], 1]
        pvals[i]  <- coefs[target_idx[1], 4]
      } else {
        pvals[i] <- 1
      }
    } else {
      pvals[i] <- 1
    }
  }

  # Mode/median compositional bias correction
  bias_estimate <- stats::median(log2fc, na.rm = TRUE)
  log2fc_corrected <- log2fc - bias_estimate

  # Re-compute p-values with corrected effect size
  for (i in seq_len(n_taxa)) {
    fit <- fits[[i]]
    if (!is.null(fit)) {
      coefs <- summary(fit)$coefficients
      pattern <- paste0("^", target_var)
      target_idx <- grep(pattern, rownames(coefs))
      if (length(target_idx) >= 1) {
        se <- coefs[target_idx[1], 2]
        if (!is.na(se) && se > 0) {
          t_val <- log2fc_corrected[i] / se
          pvals[i] <- 2 * stats::pt(-abs(t_val), df = fit$df.residual)
        }
      }
    }
  }

  list(log2fc = log2fc_corrected, pvalue = pvals)
}

run_clr_linear_engine <- function(counts, pred_df, target_var, covariate_vars, is_continuous, pseudocount = 0.5) {
  clr_mat <- calc_clr_matrix(counts, pseudocount = pseudocount)
  n_taxa <- nrow(counts)
  log2fc <- numeric(n_taxa)
  pvals  <- numeric(n_taxa)

  f_str <- paste("y ~", paste(c(target_var, covariate_vars), collapse = " + "))
  fit_formula <- stats::as.formula(f_str)

  for (i in seq_len(n_taxa)) {
    pred_df$y <- clr_mat[i, ]
    fit <- tryCatch(stats::lm(fit_formula, data = pred_df), error = function(e) NULL)
    if (!is.null(fit)) {
      coefs <- summary(fit)$coefficients
      pattern <- paste0("^", target_var)
      target_idx <- grep(pattern, rownames(coefs))
      if (length(target_idx) == 1) {
        log2fc[i] <- coefs[target_idx, 1] / log(2)
        pvals[i]  <- coefs[target_idx, 4]
      } else if (length(target_idx) > 1) {
        log2fc[i] <- coefs[target_idx[which.max(abs(coefs[target_idx, 1]))], 1] / log(2)
        drop_res <- tryCatch(stats::drop1(fit, scope = stats::as.formula(paste("~", target_var)), test = "F"), error = function(e) NULL)
        pvals[i] <- if (!is.null(drop_res) && "Pr(>F)" %in% colnames(drop_res)) drop_res$`Pr(>F)`[2] else coefs[target_idx[1], 4]
      } else {
        log2fc[i] <- 0
        pvals[i] <- 1
      }
    } else {
      log2fc[i] <- 0
      pvals[i] <- 1
    }
  }
  list(log2fc = log2fc, pvalue = pvals)
}

run_wilcoxon_engine <- function(counts, pred_df, target_var, covariate_vars, is_continuous) {
  rel_mat <- calc_relabundance_matrix(counts)
  n_taxa <- nrow(counts)
  log2fc <- numeric(n_taxa)
  pvals  <- numeric(n_taxa)

  target_vals <- pred_df[[target_var]]
  has_covariates <- length(covariate_vars) > 0

  if (!is_continuous && length(unique(target_vals)) == 2 && !has_covariates) {
    g_fac <- as.factor(target_vals)
    lvls <- levels(g_fac)
    g1_idx <- which(g_fac == lvls[1])
    g2_idx <- which(g_fac == lvls[2])
    for (i in seq_len(n_taxa)) {
      v1 <- rel_mat[i, g1_idx]
      v2 <- rel_mat[i, g2_idx]
      m1 <- mean(v1, na.rm = TRUE)
      m2 <- mean(v2, na.rm = TRUE)
      log2fc[i] <- log2((m2 + 1e-6) / (m1 + 1e-6))
      wt <- tryCatch(stats::wilcox.test(v2, v1, exact = FALSE), error = function(e) NULL)
      pvals[i] <- if (!is.null(wt)) wt$p.value else 1
    }
  } else {
    for (i in seq_len(n_taxa)) {
      y <- rel_mat[i, ]
      pred_df$rank_y <- rank(y)
      f_str <- paste("rank_y ~", paste(c(target_var, covariate_vars), collapse = " + "))
      fit <- tryCatch(stats::lm(stats::as.formula(f_str), data = pred_df), error = function(e) NULL)
      if (!is.null(fit)) {
        coefs <- summary(fit)$coefficients
        pattern <- paste0("^", target_var)
        t_idx <- grep(pattern, rownames(coefs))
        if (length(t_idx) >= 1) {
          log2fc[i] <- coefs[t_idx[1], 1] / (n_taxa + 1)
          pvals[i]  <- coefs[t_idx[1], 4]
        } else {
          pvals[i] <- 1
        }
      } else {
        pvals[i] <- 1
      }
    }
  }
  list(log2fc = log2fc, pvalue = pvals)
}

# CAFT: Compositional Log-Linear Model with Zero Cells (bioRxiv Dec 2025)
run_caft_engine <- function(counts, pred_df, target_var, covariate_vars, is_continuous) {
  rel_mat <- calc_relabundance_matrix(counts)
  n_taxa <- nrow(counts)
  log2fc <- numeric(n_taxa)
  pvals  <- numeric(n_taxa)

  null_f <- if (length(covariate_vars) > 0) {
    paste("z ~", paste(covariate_vars, collapse = " + "))
  } else {
    "z ~ 1"
  }
  full_f <- paste("z ~", paste(c(target_var, covariate_vars), collapse = " + "))
  abund_f <- paste("log_y ~", paste(c(target_var, covariate_vars), collapse = " + "))

  for (i in seq_len(n_taxa)) {
    y <- rel_mat[i, ]
    z <- as.numeric(y > 0)
    pred_df$z <- z

    # Occurrence LRT
    w_zero <- 0
    if (length(unique(z)) == 2) {
      fit_null <- suppressWarnings(tryCatch(stats::glm(stats::as.formula(null_f), data = pred_df, family = stats::binomial()), error = function(e) NULL))
      fit_full <- suppressWarnings(tryCatch(stats::glm(stats::as.formula(full_f), data = pred_df, family = stats::binomial()), error = function(e) NULL))
      if (!is.null(fit_null) && !is.null(fit_full)) {
        lrt_val <- as.numeric(2 * (stats::logLik(fit_full) - stats::logLik(fit_null)))
        w_zero <- max(0, lrt_val)
      }
    }

    # Abundance part
    pos_idx <- which(z == 1)
    w_abund <- 0
    lfc_est <- 0
    if (length(pos_idx) >= 4) {
      sub_df <- pred_df[pos_idx, , drop = FALSE]
      sub_df$log_y <- log(y[pos_idx])
      fit_abund <- tryCatch(stats::lm(stats::as.formula(abund_f), data = sub_df), error = function(e) NULL)
      if (!is.null(fit_abund)) {
        c_abund <- summary(fit_abund)$coefficients
        pattern <- paste0("^", target_var)
        t_idx <- grep(pattern, rownames(c_abund))
        if (length(t_idx) >= 1 && !is.na(c_abund[t_idx[1], 2]) && c_abund[t_idx[1], 2] > 0) {
          w_abund <- (c_abund[t_idx[1], 1] / c_abund[t_idx[1], 2])^2
          lfc_est <- c_abund[t_idx[1], 1] / log(2)
        }
      }
    }

    if (lfc_est == 0 && !is_continuous && length(unique(pred_df[[target_var]])) == 2) {
      lvls <- levels(as.factor(pred_df[[target_var]]))
      m0 <- mean(y[pred_df[[target_var]] == lvls[1]], na.rm = TRUE)
      m1 <- mean(y[pred_df[[target_var]] == lvls[2]], na.rm = TRUE)
      lfc_est <- log2((m1 + 1e-6) / (m0 + 1e-6))
    }
    log2fc[i] <- lfc_est

    s_stat <- w_zero + w_abund
    df_test <- (w_zero > 0) + (w_abund > 0)
    if (df_test > 0) {
      pvals[i] <- stats::pchisq(s_stat, df = df_test, lower.tail = FALSE)
    } else {
      pvals[i] <- 1
    }
  }
  list(log2fc = log2fc, pvalue = pvals)
}
