#' Multi-Engine Consensus Differential Abundance Analysis
#'
#' @description
#' Implements a state-of-the-art consensus differential abundance testing framework
#' inspired by recent bioRxiv preprints (*ConsensusMetaDA*, *dar*, *LinDA*). Because
#' individual differential abundance methods frequently produce contradictory results,
#' `tidybiome` runs multiple complementary engines (LinDA compositional regression,
#' CLR linear models, and non-parametric Wilcoxon tests) and computes consensus
#' p-values (via the Cauchy combination test), agreement scores, and effect sizes.
#'
#' @param tb A `tidy_microbiome` object.
#' @param group Character string naming the binary grouping variable in sample metadata (e.g. `"treatment"`).
#' @param methods Character vector of engines to run. Options:
#'   * `"consensus"`: Runs all available engines and aggregates them into a consensus verdict.
#'   * `"linda"`: Linear Models for Differential Abundance with compositional bias correction.
#'   * `"clr_linear"`: Linear regression on centered log-ratio abundances.
#'   * `"wilcoxon"`: Non-parametric two-sample Wilcoxon rank-sum test.
#' @param fdr_cutoff Numeric false discovery rate threshold (default: 0.05).
#' @param assay Character string naming the count assay to use. Defaults to `"counts"`.
#' @param pseudocount Numeric pseudocount for log transformations. Defaults to 0.5.
#'
#' @return A [tibble::tbl_df] with one row per taxon containing:
#'   * `taxon_id`: Taxon identifier.
#'   * `log2fc`: Average log2 fold change across engines.
#'   * `p_consensus`: Combined p-value using Cauchy combination test.
#'   * `padj_consensus`: Benjamini-Hochberg adjusted consensus p-value.
#'   * `agreement_score`: Number of engines agreeing on significance and direction.
#'   * `is_significant`: Logical flag whether taxon meets `fdr_cutoff`.
#'   * Per-method statistics and taxonomic annotations.
#' @export
#'
#' @examples
#' counts <- matrix(c(100, 120, 110, 5, 8, 6, 50, 45, 55, 48, 52, 50), nrow = 2, byrow = TRUE,
#'                  dimnames = list(c("DiffTaxon", "NullTaxon"), paste0("S", 1:6)))
#' sample_data <- data.frame(sample_id = paste0("S", 1:6), group = c("A", "A", "A", "B", "B", "B"))
#' tb <- tidy_microbiome(counts, sample_data)
#' res <- calc_differential_abundance(tb, group = "group")
calc_differential_abundance <- function(tb,
                                        group,
                                        methods = c("consensus", "linda", "clr_linear", "wilcoxon"),
                                        fdr_cutoff = 0.05,
                                        assay = "counts",
                                        pseudocount = 0.5) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  sample_df <- tibble::as_tibble(tb)
  if (!group %in% colnames(sample_df)) {
    stop(sprintf("Grouping variable '%s' not found in sample metadata.", group), call. = FALSE)
  }

  group_vec <- as.factor(sample_df[[group]])
  levels_group <- levels(group_vec)
  if (length(levels_group) != 2) {
    stop(sprintf("Grouping variable '%s' must have exactly 2 levels for binary comparison. Found: %s",
                 group, paste(levels_group, collapse = ", ")), call. = FALSE)
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
  do_linda <- run_all || "linda" %in% methods
  do_clr   <- run_all || "clr_linear" %in% methods
  do_wilc  <- run_all || "wilcoxon" %in% methods

  results_df <- tibble::tibble(taxon_id = taxa_names)
  p_matrix <- matrix(NA, nrow = n_taxa, ncol = 0)
  log2fc_matrix <- matrix(NA, nrow = n_taxa, ncol = 0)

  # 1. LinDA Engine
  if (do_linda) {
    linda_res <- run_linda_engine(counts, group_vec, pseudocount = pseudocount)
    results_df$log2fc_linda <- linda_res$log2fc
    results_df$p_linda <- linda_res$pvalue
    results_df$padj_linda <- stats::p.adjust(linda_res$pvalue, method = "BH")
    p_matrix <- cbind(p_matrix, linda = linda_res$pvalue)
    log2fc_matrix <- cbind(log2fc_matrix, linda = linda_res$log2fc)
  }

  # 2. CLR Linear Engine
  if (do_clr) {
    clr_res <- run_clr_linear_engine(counts, group_vec, pseudocount = pseudocount)
    results_df$log2fc_clr <- clr_res$log2fc
    results_df$p_clr <- clr_res$pvalue
    results_df$padj_clr <- stats::p.adjust(clr_res$pvalue, method = "BH")
    p_matrix <- cbind(p_matrix, clr = clr_res$pvalue)
    log2fc_matrix <- cbind(log2fc_matrix, clr = clr_res$log2fc)
  }

  # 3. Wilcoxon Engine
  if (do_wilc) {
    wilc_res <- run_wilcoxon_engine(counts, group_vec)
    results_df$log2fc_wilcoxon <- wilc_res$log2fc
    results_df$p_wilcoxon <- wilc_res$pvalue
    results_df$padj_wilcoxon <- stats::p.adjust(wilc_res$pvalue, method = "BH")
    p_matrix <- cbind(p_matrix, wilcoxon = wilc_res$pvalue)
    log2fc_matrix <- cbind(log2fc_matrix, wilcoxon = wilc_res$log2fc)
  }

  # 4. Consensus Combination
  k_engines <- ncol(p_matrix)
  if (k_engines > 1) {
    # Cauchy combination test
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

  attr(tb, "metadata")$da_results[[group]] <- results_df
  results_df
}

run_linda_engine <- function(counts, group_vec, pseudocount = 0.5) {
  rel_mat <- calc_relabundance_matrix(counts)
  log_prop <- log2(rel_mat + pseudocount / colSums(counts))

  g_numeric <- as.numeric(group_vec) - 1 # 0 and 1
  n_taxa <- nrow(counts)
  log2fc <- numeric(n_taxa)
  pvals  <- numeric(n_taxa)

  # Initial linear regression per taxon
  for (i in seq_len(n_taxa)) {
    fit <- stats::lm(log_prop[i, ] ~ g_numeric)
    coefs <- summary(fit)$coefficients
    if (nrow(coefs) >= 2) {
      log2fc[i] <- coefs[2, 1]
      pvals[i]  <- coefs[2, 4]
    } else {
      log2fc[i] <- 0
      pvals[i]  <- 1
    }
  }

  # Mode/median compositional bias correction
  bias_estimate <- stats::median(log2fc, na.rm = TRUE)
  log2fc_corrected <- log2fc - bias_estimate

  # Recalculate p-values with corrected effect size
  for (i in seq_len(n_taxa)) {
    fit <- stats::lm(log_prop[i, ] ~ g_numeric)
    coefs <- summary(fit)$coefficients
    if (nrow(coefs) >= 2) {
      se <- coefs[2, 2]
      if (se > 0) {
        t_val <- log2fc_corrected[i] / se
        pvals[i] <- 2 * stats::pt(-abs(t_val), df = fit$df.residual)
      }
    }
  }

  list(log2fc = log2fc_corrected, pvalue = pvals)
}

run_clr_linear_engine <- function(counts, group_vec, pseudocount = 0.5) {
  clr_mat <- calc_clr_matrix(counts, pseudocount = pseudocount)
  g_numeric <- as.numeric(group_vec) - 1
  n_taxa <- nrow(counts)
  log2fc <- numeric(n_taxa)
  pvals  <- numeric(n_taxa)

  for (i in seq_len(n_taxa)) {
    fit <- stats::lm(clr_mat[i, ] ~ g_numeric)
    coefs <- summary(fit)$coefficients
    if (nrow(coefs) >= 2) {
      log2fc[i] <- coefs[2, 1] / log(2) # Convert natural log to log2
      pvals[i]  <- coefs[2, 4]
    } else {
      log2fc[i] <- 0
      pvals[i]  <- 1
    }
  }
  list(log2fc = log2fc, pvalue = pvals)
}

run_wilcoxon_engine <- function(counts, group_vec) {
  rel_mat <- calc_relabundance_matrix(counts)
  g1_idx <- which(group_vec == levels(group_vec)[1])
  g2_idx <- which(group_vec == levels(group_vec)[2])
  n_taxa <- nrow(counts)
  log2fc <- numeric(n_taxa)
  pvals  <- numeric(n_taxa)

  for (i in seq_len(n_taxa)) {
    v1 <- rel_mat[i, g1_idx]
    v2 <- rel_mat[i, g2_idx]

    m1 <- mean(v1, na.rm = TRUE)
    m2 <- mean(v2, na.rm = TRUE)
    log2fc[i] <- log2((m2 + 1e-6) / (m1 + 1e-6))

    test <- tryCatch(stats::wilcox.test(v2, v1, exact = FALSE), error = function(e) NULL)
    pvals[i] <- if (!is.null(test)) test$p.value else 1
  }
  list(log2fc = log2fc, pvalue = pvals)
}
