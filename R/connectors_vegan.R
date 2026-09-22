#' Run PERMANOVA via vegan::adonis2
#'
#' Executes Permutational Multivariate Analysis of Variance using `vegan::adonis2`
#' on a `tidy_microbiome` object and returns a clean, tidy tibble.
#'
#' @param tb A `tidy_microbiome` object.
#' @param formula A model formula specifying predictors from sample metadata (e.g. `~ treatment + diet`).
#'   If the left-hand side is omitted, the community matrix is automatically supplied.
#' @param assay Assay name to use. Defaults to `"counts"`.
#' @param method Distance metric to pass to `vegan::adonis2` (e.g. `"bray"`, `"jaccard"`). Defaults to `"bray"`.
#' @param permutations Number of permutations. Defaults to `999`.
#' @param by How terms are assessed: `"terms"`, `"margin"`, or `NULL`. Defaults to `"terms"`.
#' @param ... Additional arguments passed to `vegan::adonis2`.
#' @return A tidy tibble with columns:
#'   \itemize{
#'     \item `term`: Predictor name, Residual, or Total.
#'     \item `df`: Degrees of freedom.
#'     \item `sum_sq`: Sum of squares.
#'     \item `r2`: Coefficient of determination (\eqn{R^2}).
#'     \item `f_stat`: Pseudo-F statistic.
#'     \item `p_value`: Permutation p-value.
#'   }
#' @export
run_permanova <- function(tb, formula, assay = "counts", method = "bray", permutations = 999, by = "terms", ...) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for `run_permanova()`. Please install it.", call. = FALSE)
  }

  veg <- to_vegan(tb, assay = assay)
  comm <- veg$comm
  env <- as.data.frame(veg$env)
  rownames(env) <- env$sample_id

  # Ensure left-hand side of formula references the community matrix
  f_char <- deparse(formula)
  if (startsWith(f_char, "~")) {
    formula <- stats::as.formula(paste("comm", f_char))
  }

  ad_res <- vegan::adonis2(
    formula = formula,
    data = env,
    method = method,
    permutations = permutations,
    by = by,
    ...
  )

  # Convert adonis2 table to tidy tibble
  res_df <- as.data.frame(ad_res)
  terms <- rownames(res_df)

  df_out <- tibble::tibble(
    term = terms,
    df = as.numeric(res_df$Df),
    sum_sq = as.numeric(res_df$SumOfSqs),
    r2 = if ("R2" %in% colnames(res_df)) as.numeric(res_df$R2) else NA_real_,
    f_stat = if ("F" %in% colnames(res_df)) as.numeric(res_df$F) else NA_real_,
    p_value = if ("Pr(>F)" %in% colnames(res_df)) as.numeric(res_df$`Pr(>F)`) else NA_real_
  )

  attr(df_out, "method") <- method
  attr(df_out, "permutations") <- permutations
  attr(df_out, "by") <- by

  df_out
}

#' Test for Homogeneity of Multivariate Dispersions via vegan::betadisper
#'
#' Quantifies and tests within-group multivariate dispersion using `vegan::betadisper`
#' and permutation tests, returning tidy summary tibbles.
#'
#' @param tb A `tidy_microbiome` object.
#' @param group Character string naming the grouping variable in sample metadata.
#' @param assay Assay name to use. Defaults to `"counts"`.
#' @param method Distance metric to compute (e.g. `"bray"`, `"jaccard"`). Defaults to `"bray"`.
#' @param permutations Number of permutations for dispersion test. Defaults to `999`.
#' @param ... Additional arguments passed to `vegan::betadisper`.
#' @return An S3 object of class `tidybiome_betadisper` containing:
#'   \itemize{
#'     \item `sample_distances`: Tibble with `sample_id`, `group`, and `distance_to_centroid`.
#'     \item `group_summary`: Tibble with per-group counts, mean distances, and standard errors.
#'     \item `test`: Tibble with F-statistic, degrees of freedom, and permutation p-value.
#'   }
#' @export
run_betadisper <- function(tb, group, assay = "counts", method = "bray", permutations = 999, ...) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for `run_betadisper()`. Please install it.", call. = FALSE)
  }
  if (!group %in% colnames(tb)) {
    stop(sprintf("Group '%s' not found in sample metadata.", group), call. = FALSE)
  }

  veg <- to_vegan(tb, assay = assay)
  comm <- veg$comm
  group_vec <- as.factor(veg$env[[group]])

  dist_mat <- vegan::vegdist(comm, method = method)
  dots <- list(...)
  if (!"add" %in% names(dots)) {
    dots$add <- "lingoes"
  }
  dots$d <- dist_mat
  dots$group <- group_vec

  mod <- do.call(vegan::betadisper, dots)
  perm_test <- vegan::permutest(mod, permutations = permutations)

  # Sample distances tibble
  sample_dist_df <- tibble::tibble(
    sample_id = rownames(comm),
    group = group_vec,
    distance_to_centroid = as.numeric(mod$distances)
  )

  # Group summary
  group_summary_df <- sample_dist_df %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_distance = mean(.data$distance_to_centroid, na.rm = TRUE),
      se = stats::sd(.data$distance_to_centroid, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups = "drop"
    )

  # Test table
  tab <- as.data.frame(perm_test$tab)
  test_df <- tibble::tibble(
    term = rownames(tab)[1],
    df = as.numeric(tab$Df[1]),
    sum_sq = as.numeric(tab$`Sum Sq`[1]),
    mean_sq = as.numeric(tab$`Mean Sq`[1]),
    f_stat = as.numeric(tab$F[1]),
    p_value = as.numeric(tab$`Pr(>F)`[1])
  )

  res <- list(
    sample_distances = sample_dist_df,
    group_summary = group_summary_df,
    test = test_df,
    model = mod
  )
  class(res) <- c("tidybiome_betadisper", "list")
  res
}

#' @export
print.tidybiome_betadisper <- function(x, ...) {
  cat("-- tidybiome_betadisper (Multivariate Dispersion Test) --\n")
  cat(sprintf("  * Statistic F: %.4f (p = %.4f)\n", x$test$f_stat[1], x$test$p_value[1]))
  cat("  * Group centroids:\n")
  for (i in seq_len(nrow(x$group_summary))) {
    cat(sprintf("    - %s: mean distance = %.4f (+/- %.4f)\n",
                x$group_summary$group[i],
                x$group_summary$mean_distance[i],
                x$group_summary$se[i]))
  }
  invisible(x)
}

#' Run Non-Metric Multidimensional Scaling via vegan::metaMDS
#'
#' @param tb A `tidy_microbiome` object.
#' @param assay Assay name to use. Defaults to `"counts"`.
#' @param method Dissimilarity index to pass to `vegan::metaMDS`. Defaults to `"bray"`.
#' @param k Number of dimensions. Defaults to `2`.
#' @param trymax Maximum number of random starts. Defaults to `50`.
#' @param trace Numeric trace output (0 for silent). Defaults to `0`.
#' @param ... Additional arguments passed to `vegan::metaMDS`.
#' @return A tibble of ordination coordinates (`NMDS1`, `NMDS2`, etc.) joined with sample metadata,
#'   with `stress` stored as an attribute.
#' @export
run_nmds <- function(tb, assay = "counts", method = "bray", k = 2, trymax = 50, trace = 0, ...) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for `run_nmds()`. Please install it.", call. = FALSE)
  }

  veg <- to_vegan(tb, assay = assay)
  comm <- veg$comm

  nmds_mod <- vegan::metaMDS(
    comm = comm,
    distance = method,
    k = k,
    trymax = trymax,
    trace = trace,
    ...
  )

  # Extract points
  coords <- as.data.frame(nmds_mod$points)
  colnames(coords) <- paste0("NMDS", seq_len(ncol(coords)))
  coords$sample_id <- rownames(coords)

  # Join with sample metadata
  meta_df <- tibble::as_tibble(tb)
  out <- dplyr::left_join(coords, meta_df, by = "sample_id")
  out <- tibble::as_tibble(out)

  attr(out, "stress") <- nmds_mod$stress
  attr(out, "converged") <- nmds_mod$converged
  attr(out, "method") <- method

  out
}

#' Distance-Based Redundancy Analysis via vegan::dbrda
#'
#' Performs constrained ordination (db-RDA) to model microbial community variation
#' explained by environmental predictors or experimental conditions.
#'
#' @param tb A `tidy_microbiome` object.
#' @param formula A model formula specifying predictors from sample metadata (e.g. `~ treatment + age`).
#' @param assay Assay name to use. Defaults to `"counts"`.
#' @param distance Distance metric to compute (e.g. `"bray"`, `"raitchison"`). Defaults to `"bray"`.
#' @param ... Additional arguments passed to `vegan::dbrda`.
#' @return An S3 object of class `tidybiome_dbrda` containing:
#'   \itemize{
#'     \item `samples`: Tibble of sample coordinates along db-RDA axes joined with metadata.
#'     \item `biplot`: Tibble of constraint vector loadings (continuous predictors).
#'     \item `centroids`: Tibble of factor centroids (categorical predictors).
#'     \item `variance_explained`: Proportion of variance explained by constrained axes.
#'     \item `model`: The underlying `vegan::dbrda` model object.
#'   }
#' @export
run_dbrda <- function(tb, formula, assay = "counts", distance = "bray", ...) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for `run_dbrda()`. Please install it.", call. = FALSE)
  }

  veg <- to_vegan(tb, assay = assay)
  comm <- veg$comm
  env <- as.data.frame(veg$env)
  rownames(env) <- env$sample_id

  # Left hand side of formula references the community matrix
  f_char <- deparse(formula)
  if (startsWith(f_char, "~")) {
    formula <- stats::as.formula(paste("comm", f_char))
  }

  mod <- vegan::dbrda(
    formula = formula,
    data = env,
    distance = distance,
    ...
  )

  # Variance explained
  eig_constrained <- mod$CCA$eig
  total_inertia <- mod$tot.chi
  var_exp <- if (!is.null(eig_constrained) && total_inertia > 0) {
    eig_constrained / total_inertia
  } else {
    numeric(0)
  }
  names(var_exp) <- paste0("dbRDA", seq_along(var_exp))

  # Sample coordinates (sites in CCA space)
  site_scores <- tryCatch(as.data.frame(vegan::scores(mod, display = "sites", choices = seq_along(eig_constrained))), error = function(e) NULL)
  if (is.null(site_scores) || ncol(site_scores) == 0) {
    site_scores <- as.data.frame(mod$CCA$u %*% diag(sqrt(pmax(0, mod$CCA$eig))))
  }
  colnames(site_scores) <- paste0("dbRDA", seq_len(ncol(site_scores)))
  site_scores$sample_id <- rownames(comm)

  samples_df <- dplyr::left_join(site_scores, tibble::as_tibble(tb), by = "sample_id")
  samples_df <- tibble::as_tibble(samples_df)

  # Biplot vectors for continuous terms
  biplot_scores <- tryCatch({
    bp <- as.data.frame(vegan::scores(mod, display = "bp"))
    if (!is.null(bp) && nrow(bp) > 0) {
      bp$term <- rownames(bp)
      colnames(bp)[seq_along(eig_constrained)] <- paste0("dbRDA", seq_along(eig_constrained))
      tibble::as_tibble(bp)
    } else {
      tibble::tibble()
    }
  }, error = function(e) tibble::tibble())

  # Centroids for categorical terms
  centroids_scores <- tryCatch({
    cnt <- as.data.frame(vegan::scores(mod, display = "cn"))
    if (!is.null(cnt) && nrow(cnt) > 0) {
      cnt$term <- rownames(cnt)
      colnames(cnt)[seq_along(eig_constrained)] <- paste0("dbRDA", seq_along(eig_constrained))
      tibble::as_tibble(cnt)
    } else {
      tibble::tibble()
    }
  }, error = function(e) tibble::tibble())

  res <- list(
    samples = samples_df,
    biplot = biplot_scores,
    centroids = centroids_scores,
    variance_explained = var_exp,
    model = mod
  )
  class(res) <- c("tidybiome_dbrda", "list")
  res
}

#' @export
print.tidybiome_dbrda <- function(x, ...) {
  n_ax <- length(x$variance_explained)
  cat("-- tidybiome_dbrda (Distance-Based Redundancy Analysis) --\n")
  cat(sprintf("  * Constrained axes: %d\n", n_ax))
  if (n_ax >= 2) {
    cat(sprintf("  * Inertia explained: dbRDA1 = %.2f%%, dbRDA2 = %.2f%%\n",
                x$variance_explained[1] * 100, x$variance_explained[2] * 100))
  }
  cat(sprintf("  * Samples: %d\n", nrow(x$samples)))
  invisible(x)
}

