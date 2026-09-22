#' Calculate Cross-Associations Between Taxa and Metadata Variables
#'
#' Computes pairwise associations (Spearman or Pearson correlations and significance)
#' between microbial feature abundances and continuous metadata covariates (e.g. host age, BMI, clinical metrics).
#'
#' @param tb A `tidy_microbiome` object.
#' @param variables Character vector of numeric metadata columns in `tb` to correlate with taxa.
#'   If `NULL`, automatically selects all numeric columns in sample metadata. Defaults to `NULL`.
#' @param assay Name of assay to evaluate. Defaults to `"counts"`.
#' @param method Correlation method: `"spearman"` or `"pearson"`. Defaults to `"spearman"`.
#' @param p_adj_method Multiple testing correction method passed to `stats::p.adjust` (e.g. `"BH"`, `"bonferroni"`).
#'   Defaults to `"BH"`.
#' @param sort Logical; if `TRUE`, sorts output by adjusted p-value (`padj`). Defaults to `TRUE`.
#' @return A tidy tibble with columns:
#'   \itemize{
#'     \item `taxon_id`: Taxon identifier.
#'     \item `variable`: Metadata covariate name.
#'     \item `correlation`: Estimated correlation coefficient.
#'     \item `p_value`: Raw association p-value from `cor.test`.
#'     \item `padj`: Multiple-testing adjusted p-value.
#'   }
#'   Joined with taxonomic lineages from `tax_table`.
#' @export
calc_cross_association <- function(tb,
                                    variables = NULL,
                                    assay = "counts",
                                    method = c("spearman", "pearson"),
                                    p_adj_method = "BH",
                                    sort = TRUE) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  method <- match.arg(method)

  assays <- attr(tb, "assays")
  if (!assay %in% names(assays)) {
    stop(sprintf("Assay '%s' not found in `tb`.", assay), call. = FALSE)
  }

  mat <- assays[[assay]]
  meta_df <- as.data.frame(tibble::as_tibble(tb))

  # Select numeric variables
  if (is.null(variables)) {
    is_num <- vapply(meta_df, is.numeric, logical(1))
    variables <- names(meta_df)[is_num]
    variables <- setdiff(variables, c("sample_id"))
  } else {
    missing_vars <- setdiff(variables, colnames(meta_df))
    if (length(missing_vars) > 0) {
      stop(sprintf("Variables not found in sample metadata: %s", paste(missing_vars, collapse = ", ")), call. = FALSE)
    }
    non_num <- variables[!vapply(meta_df[variables], is.numeric, logical(1))]
    if (length(non_num) > 0) {
      stop(sprintf("Selected variables must be numeric: %s", paste(non_num, collapse = ", ")), call. = FALSE)
    }
  }

  if (length(variables) == 0) {
    stop("No continuous/numeric metadata variables available for cross-association.", call. = FALSE)
  }

  n_taxa <- nrow(mat)
  taxa_names <- rownames(mat)
  results_list <- vector("list", length(variables) * n_taxa)
  idx <- 1

  for (v in variables) {
    var_vals <- meta_df[[v]]
    for (i in seq_len(n_taxa)) {
      tax_vals <- mat[i, ]

      # Safe cor.test
      ct <- tryCatch({
        stats::cor.test(tax_vals, var_vals, method = method, exact = FALSE)
      }, error = function(e) {
        list(estimate = NA_real_, p.value = NA_real_)
      })

      results_list[[idx]] <- tibble::tibble(
        taxon_id = taxa_names[i],
        variable = v,
        correlation = as.numeric(ct$estimate),
        p_value = as.numeric(ct$p.value)
      )
      idx <- idx + 1
    }
  }

  res <- dplyr::bind_rows(results_list)

  # Adjust p-values per variable or cohort-wide
  res$padj <- stats::p.adjust(res$p_value, method = p_adj_method)

  # Attach tax_table
  tax_df <- attr(tb, "tax_table")
  if (!is.null(tax_df)) {
    res <- dplyr::left_join(res, tax_df, by = "taxon_id")
  }

  if (sort) {
    res <- dplyr::arrange(res, .data$padj, .data$p_value)
  }

  res
}

#' Calculate Microbial Co-Occurrence Network
#'
#' Evaluates pairwise co-occurrence correlations across taxa and constructs
#' network nodes and edges filtered by correlation strength and statistical significance.
#'
#' @param tb A `tidy_microbiome` object.
#' @param assay Assay to evaluate. Defaults to `"counts"`.
#' @param method Correlation method: `"spearman"` or `"pearson"`. Defaults to `"spearman"`.
#' @param min_prevalence Minimum proportion of samples where taxon must be detected (default: 0.20).
#' @param r_cutoff Minimum absolute correlation coefficient threshold (default: 0.40).
#' @param p_cutoff Maximum multiple testing adjusted p-value threshold (default: 0.05).
#' @param p_adj_method Multiple testing adjustment method (default: `"BH"`).
#'
#' @return An S3 object of class `tidybiome_network` containing:
#'   \itemize{
#'     \item `nodes`: Tibble of network nodes with `taxon_id`, `degree`, `mean_abundance`, and taxonomy.
#'     \item `edges`: Tibble of network edges with `from`, `to`, `correlation`, `p_value`, `padj`, `weight`, and `direction`.
#'     \item `r_cutoff`: Filtering cutoff for r.
#'     \item `p_cutoff`: Filtering cutoff for p.
#'   }
#' @export
calc_network <- function(tb,
                         assay = "counts",
                         method = c("spearman", "pearson"),
                         min_prevalence = 0.20,
                         r_cutoff = 0.40,
                         p_cutoff = 0.05,
                         p_adj_method = "BH") {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  method <- match.arg(method)

  tb_prev <- filter_prevalent(tb, min_prevalence = min_prevalence, assay = assay)
  assays <- attr(tb_prev, "assays")
  mat <- assays[[assay]]
  n_taxa <- nrow(mat)
  taxa_names <- rownames(mat)

  if (n_taxa < 2) {
    stop("Fewer than 2 taxa passed the prevalence threshold for network construction.", call. = FALSE)
  }

  rel_mat <- calc_relabundance_matrix(mat)
  cor_mat <- stats::cor(t(rel_mat), method = method)

  pairs <- which(upper.tri(cor_mat), arr.ind = TRUE)
  n_pairs <- nrow(pairs)

  from_vec <- taxa_names[pairs[, 1]]
  to_vec   <- taxa_names[pairs[, 2]]
  r_vec    <- cor_mat[upper.tri(cor_mat)]

  p_vec <- numeric(n_pairs)
  n_samp <- ncol(mat)
  for (idx in seq_len(n_pairs)) {
    r_val <- r_vec[idx]
    if (is.na(r_val) || abs(r_val) >= 1) {
      p_vec[idx] <- 0
    } else {
      df <- n_samp - 2
      t_stat <- r_val * sqrt(df / max(1e-15, (1 - r_val^2)))
      p_vec[idx] <- 2 * stats::pt(-abs(t_stat), df = df)
    }
  }

  padj_vec <- stats::p.adjust(p_vec, method = p_adj_method)

  edges_df <- tibble::tibble(
    from = from_vec,
    to = to_vec,
    correlation = as.numeric(r_vec),
    p_value = as.numeric(p_vec),
    padj = as.numeric(padj_vec),
    weight = abs(as.numeric(r_vec)),
    direction = ifelse(r_vec > 0, "positive", "negative")
  )

  sig_edges <- edges_df %>%
    dplyr::filter(!is.na(.data$correlation) & .data$weight >= r_cutoff & .data$padj <= p_cutoff)

  degree_from <- table(sig_edges$from)
  degree_to   <- table(sig_edges$to)
  all_degrees <- stats::setNames(integer(n_taxa), taxa_names)

  if (length(degree_from) > 0) {
    all_degrees[names(degree_from)] <- all_degrees[names(degree_from)] + as.integer(degree_from)
  }
  if (length(degree_to) > 0) {
    all_degrees[names(degree_to)] <- all_degrees[names(degree_to)] + as.integer(degree_to)
  }

  nodes_df <- tibble::tibble(
    taxon_id = taxa_names,
    degree = as.integer(all_degrees[taxa_names]),
    mean_abundance = as.numeric(rowMeans(rel_mat))
  )

  tax_df <- attr(tb_prev, "tax_table")
  if (!is.null(tax_df)) {
    nodes_df <- dplyr::left_join(nodes_df, tax_df, by = "taxon_id")
  }

  res <- list(
    nodes = nodes_df,
    edges = sig_edges,
    all_edges = edges_df,
    r_cutoff = r_cutoff,
    p_cutoff = p_cutoff
  )
  class(res) <- c("tidybiome_network", "list")
  res
}

#' @export
print.tidybiome_network <- function(x, ...) {
  cat("-- tidybiome_network (Microbial Co-Occurrence Network) --\n")
  cat(sprintf("  * Nodes (Taxa): %d\n", nrow(x$nodes)))
  cat(sprintf("  * Significant edges: %d (|r| >= %.2f, padj <= %.2f)\n",
              nrow(x$edges), x$r_cutoff, x$p_cutoff))
  pos_n <- sum(x$edges$direction == "positive")
  neg_n <- sum(x$edges$direction == "negative")
  cat(sprintf("  * Positive: %d, Negative: %d\n", pos_n, neg_n))
  invisible(x)
}

