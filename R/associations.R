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
