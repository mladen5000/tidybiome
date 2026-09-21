#' Calculate Prevalence and Abundance Metrics Across Taxa
#'
#' Evaluates the proportion and number of samples where each taxon exceeds a specified
#' detection threshold, along with mean and median abundances.
#'
#' @param tb A `tidy_microbiome` object.
#' @param assay Assay to evaluate. Defaults to `"counts"`.
#' @param detection Abundance threshold above which a taxon is considered detected. Defaults to `0`.
#' @param sort Logical; if `TRUE`, sorts taxa in descending order of prevalence. Defaults to `TRUE`.
#' @return A tibble with columns:
#'   \itemize{
#'     \item `taxon_id`: Unique identifier for each taxon.
#'     \item `prevalence`: Proportion of samples where taxon was detected (\eqn{[0, 1]}).
#'     \item `prevalence_n`: Count of samples where taxon was detected.
#'     \item `mean_abundance`: Average abundance across all samples.
#'     \item `median_abundance`: Median abundance across all samples.
#'   }
#'   Joined with available taxonomic ranks from `tax_table`.
#' @export
calc_prevalence <- function(tb, assay = "counts", detection = 0, sort = TRUE) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  assays <- attr(tb, "assays")
  if (!assay %in% names(assays)) {
    stop(sprintf("Assay '%s' not found in `tb`.", assay), call. = FALSE)
  }

  mat <- assays[[assay]]
  n_samples <- ncol(mat)

  detected_counts <- rowSums(mat > detection)
  prev_prop <- detected_counts / n_samples
  mean_abund <- rowMeans(mat)
  med_abund <- apply(mat, 1, stats::median)

  res <- tibble::tibble(
    taxon_id = rownames(mat),
    prevalence = as.numeric(prev_prop),
    prevalence_n = as.integer(detected_counts),
    mean_abundance = as.numeric(mean_abund),
    median_abundance = as.numeric(med_abund)
  )

  tax_df <- attr(tb, "tax_table")
  if (!is.null(tax_df)) {
    res <- dplyr::left_join(res, tax_df, by = "taxon_id")
  }

  if (sort) {
    res <- dplyr::arrange(res, dplyr::desc(.data$prevalence), dplyr::desc(.data$mean_abundance))
  }

  res
}

#' Identify Dominant Taxa per Sample and Cohort-Wide
#'
#' Finds the taxon with highest relative abundance in each sample and computes cohort-wide dominance.
#'
#' @param tb A `tidy_microbiome` object.
#' @param assay Assay to evaluate. Defaults to `"counts"`.
#' @param rank Optional taxonomic rank to agglomerate before determining dominance. Defaults to `NULL`.
#' @param add_to_metadata Logical; if `TRUE`, adds column `dominant_taxa` directly to `tb`'s sample metadata.
#'   If `FALSE`, returns a summary tibble of dominant taxa across the cohort. Defaults to `TRUE`.
#' @return A `tidy_microbiome` object with `dominant_taxa` in metadata (if `add_to_metadata = TRUE`),
#'   or a tibble summarizing dominance frequencies across samples.
#' @export
calc_dominant <- function(tb, assay = "counts", rank = NULL, add_to_metadata = TRUE) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  target_tb <- tb
  if (!is.null(rank)) {
    target_tb <- aggregate_taxa(tb, rank = rank)
  }

  assays <- attr(target_tb, "assays")
  if (!assay %in% names(assays)) {
    stop(sprintf("Assay '%s' not found in `tb`.", assay), call. = FALSE)
  }

  mat <- assays[[assay]]
  # Find row index with maximum value for each column (sample)
  max_indices <- apply(mat, 2, which.max)
  dominant_vec <- rownames(mat)[max_indices]

  if (add_to_metadata) {
    tb$dominant_taxa <- dominant_vec
    return(tb)
  }

  # Return summary table
  dom_summary <- tibble::tibble(dominant_taxa = dominant_vec) %>%
    dplyr::group_by(.data$dominant_taxa) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      frequency = dplyr::n() / ncol(mat),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$n_samples))

  dom_summary
}

#' Filter Microbiome by Prevalence and Abundance Thresholds
#'
#' Keeps only taxa exceeding the specified minimum prevalence (fraction of samples)
#' and mean relative abundance, keeping count matrices and taxonomy in sync.
#'
#' @param tb A `tidy_microbiome` object.
#' @param min_prevalence Minimum proportion of samples where taxon must be detected (0 to 1). Defaults to `0.10`.
#' @param min_abundance Minimum mean abundance across samples. Defaults to `0`.
#' @param assay Assay to evaluate. Defaults to `"counts"`.
#' @param detection Abundance threshold for presence. Defaults to `0`.
#' @return A filtered `tidy_microbiome` object.
#' @export
filter_prevalent <- function(tb, min_prevalence = 0.10, min_abundance = 0, assay = "counts", detection = 0) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  prev_df <- calc_prevalence(tb, assay = assay, detection = detection, sort = FALSE)

  passing_taxa <- prev_df %>%
    dplyr::filter(.data$prevalence >= min_prevalence, .data$mean_abundance >= min_abundance) %>%
    dplyr::pull(.data$taxon_id)

  if (length(passing_taxa) == 0) {
    stop("No taxa passed the specified prevalence and abundance thresholds.", call. = FALSE)
  }

  # Subset all assays
  assays <- attr(tb, "assays")
  for (nm in names(assays)) {
    assays[[nm]] <- assays[[nm]][passing_taxa, , drop = FALSE]
  }
  attr(tb, "assays") <- assays

  # Subset tax_table
  tax_df <- attr(tb, "tax_table")
  if (!is.null(tax_df)) {
    attr(tb, "tax_table") <- tax_df[tax_df$taxon_id %in% passing_taxa, , drop = FALSE]
  }

  # Prune tree if attached
  tree <- attr(tb, "phy_tree")
  if (!is.null(tree) && requireNamespace("ape", quietly = TRUE)) {
    common_tips <- intersect(tree$tip.label, passing_taxa)
    if (length(common_tips) > 1) {
      attr(tb, "phy_tree") <- ape::keep.tip(tree, common_tips)
    } else {
      attr(tb, "phy_tree") <- NULL
    }
  }

  tb
}
