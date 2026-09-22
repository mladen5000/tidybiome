#' Calculate Sample-Level Quality Control (QC) Metrics
#'
#' Computes key quality control and diagnostic statistics for each sample in a
#' `tidy_microbiome` object, including total library size, observed feature richness,
#' sample-level sparsity, dominance of the top taxon, and Shannon entropy.
#'
#' @param tb A `tidy_microbiome` object.
#' @param assay Name of the assay to compute metrics on (defaults to `"counts"`).
#' @param augment Logical; if `TRUE` (default), the QC metrics are appended directly as columns
#'   to the sample metadata within the returned `tidy_microbiome` object, enabling seamless
#'   filtering in dplyr pipelines. If `FALSE`, returns a tidy [tibble::tibble] of QC metrics.
#' @return If `augment = TRUE`, an updated `tidy_microbiome` object with added columns:
#'   \describe{
#'     \item{qc_total_reads}{Total sum of read counts in the sample.}
#'     \item{qc_n_features}{Number of features (taxa) observed with count > 0.}
#'     \item{qc_sparsity}{Proportion of unobserved features in the sample (fraction of zeros).}
#'     \item{qc_top_taxon_share}{Proportion of total sample library represented by the single most abundant taxon.}
#'     \item{qc_shannon}{Shannon diversity / entropy index calculated from sample proportions.}
#'   }
#'   If `augment = FALSE`, returns a tibble with `sample_id` and the above columns.
#' @export
#' @examples
#' data(gut_microbiome)
#' # Augment metadata with QC metrics and filter high-quality samples
#' tb_clean <- calc_qc_metrics(gut_microbiome)
#' tb_filtered <- tb_clean[tb_clean$qc_total_reads >= 5000, ]
#'
#' # Or extract QC metrics as a tibble
#' qc_df <- calc_qc_metrics(gut_microbiome, augment = FALSE)
calc_qc_metrics <- function(tb, assay = "counts", augment = TRUE) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  mat <- assay(tb, assay)
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    stop("Specified assay is empty or not found.", call. = FALSE)
  }

  sample_ids <- colnames(mat)
  n_taxa <- nrow(mat)

  total_reads <- colSums(mat)
  n_features <- colSums(mat > 0)
  sparsity <- (n_taxa - n_features) / n_taxa

  top_counts <- apply(mat, 2, max)
  top_share <- ifelse(total_reads > 0, top_counts / total_reads, 0)

  shannon <- apply(mat, 2, function(x) {
    s <- sum(x)
    if (s <= 0) return(0)
    p <- x[x > 0] / s
    -sum(p * log(p))
  })

  qc_tbl <- tibble::tibble(
    sample_id = sample_ids,
    qc_total_reads = total_reads,
    qc_n_features = n_features,
    qc_sparsity = sparsity,
    qc_top_taxon_share = top_share,
    qc_shannon = shannon
  )

  if (augment) {
    current_meta <- tibble::as_tibble(tb)
    qc_col_names <- c("qc_total_reads", "qc_n_features", "qc_sparsity", "qc_top_taxon_share", "qc_shannon")
    current_meta <- current_meta[, setdiff(colnames(current_meta), qc_col_names), drop = FALSE]

    merged_meta <- dplyr::left_join(current_meta, qc_tbl, by = "sample_id")
    reconstruct_tidy_microbiome(merged_meta, tb, sync_samples = FALSE)
  } else {
    qc_tbl
  }
}
