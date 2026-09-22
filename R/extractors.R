#' Extract Taxonomy Annotations as a Tidy Tibble
#'
#' @param tb A `tidy_microbiome` object.
#' @return A [tibble::tbl_df] with `taxon_id` and taxonomic ranks.
#' @export
tidy_taxa <- function(tb) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  attr(tb, "tax_table")
}

#' Extract Abundance Data in Tidy Format
#'
#' @param tb A `tidy_microbiome` object.
#' @param assay Character string naming the assay to extract (e.g. `"counts"`, `"relabundance"`, `"rclr"`).
#'   Defaults to `"counts"`.
#' @param rank Optional taxonomic rank to aggregate by before extraction (e.g. `"Genus"`, `"Phylum"`).
#' @param long Logical. If `TRUE` (default), returns a fully denormalized long tibble
#'   containing `sample_id`, `taxon_id`, `abundance`, sample metadata, and taxonomy.
#'   If `FALSE`, returns the abundance matrix directly.
#'
#' @return A [tibble::tbl_df] if `long = TRUE`, or a numeric matrix if `long = FALSE`.
#' @export
tidy_abundance <- function(tb, assay = "counts", rank = NULL, long = TRUE) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  if (!is.null(rank)) {
    tb <- aggregate_taxa(tb, rank = rank)
  }

  assays <- attr(tb, "assays")
  if (!assay %in% names(assays)) {
    stop(sprintf("Assay '%s' not found in tidy_microbiome. Available: %s",
                 assay, paste(names(assays), collapse = ", ")), call. = FALSE)
  }

  mat <- assays[[assay]]

  if (!long) {
    return(mat)
  }

  # Build long format efficiently without expand.grid memory overhead
  sample_names <- colnames(mat)
  taxon_names  <- rownames(mat)
  n_taxa <- length(taxon_names)
  n_samp <- length(sample_names)

  df_long <- tibble::tibble(
    taxon_id  = rep(taxon_names, times = n_samp),
    sample_id = rep(sample_names, each = n_taxa),
    abundance = as.vector(mat)
  )

  # Join with sample metadata
  sample_df <- tibble::as_tibble(tb)
  res <- dplyr::left_join(df_long, sample_df, by = "sample_id")

  # Join with taxonomy
  tax_df <- attr(tb, "tax_table")
  if (!is.null(tax_df) && nrow(tax_df) > 0) {
    res <- dplyr::left_join(res, tax_df, by = "taxon_id")
  }

  res
}

#' Aggregate Taxa to a Higher Taxonomic Rank
#'
#' @description
#' Agglomerates/merges features (ASVs/OTUs) by a specified taxonomic rank (e.g., Phylum, Family, Genus)
#' by summing abundances across all assays.
#'
#' @param tb A `tidy_microbiome` object.
#' @param rank Character string specifying the taxonomic column to aggregate by.
#' @param na.rm Logical. If `TRUE`, removes unassigned taxa. If `FALSE` (default), groups them as `"Unclassified"`.
#'
#' @return A new `tidy_microbiome` object aggregated at the specified rank.
#' @export
aggregate_taxa <- function(tb, rank, na.rm = FALSE) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  tax_df <- attr(tb, "tax_table")
  if (is.null(tax_df) || !rank %in% colnames(tax_df)) {
    stop(sprintf("Rank '%s' not found in taxonomy table.", rank), call. = FALSE)
  }

  group_vec <- as.character(tax_df[[rank]])
  group_vec[is.na(group_vec) | group_vec == ""] <- paste0("Unclassified_", rank)

  if (na.rm) {
    keep_idx <- !grepl("^Unclassified_", group_vec)
    group_vec <- group_vec[keep_idx]
    tax_df <- tax_df[keep_idx, , drop = FALSE]
  }

  unique_ranks <- unique(group_vec)
  assays <- attr(tb, "assays")

  new_assays <- lapply(assays, function(mat) {
    if (na.rm) mat <- mat[keep_idx, , drop = FALSE]
    res_mat <- matrix(0, nrow = length(unique_ranks), ncol = ncol(mat),
                      dimnames = list(unique_ranks, colnames(mat)))
    for (ur in unique_ranks) {
      idx <- which(group_vec == ur)
      if (length(idx) == 1) {
        res_mat[ur, ] <- mat[idx, ]
      } else {
        res_mat[ur, ] <- colSums(mat[idx, , drop = FALSE], na.rm = TRUE)
      }
    }
    res_mat
  })

  # Aggregate taxonomy table
  rank_cols <- colnames(tax_df)
  target_idx <- which(rank_cols == rank)
  kept_cols <- rank_cols[1:target_idx]

  new_tax <- data.frame(taxon_id = unique_ranks, stringsAsFactors = FALSE)
  for (col in kept_cols) {
    if (col != "taxon_id") {
      val_map <- tapply(tax_df[[col]], group_vec, function(v) {
        uv <- unique(v[!is.na(v) & v != ""])
        if (length(uv) == 1) uv else NA_character_
      })
      new_tax[[col]] <- as.character(val_map[unique_ranks])
    }
  }

  orig_tree <- attr(tb, "phy_tree")
  new_tree <- NULL
  if (!is.null(orig_tree) && requireNamespace("ape", quietly = TRUE)) {
    rep_taxa <- character(length(unique_ranks))
    names(rep_taxa) <- unique_ranks
    for (ur in unique_ranks) {
      members <- tax_df$taxon_id[group_vec == ur]
      in_tree <- intersect(members, orig_tree$tip.label)
      if (length(in_tree) > 0) {
        rep_taxa[ur] <- in_tree[1]
      }
    }
    valid_reps <- rep_taxa[rep_taxa != ""]
    if (length(valid_reps) >= 2) {
      pruned <- ape::keep.tip(orig_tree, valid_reps)
      name_map <- stats::setNames(names(valid_reps), valid_reps)
      pruned$tip.label <- as.character(name_map[pruned$tip.label])
      new_tree <- pruned
    }
  }

  sample_df <- tibble::as_tibble(tb)
  out <- sample_df
  attr(out, "assays")    <- new_assays
  attr(out, "tax_table") <- tibble::as_tibble(new_tax)
  attr(out, "phy_tree")  <- new_tree
  attr(out, "metadata")  <- attr(tb, "metadata")
  class(out) <- c("tidy_microbiome", "tbl_df", "tbl", "data.frame")
  out
}
