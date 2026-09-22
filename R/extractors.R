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
#' @param taxa Optional character vector of taxon IDs to filter down before long-format expansion.
#' @param samples Optional character vector of sample IDs to filter down before long-format expansion.
#' @param include_metadata Logical; whether to join sample metadata in long format (default: `TRUE`).
#' @param include_taxonomy Logical; whether to join taxonomy metadata in long format (default: `TRUE`).
#'
#' @return A [tibble::tbl_df] if `long = TRUE`, or a numeric matrix if `long = FALSE`.
#' @export
tidy_abundance <- function(tb,
                           assay = "counts",
                           rank = NULL,
                           long = TRUE,
                           taxa = NULL,
                           samples = NULL,
                           include_metadata = TRUE,
                           include_taxonomy = TRUE) {
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

  # Optional pre-filtering before expansion
  if (!is.null(taxa)) {
    common_taxa <- intersect(taxa, rownames(mat))
    mat <- mat[common_taxa, , drop = FALSE]
  }
  if (!is.null(samples)) {
    common_samples <- intersect(samples, colnames(mat))
    mat <- mat[, common_samples, drop = FALSE]
  }

  if (!long) {
    return(mat)
  }

  sample_names <- colnames(mat)
  taxon_names  <- rownames(mat)
  n_taxa <- length(taxon_names)
  n_samp <- length(sample_names)

  df_long <- tibble::tibble(
    taxon_id  = rep(taxon_names, times = n_samp),
    sample_id = rep(sample_names, each = n_taxa),
    abundance = as.vector(mat)
  )

  # Join with sample metadata if requested
  if (include_metadata) {
    sample_df <- tibble::as_tibble(tb)
    if (!is.null(samples)) {
      sample_df <- sample_df[sample_df$sample_id %in% sample_names, , drop = FALSE]
    }
    df_long <- dplyr::left_join(df_long, sample_df, by = "sample_id")
  }

  # Join with taxonomy if requested
  if (include_taxonomy) {
    tax_df <- attr(tb, "tax_table")
    if (!is.null(tax_df) && nrow(tax_df) > 0) {
      if (!is.null(taxa)) {
        tax_df <- tax_df[tax_df$taxon_id %in% taxon_names, , drop = FALSE]
      }
      df_long <- dplyr::left_join(df_long, tax_df, by = "taxon_id")
    }
  }

  df_long
}

#' Aggregate Taxa to a Higher Taxonomic Rank
#'
#' @description
#' Agglomerates/merges features (ASVs/OTUs) by a specified taxonomic rank (e.g., Phylum, Family, Genus)
#' by summing abundances across all assays using compiled C-level primitives.
#'
#' @param tb A `tidy_microbiome` object.
#' @param rank Character string specifying the taxonomic column to aggregate by.
#' @param na.rm Logical. If `TRUE`, removes unassigned taxa. If `FALSE` (default), groups them as `"Unclassified"`.
#'
#' @return An updated `tidy_microbiome` object aggregated at the specified rank.
#' @export
#' @examples
#' data(gut_microbiome)
#' tb_phylum <- aggregate_taxa(gut_microbiome, rank = "Phylum")
#' dim(assay(tb_phylum))
aggregate_taxa <- function(tb, rank, na.rm = FALSE) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  tax_df <- attr(tb, "tax_table")
  if (is.null(tax_df)) {
    stop("Taxonomy table is missing; cannot aggregate taxa.", call. = FALSE)
  }

  if (!rank %in% colnames(tax_df)) {
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

  # Vectorized C-level rowsum across all assays
  new_assays <- lapply(assays, function(mat) {
    if (na.rm) mat <- mat[keep_idx, , drop = FALSE]
    base::rowsum(mat, group = group_vec, reorder = FALSE, na.rm = TRUE)
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

#' Extract Assay Matrix from tidy_microbiome
#'
#' @param x A `tidy_microbiome` object.
#' @param name Character string naming the assay to extract (e.g. `"counts"`, `"rclr"`, `"tss"`).
#'   Defaults to `"counts"`.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric matrix with taxa as rows and samples as columns.
#' @export
assay <- function(x, name = "counts", ...) {
  UseMethod("assay")
}

#' @export
assay.tidy_microbiome <- function(x, name = "counts", ...) {
  assays <- attr(x, "assays")
  if (is.null(assays) || !name %in% names(assays)) {
    stop(sprintf("Assay '%s' not found in tidy_microbiome. Available: %s",
                 name, paste(names(assays), collapse = ", ")), call. = FALSE)
  }
  assays[[name]]
}

#' Extract Taxonomy Table from tidy_microbiome
#'
#' @param x A `tidy_microbiome` object.
#' @param ... Additional arguments (not used).
#'
#' @return A [tibble::tbl_df] with taxonomic lineage annotations.
#' @export
tax_table <- function(x, ...) {
  UseMethod("tax_table")
}

#' @export
tax_table.tidy_microbiome <- function(x, ...) {
  tidy_taxa(x)
}
