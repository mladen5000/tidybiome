#' Create a tidy_microbiome Object
#'
#' @description
#' Creates a `tidy_microbiome` S3 object that inherits from [tibble::tbl_df].
#' The object behaves as a sample metadata tibble on the outside, while keeping
#' count matrices (assays), taxonomy annotations, phylogenetic trees, and
#' experiment metadata synchronized under the hood.
#'
#' @param counts A numeric matrix of counts (taxa in rows, samples in columns).
#' @param sample_data Optional `data.frame` or `tbl_df` of sample metadata.
#'   Must include a `sample_id` column or rownames matching `colnames(counts)`.
#' @param tax_table Optional `data.frame` or `tbl_df` of taxonomy annotations.
#'   Must include a `taxon_id` column or rownames matching `rownames(counts)`.
#' @param phy_tree Optional phylogenetic tree (e.g. of class `phylo`).
#' @param metadata Optional list of arbitrary experiment metadata.
#'
#' @return An object of class `tidy_microbiome`, inheriting from `tbl_df`, `tbl`, and `data.frame`.
#' @export
#'
#' @examples
#' counts <- matrix(c(10, 0, 5, 20), nrow = 2,
#'                  dimnames = list(c("ASV1", "ASV2"), c("S1", "S2")))
#' tb <- tidy_microbiome(counts)
#' tb
tidy_microbiome <- function(counts,
                            sample_data = NULL,
                            tax_table = NULL,
                            phy_tree = NULL,
                            metadata = list()) {
  # 1. Validate counts
  if (!is.matrix(counts)) {
    counts <- as.matrix(counts)
  }
  if (!is.numeric(counts)) {
    stop("`counts` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    stop("`counts` matrix must have at least 1 taxon (row) and 1 sample (column).", call. = FALSE)
  }
  if (any(is.na(counts))) {
    stop("`counts` matrix contains NA or NaN values. Missing values should be imputed or replaced with 0.", call. = FALSE)
  }
  if (any(is.infinite(counts))) {
    stop("`counts` matrix cannot contain Infinite values.", call. = FALSE)
  }
  if (any(counts < 0, na.rm = TRUE)) {
    stop("`counts` matrix cannot contain negative values.", call. = FALSE)
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("`counts` matrix must have both rownames (taxa) and colnames (samples).", call. = FALSE)
  }

  colnames(counts) <- trimws(colnames(counts))
  rownames(counts) <- trimws(rownames(counts))

  sample_names <- colnames(counts)
  taxon_names  <- rownames(counts)

  if (anyDuplicated(sample_names)) {
    stop("Duplicate column names (sample IDs) found in `counts` matrix.", call. = FALSE)
  }
  if (anyDuplicated(taxon_names)) {
    stop("Duplicate row names (taxon IDs) found in `counts` matrix.", call. = FALSE)
  }

  # 2. Process sample_data
  if (is.null(sample_data)) {
    sample_df <- tibble::tibble(
      sample_id = sample_names,
      depth     = colSums(counts, na.rm = TRUE)
    )
  } else {
    sample_df <- as.data.frame(sample_data)
    if (!"sample_id" %in% colnames(sample_df)) {
      if (!is.null(rownames(sample_df)) && all(rownames(sample_df) %in% sample_names)) {
        sample_df$sample_id <- rownames(sample_df)
      } else if (nrow(sample_df) == length(sample_names)) {
        sample_df$sample_id <- sample_names
      } else {
        stop("`sample_data` must contain a `sample_id` column matching `colnames(counts)`.", call. = FALSE)
      }
    }
    sample_df$sample_id <- as.character(sample_df$sample_id)
    if (!all(sample_names %in% sample_df$sample_id)) {
      stop("Not all sample names in `counts` are present in `sample_data$sample_id`.", call. = FALSE)
    }
    # Align order
    sample_df <- sample_df[match(sample_names, sample_df$sample_id), , drop = FALSE]
    if (!"depth" %in% colnames(sample_df)) {
      sample_df$depth <- colSums(counts, na.rm = TRUE)
    }
    sample_df <- tibble::as_tibble(sample_df)
  }

  # 3. Process tax_table
  if (is.null(tax_table)) {
    tax_df <- tibble::tibble(taxon_id = taxon_names)
  } else {
    tax_df <- as.data.frame(tax_table)
    if (!"taxon_id" %in% colnames(tax_df)) {
      if (!is.null(rownames(tax_df)) && all(rownames(tax_df) %in% taxon_names)) {
        tax_df$taxon_id <- rownames(tax_df)
      } else if (nrow(tax_df) == length(taxon_names)) {
        tax_df$taxon_id <- taxon_names
      } else {
        stop("`tax_table` must contain a `taxon_id` column matching `rownames(counts)`.", call. = FALSE)
      }
    }
    tax_df$taxon_id <- as.character(tax_df$taxon_id)
    if (!all(taxon_names %in% tax_df$taxon_id)) {
      stop("Not all taxon names in `counts` are present in `tax_table$taxon_id`.", call. = FALSE)
    }
    tax_df <- tax_df[match(taxon_names, tax_df$taxon_id), , drop = FALSE]
    tax_df <- tibble::as_tibble(tax_df)
  }

  # 4. Construct container
  if (!is.list(metadata)) {
    metadata <- list()
  }
  if (!"distances" %in% names(metadata)) {
    metadata$distances <- list()
  }
  if (!"ordinations" %in% names(metadata)) {
    metadata$ordinations <- list()
  }
  if (!"da_results" %in% names(metadata)) {
    metadata$da_results <- list()
  }

  out <- sample_df
  attr(out, "assays")    <- list(counts = counts)
  attr(out, "tax_table") <- tax_df
  attr(out, "phy_tree")  <- phy_tree
  attr(out, "metadata")  <- metadata

  class(out) <- c("tidy_microbiome", "tbl_df", "tbl", "data.frame")
  out
}

# Helper to restore attributes after tibble transformation
reconstruct_tidy_microbiome <- function(new_df, old_obj, sync_samples = FALSE) {
  assays    <- attr(old_obj, "assays")
  tax_table <- attr(old_obj, "tax_table")
  phy_tree  <- attr(old_obj, "phy_tree")
  metadata  <- attr(old_obj, "metadata")

  if (sync_samples) {
    if (!"sample_id" %in% colnames(new_df)) {
      stop("Cannot subset `tidy_microbiome`: `sample_id` column was dropped.", call. = FALSE)
    }
    new_samples <- as.character(new_df$sample_id)
    # Sync all assay columns
    assays <- lapply(assays, function(m) {
      m[, new_samples, drop = FALSE]
    })
  }

  out <- new_df
  attr(out, "assays")    <- assays
  attr(out, "tax_table") <- tax_table
  attr(out, "phy_tree")  <- phy_tree
  attr(out, "metadata")  <- metadata
  class(out) <- c("tidy_microbiome", "tbl_df", "tbl", "data.frame")
  out
}

#' @importFrom dplyr filter slice arrange mutate select
#' @rawNamespace S3method(dplyr::filter, tidy_microbiome)
#' @rawNamespace S3method(dplyr::slice, tidy_microbiome)
#' @rawNamespace S3method(dplyr::arrange, tidy_microbiome)
#' @rawNamespace S3method(dplyr::mutate, tidy_microbiome)
#' @rawNamespace S3method(dplyr::select, tidy_microbiome)
NULL

#' @method filter tidy_microbiome
#' @export
filter.tidy_microbiome <- function(.data, ...) {
  df_filtered <- NextMethod()
  reconstruct_tidy_microbiome(df_filtered, .data, sync_samples = TRUE)
}

#' @method slice tidy_microbiome
#' @export
slice.tidy_microbiome <- function(.data, ...) {
  df_sliced <- NextMethod()
  reconstruct_tidy_microbiome(df_sliced, .data, sync_samples = TRUE)
}

#' @method arrange tidy_microbiome
#' @export
arrange.tidy_microbiome <- function(.data, ...) {
  df_arranged <- NextMethod()
  reconstruct_tidy_microbiome(df_arranged, .data, sync_samples = TRUE)
}

#' @method mutate tidy_microbiome
#' @export
mutate.tidy_microbiome <- function(.data, ...) {
  df_mutated <- NextMethod()
  reconstruct_tidy_microbiome(df_mutated, .data, sync_samples = FALSE)
}

#' @method select tidy_microbiome
#' @export
select.tidy_microbiome <- function(.data, ...) {
  df_selected <- NextMethod()
  if (!"sample_id" %in% colnames(df_selected)) {
    df_selected$sample_id <- .data$sample_id
    df_selected <- df_selected[, c("sample_id", setdiff(colnames(df_selected), "sample_id")), drop = FALSE]
  }
  reconstruct_tidy_microbiome(df_selected, .data, sync_samples = FALSE)
}

#' @export
`[.tidy_microbiome` <- function(x, i, j, drop = FALSE) {
  is_single_arg <- nargs() <= 2
  res <- NextMethod()
  
  if (!is.data.frame(res)) {
    return(res)
  }
  
  if (is_single_arg) {
    if ("sample_id" %in% colnames(res)) {
      return(reconstruct_tidy_microbiome(tibble::as_tibble(res), x, sync_samples = FALSE))
    } else {
      return(res)
    }
  }
  
  if (!missing(j)) {
    if (!"sample_id" %in% colnames(res)) {
      orig_sample_id <- x$sample_id
      if (!missing(i)) orig_sample_id <- orig_sample_id[i]
      res$sample_id <- orig_sample_id
      res <- res[, c("sample_id", setdiff(colnames(res), "sample_id")), drop = FALSE]
    }
  }
  
  sync_samples <- !missing(i)
  reconstruct_tidy_microbiome(tibble::as_tibble(res), x, sync_samples = sync_samples)
}

#' @export
`[<-.tidy_microbiome` <- function(x, i, j, value) {
  res <- NextMethod()
  reconstruct_tidy_microbiome(tibble::as_tibble(res), x, sync_samples = FALSE)
}

#' Filter Taxa Based on Taxonomy Table Attributes
#'
#' @param tb A `tidy_microbiome` object.
#' @param ... Logical predicates passed to [dplyr::filter] on the taxonomy table.
#'
#' @return A filtered `tidy_microbiome` object with synchronized assays, taxonomy, and tree.
#' @export
filter_taxa <- function(tb, ...) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  
  tax_df <- attr(tb, "tax_table")
  if (is.null(tax_df) || nrow(tax_df) == 0) {
    stop("No taxonomy table found in `tb`.", call. = FALSE)
  }
  
  filtered_tax <- dplyr::filter(tax_df, ...)
  keep_taxa <- as.character(filtered_tax$taxon_id)
  
  if (length(keep_taxa) == 0) {
    stop("No taxa match the specified filter conditions.", call. = FALSE)
  }
  
  assays <- attr(tb, "assays")
  new_assays <- lapply(assays, function(mat) {
    mat[keep_taxa, , drop = FALSE]
  })
  
  tree <- attr(tb, "phy_tree")
  new_tree <- NULL
  if (!is.null(tree) && requireNamespace("ape", quietly = TRUE)) {
    common_tips <- intersect(tree$tip.label, keep_taxa)
    if (length(common_tips) > 1) {
      new_tree <- ape::keep.tip(tree, common_tips)
    } else if (length(common_tips) == 1) {
      new_tree <- tree
      new_tree$tip.label <- common_tips
      new_tree$edge <- matrix(c(2L, 1L), 1, 2)
      new_tree$edge.length <- 1
      new_tree$Nnode <- 1L
      class(new_tree) <- "phylo"
    }
  }
  
  out <- tb
  attr(out, "assays")    <- new_assays
  attr(out, "tax_table") <- filtered_tax
  attr(out, "phy_tree")  <- new_tree
  out
}

#' Validate Integrity of a tidy_microbiome Object
#'
#' @param tb An object to validate.
#' @param verbose Logical. If `TRUE`, prints a confirmation message upon successful validation.
#'   Defaults to `FALSE`.
#' @return `TRUE` if valid, otherwise throws an error.
#' @export
validate_tidy_microbiome <- function(tb, verbose = FALSE) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("Object must inherit from 'tidy_microbiome'.", call. = FALSE)
  }
  if (!"sample_id" %in% colnames(tb)) {
    stop("Sample metadata must contain a 'sample_id' column.", call. = FALSE)
  }
  
  sample_ids <- as.character(tb$sample_id)
  if (anyDuplicated(sample_ids)) {
    stop("Duplicate 'sample_id' values found in sample metadata.", call. = FALSE)
  }
  
  assays <- attr(tb, "assays")
  if (!is.list(assays) || length(assays) == 0) {
    stop("`tb` must contain an 'assays' list attribute with at least one assay.", call. = FALSE)
  }
  
  for (nm in names(assays)) {
    mat <- assays[[nm]]
    if (!is.matrix(mat)) {
      stop(sprintf("Assay '%s' must be a matrix.", nm), call. = FALSE)
    }
    if (ncol(mat) != length(sample_ids) || !identical(colnames(mat), sample_ids)) {
      stop(sprintf("Sample names in assay '%s' do not match metadata 'sample_id's.", nm), call. = FALSE)
    }
  }
  
  tax_df <- attr(tb, "tax_table")
  if (!is.null(tax_df)) {
    if (!"taxon_id" %in% colnames(tax_df)) {
      stop("`tax_table` must contain a 'taxon_id' column.", call. = FALSE)
    }
    taxon_ids <- as.character(tax_df$taxon_id)
    first_mat <- assays[[1]]
    if (nrow(first_mat) != length(taxon_ids) || !identical(rownames(first_mat), taxon_ids)) {
      stop("Taxon names in assay do not match 'taxon_id' in tax_table.", call. = FALSE)
    }
  }
  
  tree <- attr(tb, "phy_tree")
  if (!is.null(tree)) {
    if (!inherits(tree, "phylo")) {
      stop("`phy_tree` attribute must be of class 'phylo'.", call. = FALSE)
    }
  }
  
  if (verbose) {
    cat("Container validation passed: all assays, taxonomy, and tree tips synchronized.\n")
  }
  invisible(TRUE)
}


#' @export
print.tidy_microbiome <- function(x, ...) {
  n_samp <- nrow(x)
  assays <- attr(x, "assays")
  counts <- assays$counts
  n_taxa <- if (!is.null(counts)) nrow(counts) else 0
  assay_names <- paste(names(assays), collapse = ", ")
  tax_df <- attr(x, "tax_table")
  tax_ranks <- setdiff(colnames(tax_df), "taxon_id")
  rank_str <- if (length(tax_ranks) > 0) paste(tax_ranks, collapse = " > ") else "None"
  has_tree <- !is.null(attr(x, "phy_tree"))

  cat(cli_rule(paste0("tidy_microbiome [", n_samp, " samples \u00d7 ", n_taxa, " taxa]")), "\n")
  cat("  \u2022 Assays:        ", assay_names, "\n")
  cat("  \u2022 Taxonomy ranks:", rank_str, "\n")
  cat("  \u2022 Tree:          ", if (has_tree) "Present" else "None", "\n")
  cat(cli_rule("Sample Metadata"), "\n")

  # Print as tibble
  print(tibble::as_tibble(x), ...)
  invisible(x)
}

cli_rule <- function(title = "") {
  chars <- paste(rep("\u2500", max(2, 50 - nchar(title))), collapse = "")
  paste0("\u2500\u2500 ", title, " ", chars)
}

#' Summary of a tidy_microbiome Object
#'
#' Computes comprehensive ecological, sequencing, and taxonomic summary statistics
#' for a `tidy_microbiome` object.
#'
#' @param object A `tidy_microbiome` object.
#' @param ... Additional arguments passed to methods (currently unused).
#' @return An object of class `summary_tidy_microbiome` containing:
#'   \item{n_samples}{Total number of samples.}
#'   \item{n_taxa}{Total number of taxa.}
#'   \item{depth_summary}{Named vector of sequencing depth statistics.}
#'   \item{sparsity}{Global sparsity proportion (fraction of zeros in counts matrix).}
#'   \item{assays}{Vector of assay names present and dimensions.}
#'   \item{tax_ranks}{Named integer vector of distinct taxa at each taxonomic rank.}
#'   \item{tree_info}{List containing tree presence, number of tips, internal nodes, and rooted status.}
#'   \item{metadata_cols}{Character vector of sample metadata column names and data types.}
#' @export
summary.tidy_microbiome <- function(object, ...) {
  n_samp <- nrow(object)
  assays <- attr(object, "assays")
  counts <- assays$counts
  n_taxa <- if (!is.null(counts)) nrow(counts) else 0

  if (!is.null(counts) && ncol(counts) > 0) {
    depths <- colSums(counts)
    depth_stats <- c(
      Min = min(depths),
      Q25 = stats::quantile(depths, 0.25, names = FALSE),
      Median = stats::median(depths),
      Mean = mean(depths),
      Q75 = stats::quantile(depths, 0.75, names = FALSE),
      Max = max(depths)
    )
    sparsity <- sum(counts == 0) / (nrow(counts) * ncol(counts))
  } else {
    depth_stats <- NULL
    sparsity <- NA_real_
  }

  tax_df <- attr(object, "tax_table")
  tax_ranks <- setdiff(colnames(tax_df), "taxon_id")
  rank_counts <- if (length(tax_ranks) > 0 && !is.null(tax_df)) {
    vapply(tax_ranks, function(r) length(unique(stats::na.omit(tax_df[[r]]))), integer(1))
  } else {
    integer(0)
  }

  tree <- attr(object, "phy_tree")
  tree_info <- if (!is.null(tree) && inherits(tree, "phylo")) {
    is_rooted <- if (requireNamespace("ape", quietly = TRUE)) ape::is.rooted(tree) else TRUE
    list(
      present = TRUE,
      n_tips = length(tree$tip.label),
      n_nodes = tree$Nnode,
      is_rooted = is_rooted
    )
  } else {
    list(present = FALSE)
  }

  meta_df <- tibble::as_tibble(object)
  meta_dict <- vapply(meta_df, function(col) paste(class(col), collapse = "/"), character(1))

  res <- list(
    n_samples = n_samp,
    n_taxa = n_taxa,
    depth_summary = depth_stats,
    sparsity = sparsity,
    assays = names(assays),
    tax_ranks = rank_counts,
    tree_info = tree_info,
    metadata_cols = meta_dict
  )
  class(res) <- "summary_tidy_microbiome"
  res
}

#' @export
print.summary_tidy_microbiome <- function(x, ...) {
  cat(cli_rule(paste0("tidy_microbiome Summary [", x$n_samples, " samples \u00d7 ", x$n_taxa, " taxa]")), "\n")
  cat("  \u2022 Assays:        ", paste(x$assays, collapse = ", "), "\n")
  cat("  \u2022 Sparsity:      ", sprintf("%.1f%% (zeros in count matrix)", x$sparsity * 100), "\n")

  if (!is.null(x$depth_summary)) {
    cat(cli_rule("Sequencing Depth"), "\n")
    ds <- x$depth_summary
    cat(sprintf("    Min: %-8.0f  Q25: %-8.0f  Median: %-8.0f\n", ds["Min"], ds["Q25"], ds["Median"]))
    cat(sprintf("    Mean: %-7.1f  Q75: %-8.0f  Max:    %-8.0f\n", ds["Mean"], ds["Q75"], ds["Max"]))
  }

  if (length(x$tax_ranks) > 0) {
    cat(cli_rule("Taxonomic Hierarchy"), "\n")
    rank_str <- paste(sprintf("%s: %d", names(x$tax_ranks), x$tax_ranks), collapse = " | ")
    cat("    ", rank_str, "\n")
  }

  if (x$tree_info$present) {
    cat(cli_rule("Phylogenetic Tree"), "\n")
    cat(sprintf("    %d tips, %d internal nodes (%s)\n",
                x$tree_info$n_tips, x$tree_info$n_nodes,
                if (x$tree_info$is_rooted) "rooted" else "unrooted"))
  }

  cat(cli_rule(paste0("Sample Metadata (", length(x$metadata_cols), " variables)")), "\n")
  vars_str <- paste(sprintf("%s <%s>", names(x$metadata_cols), x$metadata_cols), collapse = ", ")
  cat("    ", strwrap(vars_str, width = 60, exdent = 4), sep = "\n")

  invisible(x)
}
