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
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("`counts` matrix must have both rownames (taxa) and colnames (samples).", call. = FALSE)
  }

  sample_names <- colnames(counts)
  taxon_names  <- rownames(counts)

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
