#' Coerce Objects to tidy_microbiome
#'
#' @param x An object to coerce (e.g. from `phyloseq`, `TreeSummarizedExperiment`, or `matrix`).
#' @param ... Additional arguments passed to methods.
#' @return A `tidy_microbiome` object.
#' @export
as_tidybiome <- function(x, ...) {
  UseMethod("as_tidybiome")
}

#' @export
as_tidybiome.default <- function(x, sample_data = NULL, tax_table = NULL, phy_tree = NULL, ...) {
  if (inherits(x, "tidybiome_vegan")) {
    return(from_vegan(comm = x$comm, env = x$env, tax_table = tax_table, phy_tree = phy_tree, ...))
  }
  if (inherits(x, "phyloseq")) {
    return(from_phyloseq(x))
  }
  if (inherits(x, "TreeSummarizedExperiment") || inherits(x, "SummarizedExperiment")) {
    return(from_tse(x))
  }
  if (is.matrix(x) || is.data.frame(x)) {
    return(tidy_microbiome(as.matrix(x), sample_data = sample_data, tax_table = tax_table, phy_tree = phy_tree, ...))
  }
  stop(sprintf("Cannot coerce object of class '%s' to tidy_microbiome.", class(x)[1]), call. = FALSE)
}

from_phyloseq <- function(ps) {
  # Extract otu_table
  otu <- as(ps@otu_table, "matrix")
  if (!ps@otu_table@taxa_are_rows) {
    otu <- t(otu)
  }

  # Extract sample_data
  sdata <- if (!is.null(ps@sam_data)) as.data.frame(ps@sam_data) else NULL

  # Extract tax_table
  ttable <- if (!is.null(ps@tax_table)) as.data.frame(ps@tax_table) else NULL

  # Extract tree
  ptree <- if (!is.null(ps@phy_tree)) ps@phy_tree else NULL

  tidy_microbiome(counts = otu, sample_data = sdata, tax_table = ttable, phy_tree = ptree)
}

from_tse <- function(tse) {
  assay_names <- SummarizedExperiment::assayNames(tse)
  primary_assay <- if ("counts" %in% assay_names) "counts" else assay_names[1]
  counts <- as.matrix(SummarizedExperiment::assay(tse, primary_assay))

  sdata <- as.data.frame(SummarizedExperiment::colData(tse))
  ttable <- as.data.frame(SummarizedExperiment::rowData(tse))

  ptree <- tryCatch({
    if ("rowTree" %in% names(methods::getClass(class(tse))@slots)) {
      TreeSummarizedExperiment::rowTree(tse)
    } else {
      NULL
    }
  }, error = function(e) NULL)

  tb <- tidy_microbiome(counts = counts, sample_data = sdata, tax_table = ttable, phy_tree = ptree)

  # Transfer additional assays if present
  other_assays <- setdiff(assay_names, primary_assay)
  for (oa in other_assays) {
    attr(tb, "assays")[[oa]] <- as.matrix(SummarizedExperiment::assay(tse, oa))
  }
  tb
}

#' Convert tidy_microbiome to phyloseq Object
#'
#' @param tb A `tidy_microbiome` object.
#' @return A `phyloseq::phyloseq` object.
#' @export
to_phyloseq <- function(tb) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required for `to_phyloseq()`. Please install it.", call. = FALSE)
  }

  counts <- attr(tb, "assays")$counts
  otu <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)

  sdata_df <- as.data.frame(tibble::as_tibble(tb))
  rownames(sdata_df) <- sdata_df$sample_id
  sdata <- phyloseq::sample_data(sdata_df)

  tax_df <- as.data.frame(attr(tb, "tax_table"))
  rownames(tax_df) <- tax_df$taxon_id
  tax_df$taxon_id <- NULL
  ttable <- phyloseq::tax_table(as.matrix(tax_df))

  ptree <- attr(tb, "phy_tree")
  if (!is.null(ptree)) {
    phyloseq::phyloseq(otu, sdata, ttable, ptree)
  } else {
    phyloseq::phyloseq(otu, sdata, ttable)
  }
}

#' Convert tidy_microbiome to TreeSummarizedExperiment Object
#'
#' @param tb A `tidy_microbiome` object.
#' @return A `TreeSummarizedExperiment::TreeSummarizedExperiment` object.
#' @export
to_tse <- function(tb) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    stop("Package 'TreeSummarizedExperiment' is required for `to_tse()`. Please install it.", call. = FALSE)
  }

  assays <- attr(tb, "assays")
  sdata_df <- as.data.frame(tibble::as_tibble(tb))
  rownames(sdata_df) <- sdata_df$sample_id

  tax_df <- as.data.frame(attr(tb, "tax_table"))
  rownames(tax_df) <- tax_df$taxon_id

  ptree <- attr(tb, "phy_tree")
  TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays  = assays,
    colData = sdata_df,
    rowData = tax_df,
    rowTree = ptree
  )
}

#' Convert tidy_microbiome to vegan Community and Environmental Format
#'
#' @param tb A `tidy_microbiome` object.
#' @param assay Name of assay to extract. Defaults to `"counts"`.
#' @return An S3 object of class `tidybiome_vegan` containing:
#'   \itemize{
#'     \item `comm`: Community matrix with samples in rows and taxa in columns.
#'     \item `env`: Sample metadata tibble.
#'   }
#' @export
to_vegan <- function(tb, assay = "counts") {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  assays <- attr(tb, "assays")
  if (!assay %in% names(assays)) {
    stop(sprintf("Assay '%s' not found in `tb`.", assay), call. = FALSE)
  }

  mat <- assays[[assay]]
  # Vegan expects samples in rows and taxa in columns
  comm <- t(mat)

  env_df <- as.data.frame(tibble::as_tibble(tb))
  rownames(env_df) <- env_df$sample_id

  res <- list(comm = comm, env = tibble::as_tibble(env_df), env_df = env_df)
  class(res) <- c("tidybiome_vegan", "list")
  res
}

#' @export
print.tidybiome_vegan <- function(x, ...) {
  cat(sprintf("-- tidybiome_vegan [%d samples x %d taxa] --\n", nrow(x$comm), ncol(x$comm)))
  cat(sprintf("  * Community matrix: %d rows x %d columns\n", nrow(x$comm), ncol(x$comm)))
  cat(sprintf("  * Environmental metadata: %d variables (%s)\n", ncol(x$env), paste(head(names(x$env), 4), collapse = ", ")))
  invisible(x)
}

#' Construct tidy_microbiome from vegan Community Format
#'
#' @param comm Community matrix or data frame with samples in rows and taxa in columns.
#' @param env Optional sample metadata data frame.
#' @param tax_table Optional taxonomy data frame.
#' @param phy_tree Optional phylogenetic tree.
#' @return A `tidy_microbiome` object.
#' @export
from_vegan <- function(comm, env = NULL, tax_table = NULL, phy_tree = NULL) {
  if (!is.matrix(comm) && !is.data.frame(comm)) {
    stop("`comm` must be a matrix or data frame.", call. = FALSE)
  }
  comm_mat <- as.matrix(comm)
  # Transpose to tidy_microbiome convention: taxa in rows, samples in columns
  counts <- t(comm_mat)

  if (is.null(rownames(comm_mat))) {
    colnames(counts) <- paste0("Sample_", seq_len(ncol(counts)))
  }
  if (is.null(colnames(comm_mat))) {
    rownames(counts) <- paste0("Taxon_", seq_len(nrow(counts)))
  }

  tidy_microbiome(counts = counts, sample_data = env, tax_table = tax_table, phy_tree = phy_tree)
}

#' @export
as_tidybiome.tidybiome_vegan <- function(x, ...) {
  from_vegan(comm = x$comm, env = x$env, ...)
}

#' Convert tidy_microbiome to phyloseq Object (Alias)
#'
#' @rdname to_phyloseq
#' @export
as_phyloseq <- function(tb) {
  to_phyloseq(tb)
}

#' Convert tidy_microbiome to TreeSummarizedExperiment / mia Object (Aliases)
#'
#' @rdname to_tse
#' @export
as_mia <- function(tb) {
  to_tse(tb)
}

#' @rdname to_tse
#' @export
to_mia <- function(tb) {
  to_tse(tb)
}

#' Get Phylogenetic Tree from tidy_microbiome
#'
#' @param tb A `tidy_microbiome` object.
#' @return An `ape::phylo` object, or `NULL` if none is attached.
#' @export
get_tree <- function(tb) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  attr(tb, "phy_tree")
}

#' Attach and Validate Phylogenetic Tree to tidy_microbiome
#'
#' @param tb A `tidy_microbiome` object.
#' @param tree An `ape::phylo` object.
#' @param prune Logical; if `TRUE`, automatically prunes the tree to match taxon IDs present in `tb`. Defaults to `TRUE`.
#' @return A `tidy_microbiome` object with the attached tree.
#' @export
set_tree <- function(tb, tree, prune = TRUE) {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }
  if (!inherits(tree, "phylo")) {
    stop("`tree` must be an object of class 'phylo' (from package ape).", call. = FALSE)
  }

  taxa <- attr(tb, "tax_table")$taxon_id
  common_tips <- intersect(tree$tip.label, taxa)

  if (length(common_tips) == 0) {
    stop("No tree tip labels match the taxon IDs in `tb`.", call. = FALSE)
  }

  if (prune && length(tree$tip.label) > length(common_tips)) {
    if (requireNamespace("ape", quietly = TRUE)) {
      tree <- ape::keep.tip(tree, common_tips)
    }
  }

  attr(tb, "phy_tree") <- tree
  tb
}


