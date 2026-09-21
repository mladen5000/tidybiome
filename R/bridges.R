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
