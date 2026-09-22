# Summary of a tidy_microbiome Object

Computes comprehensive ecological, sequencing, and taxonomic summary
statistics for a `tidy_microbiome` object.

## Usage

``` r
# S3 method for class 'tidy_microbiome'
summary(object, ...)
```

## Arguments

- object:

  A `tidy_microbiome` object.

- ...:

  Additional arguments passed to methods (currently unused).

## Value

An object of class `summary_tidy_microbiome` containing:

- n_samples:

  Total number of samples.

- n_taxa:

  Total number of taxa.

- depth_summary:

  Named vector of sequencing depth statistics.

- sparsity:

  Global sparsity proportion (fraction of zeros in counts matrix).

- assays:

  Vector of assay names present and dimensions.

- tax_ranks:

  Named integer vector of distinct taxa at each taxonomic rank.

- tree_info:

  List containing tree presence, number of tips, internal nodes, and
  rooted status.

- metadata_cols:

  Character vector of sample metadata column names and data types.
