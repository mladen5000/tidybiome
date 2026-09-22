# Identify Dominant Taxa per Sample and Cohort-Wide

Finds the taxon with highest relative abundance in each sample and
computes cohort-wide dominance.

## Usage

``` r
calc_dominant(tb, assay = "counts", rank = NULL, add_to_metadata = TRUE)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- assay:

  Assay to evaluate. Defaults to `"counts"`.

- rank:

  Optional taxonomic rank to agglomerate before determining dominance.
  Defaults to `NULL`.

- add_to_metadata:

  Logical; if `TRUE`, adds column `dominant_taxa` directly to `tb`'s
  sample metadata. If `FALSE`, returns a summary tibble of dominant taxa
  across the cohort. Defaults to `TRUE`.

## Value

A `tidy_microbiome` object with `dominant_taxa` in metadata (if
`add_to_metadata = TRUE`), or a tibble summarizing dominance frequencies
across samples.
