# Extract Abundance Data in Tidy Format

Extract Abundance Data in Tidy Format

## Usage

``` r
tidy_abundance(
  tb,
  assay = "counts",
  rank = NULL,
  long = TRUE,
  taxa = NULL,
  samples = NULL,
  include_metadata = TRUE,
  include_taxonomy = TRUE
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- assay:

  Character string naming the assay to extract (e.g. `"counts"`,
  `"relabundance"`, `"rclr"`). Defaults to `"counts"`.

- rank:

  Optional taxonomic rank to aggregate by before extraction (e.g.
  `"Genus"`, `"Phylum"`).

- long:

  Logical. If `TRUE` (default), returns a fully denormalized long tibble
  containing `sample_id`, `taxon_id`, `abundance`, sample metadata, and
  taxonomy. If `FALSE`, returns the abundance matrix directly.

- taxa:

  Optional character vector of taxon IDs to filter down before
  long-format expansion.

- samples:

  Optional character vector of sample IDs to filter down before
  long-format expansion.

- include_metadata:

  Logical; whether to join sample metadata in long format (default:
  `TRUE`).

- include_taxonomy:

  Logical; whether to join taxonomy metadata in long format (default:
  `TRUE`).

## Value

A
[tibble::tbl_df](https://tibble.tidyverse.org/reference/tbl_df-class.html)
if `long = TRUE`, or a numeric matrix if `long = FALSE`.
