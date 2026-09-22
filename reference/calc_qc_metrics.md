# Calculate Sample-Level Quality Control (QC) Metrics

Computes key quality control and diagnostic statistics for each sample
in a `tidy_microbiome` object, including total library size, observed
feature richness, sample-level sparsity, dominance of the top taxon, and
Shannon entropy.

## Usage

``` r
calc_qc_metrics(tb, assay = "counts", augment = TRUE)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- assay:

  Name of the assay to compute metrics on (defaults to `"counts"`).

- augment:

  Logical; if `TRUE` (default), the QC metrics are appended directly as
  columns to the sample metadata within the returned `tidy_microbiome`
  object, enabling seamless filtering in dplyr pipelines. If `FALSE`,
  returns a tidy
  [tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
  of QC metrics.

## Value

If `augment = TRUE`, an updated `tidy_microbiome` object with added
columns:

- qc_total_reads:

  Total sum of read counts in the sample.

- qc_n_features:

  Number of features (taxa) observed with count \> 0.

- qc_sparsity:

  Proportion of unobserved features in the sample (fraction of zeros).

- qc_top_taxon_share:

  Proportion of total sample library represented by the single most
  abundant taxon.

- qc_shannon:

  Shannon diversity / entropy index calculated from sample proportions.

If `augment = FALSE`, returns a tibble with `sample_id` and the above
columns.

## Examples

``` r
data(gut_microbiome)
# Augment metadata with QC metrics and filter high-quality samples
tb_clean <- calc_qc_metrics(gut_microbiome)
tb_filtered <- tb_clean[tb_clean$qc_total_reads >= 5000, ]

# Or extract QC metrics as a tibble
qc_df <- calc_qc_metrics(gut_microbiome, augment = FALSE)
```
