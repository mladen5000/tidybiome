# Calculate Prevalence and Abundance Metrics Across Taxa

Evaluates the proportion and number of samples where each taxon exceeds
a specified detection threshold, along with mean and median abundances.

## Usage

``` r
calc_prevalence(tb, assay = "counts", detection = 0, sort = TRUE)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- assay:

  Assay to evaluate. Defaults to `"counts"`.

- detection:

  Abundance threshold above which a taxon is considered detected.
  Defaults to `0`.

- sort:

  Logical; if `TRUE`, sorts taxa in descending order of prevalence.
  Defaults to `TRUE`.

## Value

A tibble with columns:

- `taxon_id`: Unique identifier for each taxon.

- `prevalence`: Proportion of samples where taxon was detected (\\\[0,
  1\]\\).

- `prevalence_n`: Count of samples where taxon was detected.

- `mean_abundance`: Average abundance across all samples.

- `median_abundance`: Median abundance across all samples.

Joined with available taxonomic ranks from `tax_table`.
