# Create a tidy_microbiome Object

Creates a `tidy_microbiome` S3 object that inherits from
[tibble::tbl_df](https://tibble.tidyverse.org/reference/tbl_df-class.html).
The object behaves as a sample metadata tibble on the outside, while
keeping count matrices (assays), taxonomy annotations, phylogenetic
trees, and experiment metadata synchronized under the hood.

## Usage

``` r
tidy_microbiome(
  counts,
  sample_data = NULL,
  tax_table = NULL,
  phy_tree = NULL,
  metadata = list()
)
```

## Arguments

- counts:

  A numeric matrix of counts (taxa in rows, samples in columns).

- sample_data:

  Optional `data.frame` or `tbl_df` of sample metadata. Must include a
  `sample_id` column or rownames matching `colnames(counts)`.

- tax_table:

  Optional `data.frame` or `tbl_df` of taxonomy annotations. Must
  include a `taxon_id` column or rownames matching `rownames(counts)`.

- phy_tree:

  Optional phylogenetic tree (e.g. of class `phylo`).

- metadata:

  Optional list of arbitrary experiment metadata.

## Value

An object of class `tidy_microbiome`, inheriting from `tbl_df`, `tbl`,
and `data.frame`.

## Examples

``` r
counts <- matrix(c(10, 0, 5, 20), nrow = 2,
                 dimnames = list(c("ASV1", "ASV2"), c("S1", "S2")))
tb <- tidy_microbiome(counts)
tb
#> ── tidy_microbiome [2 samples × 2 taxa] ────────────── 
#>   • Assays:         counts 
#>   • Taxonomy ranks: None 
#>   • Tree:           None 
#> ── Sample Metadata ─────────────────────────────────── 
#> # A tibble: 2 × 2
#>   sample_id depth
#>   <chr>     <dbl>
#> 1 S1           10
#> 2 S2           25
```
