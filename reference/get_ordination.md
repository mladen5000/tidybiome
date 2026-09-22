# Extract Ordination Results from tidy_microbiome

Extract Ordination Results from tidy_microbiome

## Usage

``` r
get_ordination(tb, method = NULL)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- method:

  Optional character string naming the ordination method to extract. If
  `NULL`, returns the most recently calculated ordination.

## Value

A list containing:

- `samples`:
  [tibble::tbl_df](https://tibble.tidyverse.org/reference/tbl_df-class.html)
  of sample coordinates and sample metadata.

- `taxa`:
  [tibble::tbl_df](https://tibble.tidyverse.org/reference/tbl_df-class.html)
  of taxon loadings (if available) and taxonomy.

- `variance_explained`: Numeric vector of variance explained by each
  axis.
