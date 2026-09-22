# Construct tidy_microbiome from vegan Community Format

Construct tidy_microbiome from vegan Community Format

## Usage

``` r
from_vegan(comm, env = NULL, tax_table = NULL, phy_tree = NULL)
```

## Arguments

- comm:

  Community matrix or data frame with samples in rows and taxa in
  columns.

- env:

  Optional sample metadata data frame.

- tax_table:

  Optional taxonomy data frame.

- phy_tree:

  Optional phylogenetic tree.

## Value

A `tidy_microbiome` object.
