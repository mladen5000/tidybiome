# Calculate Unweighted or Weighted UniFrac Distance

Computes phylogenetic UniFrac dissimilarity matrices between samples
using an attached phylogenetic tree
([`ape::phylo`](https://rdrr.io/pkg/ape/man/read.tree.html)). Supports
both qualitative (unweighted) and quantitative (weighted normalized)
UniFrac metrics.

## Usage

``` r
calc_unifrac(tb, weighted = FALSE, normalized = TRUE, assay = "counts")
```

## Arguments

- tb:

  A `tidy_microbiome` object with an attached `phy_tree`.

- weighted:

  Logical. If `TRUE`, computes weighted UniFrac (accounting for relative
  abundance). If `FALSE` (default), computes unweighted UniFrac (based
  on presence/absence).

- normalized:

  Logical. If `TRUE` (default), normalizes weighted UniFrac by total
  tree branch length.

- assay:

  Character string naming the abundance assay to use. Defaults to
  `"counts"`.

## Value

A `dist` object of pairwise UniFrac distances.
