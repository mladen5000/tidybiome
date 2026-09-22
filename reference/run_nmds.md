# Run Non-Metric Multidimensional Scaling via vegan::metaMDS

Run Non-Metric Multidimensional Scaling via vegan::metaMDS

## Usage

``` r
run_nmds(
  tb,
  assay = "counts",
  method = "bray",
  k = 2,
  trymax = 50,
  trace = 0,
  ...
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- assay:

  Assay name to use. Defaults to `"counts"`.

- method:

  Dissimilarity index to pass to
  [`vegan::metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.html).
  Defaults to `"bray"`.

- k:

  Number of dimensions. Defaults to `2`.

- trymax:

  Maximum number of random starts. Defaults to `50`.

- trace:

  Numeric trace output (0 for silent). Defaults to `0`.

- ...:

  Additional arguments passed to
  [`vegan::metaMDS`](https://vegandevs.github.io/vegan/reference/metaMDS.html).

## Value

A tibble of ordination coordinates (`NMDS1`, `NMDS2`, etc.) joined with
sample metadata, with `stress` stored as an attribute.
