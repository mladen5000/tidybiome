# Test for Homogeneity of Multivariate Dispersions via vegan::betadisper

Quantifies and tests within-group multivariate dispersion using
[`vegan::betadisper`](https://vegandevs.github.io/vegan/reference/betadisper.html)
and permutation tests, returning tidy summary tibbles.

## Usage

``` r
run_betadisper(
  tb,
  group,
  assay = "counts",
  method = "bray",
  permutations = 999,
  ...
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- group:

  Character string naming the grouping variable in sample metadata.

- assay:

  Assay name to use. Defaults to `"counts"`.

- method:

  Distance metric to compute (e.g. `"bray"`, `"jaccard"`). Defaults to
  `"bray"`.

- permutations:

  Number of permutations for dispersion test. Defaults to `999`.

- ...:

  Additional arguments passed to
  [`vegan::betadisper`](https://vegandevs.github.io/vegan/reference/betadisper.html).

## Value

An S3 object of class `tidybiome_betadisper` containing:

- `sample_distances`: Tibble with `sample_id`, `group`, and
  `distance_to_centroid`.

- `group_summary`: Tibble with per-group counts, mean distances, and
  standard errors.

- `test`: Tibble with F-statistic, degrees of freedom, and permutation
  p-value.
