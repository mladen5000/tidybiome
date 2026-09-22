# Run PERMANOVA via vegan::adonis2

Executes Permutational Multivariate Analysis of Variance using
[`vegan::adonis2`](https://vegandevs.github.io/vegan/reference/adonis.html)
on a `tidy_microbiome` object and returns a clean, tidy tibble.

## Usage

``` r
run_permanova(
  tb,
  formula,
  assay = "counts",
  method = "bray",
  permutations = 999,
  by = "terms",
  ...
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- formula:

  A model formula specifying predictors from sample metadata (e.g.
  `~ treatment + diet`). If the left-hand side is omitted, the community
  matrix is automatically supplied.

- assay:

  Assay name to use. Defaults to `"counts"`.

- method:

  Distance metric to pass to
  [`vegan::adonis2`](https://vegandevs.github.io/vegan/reference/adonis.html)
  (e.g. `"bray"`, `"jaccard"`). Defaults to `"bray"`.

- permutations:

  Number of permutations. Defaults to `999`.

- by:

  How terms are assessed: `"terms"`, `"margin"`, or `NULL`. Defaults to
  `"terms"`.

- ...:

  Additional arguments passed to
  [`vegan::adonis2`](https://vegandevs.github.io/vegan/reference/adonis.html).

## Value

A tidy tibble with columns:

- `term`: Predictor name, Residual, or Total.

- `df`: Degrees of freedom.

- `sum_sq`: Sum of squares.

- `r2`: Coefficient of determination (\\R^2\\).

- `f_stat`: Pseudo-F statistic.

- `p_value`: Permutation p-value.
