# Distance-based Test for Homogeneity (DTH)

Implements the 2025 bioRxiv Distance-based Test for Homogeneity (DTH),
comparing the empirical distributions of within-group multivariate
distances using the 1D Wasserstein optimal transport metric and
permutation testing. Unlike `betadisper`, DTH directly evaluates full
distance distributions rather than assuming spherical variance around
centroids.

## Usage

``` r
dth_test(tb, group, metric = "raitchison", n_perm = 499)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- group:

  Character string naming the two-level grouping column in sample
  metadata.

- metric:

  Character string naming the distance metric (default: `"raitchison"`).

- n_perm:

  Number of Monte Carlo permutations (default: 499).

## Value

A
[tibble::tbl_df](https://tibble.tidyverse.org/reference/tbl_df-class.html)
with `statistic` (Wasserstein distance between dispersion curves),
`p_value`, `n_perm`, and `metric`.
