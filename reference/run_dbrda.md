# Distance-Based Redundancy Analysis via vegan::dbrda

Performs constrained ordination (db-RDA) to model microbial community
variation explained by environmental predictors or experimental
conditions.

## Usage

``` r
run_dbrda(tb, formula, assay = "counts", distance = "bray", ...)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- formula:

  A model formula specifying predictors from sample metadata (e.g.
  `~ treatment + age`).

- assay:

  Assay name to use. Defaults to `"counts"`.

- distance:

  Distance metric to compute (e.g. `"bray"`, `"raitchison"`). Defaults
  to `"bray"`.

- ...:

  Additional arguments passed to
  [`vegan::dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.html).

## Value

An S3 object of class `tidybiome_dbrda` containing:

- `samples`: Tibble of sample coordinates along db-RDA axes joined with
  metadata.

- `biplot`: Tibble of constraint vector loadings (continuous
  predictors).

- `centroids`: Tibble of factor centroids (categorical predictors).

- `variance_explained`: Proportion of variance explained by constrained
  axes.

- `model`: The underlying
  [`vegan::dbrda`](https://vegandevs.github.io/vegan/reference/dbrda.html)
  model object.
