# Calculate Ecological Divergence Relative to a Reference State

Evaluates the dissimilarity of each sample's community profile against a
reference baseline (such as a healthy control group, cohort median
profile, or designated baseline sample).

## Usage

``` r
calc_divergence(
  tb,
  reference = "median",
  method = "bray",
  assay = "counts",
  add_to_metadata = TRUE
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- reference:

  Reference baseline to compare against. Can be:

  - `"median"`: Cohort median abundance profile across all samples.

  - `"mean"`: Cohort mean abundance profile.

  - A named list specifying a metadata subset (e.g.
    `list(treatment = "Control")`).

  - A character string identifying a specific reference `sample_id`.

- method:

  Dissimilarity metric: `"bray"` (Bray-Curtis), `"jsd"` (Jensen-Shannon
  divergence), or `"euclidean"`. Defaults to `"bray"`.

- assay:

  Assay to evaluate. Defaults to `"counts"`.

- add_to_metadata:

  Logical; if `TRUE`, appends column `divergence` to `tb`'s sample
  metadata. If `FALSE`, returns a tidy tibble. Defaults to `TRUE`.

## Value

A `tidy_microbiome` object with `divergence` in metadata (if
`add_to_metadata = TRUE`), or a tibble with `sample_id` and
`divergence`.
