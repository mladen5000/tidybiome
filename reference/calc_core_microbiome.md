# Identify Core Microbiome Members

Detects the core microbiome members across samples using either modern
data-driven inflection curve analysis (identifying natural mathematical
breakpoints in the prevalence spectrum) or classic static threshold
filtering.

## Usage

``` r
calc_core_microbiome(
  tb,
  method = c("threshold", "inflection"),
  min_prevalence = 0.8,
  min_abundance = 0.001,
  assay = "counts"
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- method:

  Character string specifying the core detection method:

  - `"inflection"`: Data-driven inflection curve analysis. Evaluates the
    sorted prevalence distribution across taxa and finds the natural
    elbow/inflection point.

  - `"threshold"`: Static filter on minimum prevalence and relative
    abundance.

- min_prevalence:

  Numeric minimum prevalence threshold (between 0 and 1) for
  `"threshold"` method. Defaults to 0.8 (present in at least 80% of
  samples).

- min_abundance:

  Numeric minimum relative abundance threshold for a taxon to be
  considered present. Defaults to 0.001 (0.1%).

- assay:

  Character string naming the abundance assay. Defaults to `"counts"`.

## Value

A
[tibble::tbl_df](https://tibble.tidyverse.org/reference/tbl_df-class.html)
with one row per taxon containing `taxon_id`, `prevalence`,
`mean_abundance`, `is_core`, and taxonomic annotations.

## Examples

``` r
counts <- matrix(c(50, 40, 60, 55, 45, 0, 0, 10, 0, 0), nrow = 2, byrow = TRUE,
                 dimnames = list(c("CoreTaxon", "RareTaxon"), paste0("S", 1:5)))
tb <- tidy_microbiome(counts)
core_df <- calc_core_microbiome(tb, method = "threshold", min_prevalence = 0.8)
```
