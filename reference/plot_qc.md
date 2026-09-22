# Plot Sample Quality Control Diagnostics

Generates publication-ready diagnostic scatter plots of sequencing
library depth versus observed taxonomic richness, with optional
threshold lines, log scaling, and metadata grouping.

## Usage

``` r
plot_qc(
  tb,
  x = "total_reads",
  y = "n_features",
  color_by = NULL,
  threshold_x = NULL,
  threshold_y = NULL,
  log_x = (x == "total_reads"),
  log_y = FALSE,
  palette = "tidybiome"
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- x:

  Metric for the x-axis: `"total_reads"`, `"n_features"`, `"sparsity"`,
  or `"top_taxon_share"`. Defaults to `"total_reads"`.

- y:

  Metric for the y-axis. Defaults to `"n_features"`.

- color_by:

  Optional metadata column to color points by.

- threshold_x:

  Optional numeric threshold for vertical dashed cutoff line.

- threshold_y:

  Optional numeric threshold for horizontal dashed cutoff line.

- log_x:

  Logical; whether to log10-transform the x-axis (default: `TRUE` if
  `x == "total_reads"`).

- log_y:

  Logical; whether to log10-transform the y-axis (default: `FALSE`).

- palette:

  Palette name: `"tidybiome"` or `"nature"`.

## Value

A
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
data(gut_microbiome)
plot_qc(gut_microbiome, color_by = "treatment", threshold_x = 5000)
```
