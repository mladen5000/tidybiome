# Plot Distance-Based Redundancy Analysis (db-RDA) Biplot

Plot Distance-Based Redundancy Analysis (db-RDA) Biplot

## Usage

``` r
plot_dbrda(dbrda_res, color = NULL, shape = NULL, palette = "tidybiome")
```

## Arguments

- dbrda_res:

  An object produced by
  [`run_dbrda()`](https://tidybiome.org/reference/run_dbrda.md).

- color:

  Optional metadata variable to color sample points.

- shape:

  Optional metadata variable for point shapes.

- palette:

  Palette name: `"tidybiome"` or `"nature"`.

## Value

A
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
