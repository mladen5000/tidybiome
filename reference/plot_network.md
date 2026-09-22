# Plot Microbial Co-Occurrence Network

Plot Microbial Co-Occurrence Network

## Usage

``` r
plot_network(net, color_by = "Phylum", min_degree = 1, palette = "tidybiome")
```

## Arguments

- net:

  An object produced by
  [`calc_network()`](https://tidybiome.org/reference/calc_network.md).

- color_by:

  Taxonomy column to color nodes by (e.g. `"Phylum"`). Defaults to
  `"Phylum"`.

- min_degree:

  Minimum degree for a node to be displayed (default: 1).

- palette:

  Palette name: `"tidybiome"` or `"nature"`.

## Value

A
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
