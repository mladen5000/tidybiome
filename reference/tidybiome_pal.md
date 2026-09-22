# Curated Color Palettes for tidybiome

Curated Color Palettes for tidybiome

## Usage

``` r
tidybiome_pal(n = NULL, palette = c("tidybiome", "nature"))

scale_color_tidybiome(palette = "tidybiome", ...)

scale_fill_tidybiome(palette = "tidybiome", ...)
```

## Arguments

- n:

  Number of colors desired.

- palette:

  Palette name: `"tidybiome"` (Okabe-Ito colorblind-safe) or `"nature"`
  (Editorial journal palette).

- ...:

  Arguments passed to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html).

## Value

A character vector of hexadecimal color codes.
