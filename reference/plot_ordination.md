# Plot Microbiome Ordination and Compositional Biplots

Generates elegant 2D ordination plots (RPCA, PCoA, PCA) with automatic
variance explained axis annotations, confidence ellipses, and optional
taxon loading biplot vectors.

## Usage

``` r
plot_ordination(
  tb,
  method = c("rpca", "pcoa", "pca"),
  color = NULL,
  shape = NULL,
  ellipse = TRUE,
  biplot = FALSE,
  top_taxa = 5,
  palette = "tidybiome"
)
```

## Arguments

- tb:

  A `tidy_microbiome` object. Run
  [`calc_ordination()`](https://tidybiome.org/reference/calc_ordination.md)
  beforehand or specify `method`.

- method:

  Character string specifying ordination method (`"rpca"`, `"pcoa"`, or
  `"pca"`).

- color:

  Optional character string naming the metadata column to color points
  by.

- shape:

  Optional character string naming the metadata column for point shapes.

- ellipse:

  Logical. If `TRUE` and `color` is provided, draws 95% confidence
  ellipses around groups.

- biplot:

  Logical. If `TRUE` (and `method = "rpca"`), overlays top taxon loading
  vectors as arrows.

- top_taxa:

  Integer number of top driving taxa to annotate if `biplot = TRUE`
  (default: 5).

- palette:

  Palette name: `"tidybiome"` or `"nature"`.

## Value

A
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
counts <- matrix(c(100, 20, 5, 2, 5, 10, 80, 70, 1, 1, 50, 40), nrow = 3, byrow = TRUE,
                 dimnames = list(c("Tax1", "Tax2", "Tax3"), c("S1", "S2", "S3", "S4")))
sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"))
tb <- tidy_microbiome(counts, sample_data)
tb <- calc_ordination(tb, method = "rpca")
p <- plot_ordination(tb, method = "rpca", color = "group")
```
