# Plot Compositional Microbiome Heatmap

Generates an aesthetic, publication-ready abundance heatmap for top taxa
across samples, featuring hierarchical clustering of samples and taxa,
compositional scaling, and optional metadata grouping.

## Usage

``` r
plot_heatmap(
  tb,
  rank = NULL,
  top_n = 25,
  assay = "counts",
  scale = c("log10", "relabundance", "rclr", "none"),
  cluster_samples = TRUE,
  cluster_taxa = TRUE,
  annotation_col = NULL,
  palette = "viridis"
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- rank:

  Optional taxonomic rank to aggregate by before plotting (e.g.
  `"Genus"`).

- top_n:

  Number of top abundant taxa to display (default: 25).

- assay:

  Assay to extract and plot. Defaults to `"counts"`.

- scale:

  Scaling applied to values: `"log10"`, `"relabundance"`, `"rclr"`, or
  `"none"`.

- cluster_samples:

  Logical. If `TRUE` (default), clusters samples via hierarchical
  clustering.

- cluster_taxa:

  Logical. If `TRUE` (default), clusters taxa via hierarchical
  clustering.

- annotation_col:

  Optional sample metadata column to group or facet samples by.

- palette:

  Palette option: `"viridis"`, `"magma"`, or `"plasma"`.

## Value

A
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
counts <- matrix(c(100, 20, 5, 2, 5, 10, 80, 70), nrow = 2, byrow = TRUE,
                 dimnames = list(c("Tax1", "Tax2"), c("S1", "S2", "S3", "S4")))
sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"))
tax_table <- data.frame(taxon_id = c("Tax1", "Tax2"), Genus = c("Bacteroides", "Prevotella"))
tb <- tidy_microbiome(counts, sample_data, tax_table)
p <- plot_heatmap(tb, top_n = 2)
```
