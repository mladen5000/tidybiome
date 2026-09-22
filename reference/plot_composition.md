# Plot Taxonomic Composition with Intelligent Top-N Grouping

Generates publication-ready stacked bar charts of taxonomic composition.
Automatically aggregates less abundant taxa into an aesthetically styled
`"Other"` category at the base of the plot, eliminating cluttered
legends.

## Usage

``` r
plot_composition(
  tb,
  rank = "Genus",
  top_n = 8,
  group_by = NULL,
  palette = "tidybiome"
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- rank:

  Character string specifying taxonomic rank (e.g. `"Genus"`,
  `"Phylum"`). If `NULL`, uses ASV/feature IDs.

- top_n:

  Number of top most abundant taxa to display individually (default: 8).

- group_by:

  Optional character string specifying a sample metadata column to facet
  by.

- palette:

  Palette name: `"tidybiome"` or `"nature"`.

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
p <- plot_composition(tb, rank = "Genus", top_n = 2)
```
