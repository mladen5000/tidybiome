# Plot Consensus Differential Abundance Volcano Plot

Plots effect size (\\\log_2\text{FC}\\) versus statistical significance
(\\-\log\_{10}(p\_{\text{adj}})\\), sizing points by their multi-engine
consensus agreement score and labeling top biomarker taxa.

## Usage

``` r
plot_da_volcano(da_res, fdr_cutoff = 0.05, fc_cutoff = 1, label_top = 8)
```

## Arguments

- da_res:

  A data frame or tibble produced by
  [`calc_differential_abundance()`](https://tidybiome.org/reference/calc_differential_abundance.md).

- fdr_cutoff:

  False discovery rate significance threshold (default: 0.05).

- fc_cutoff:

  Absolute log2 fold-change cutoff (default: 1.0).

- label_top:

  Number of top significant taxa to label with names (default: 8).

## Value

A
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
counts <- matrix(c(100, 120, 110, 5, 8, 6, 50, 45, 55, 48, 52, 50), nrow = 2, byrow = TRUE,
                 dimnames = list(c("DiffTaxon", "NullTaxon"), paste0("S", 1:6)))
sample_data <- data.frame(sample_id = paste0("S", 1:6), group = c("A", "A", "A", "B", "B", "B"))
tb <- tidy_microbiome(counts, sample_data)
da_df <- calc_differential_abundance(tb, group = "group")
p <- plot_da_volcano(da_df)
```
