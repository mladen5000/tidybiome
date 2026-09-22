# Plot Alpha Diversity with Statistical Comparisons

Generates publication-grade hybrid raincloud/boxplot visuals for alpha
diversity metrics (including Hill numbers \\q = 0, 1, 2\\), with
overlaid sample points and automated two-group statistical significance
annotations.

## Usage

``` r
plot_alpha(
  tb,
  metric = "hill_1",
  x = NULL,
  color = NULL,
  test = TRUE,
  palette = "tidybiome"
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- metric:

  Character string naming the alpha diversity metric (e.g. `"hill_1"`,
  `"hill_0"`, `"shannon"`).

- x:

  Character string naming the metadata column for the horizontal
  grouping axis.

- color:

  Optional character string naming the metadata column to color by.

- test:

  Logical. If `TRUE` and `x` has 2 groups, runs a Wilcoxon rank-sum test
  and displays the p-value.

- palette:

  Palette name: `"tidybiome"` or `"nature"`.

## Value

A
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Examples

``` r
counts <- matrix(c(25, 25, 25, 25, 10, 5, 2, 80), nrow = 4, ncol = 2,
                 dimnames = list(paste0("T", 1:4), c("S1", "S2")))
sample_data <- data.frame(sample_id = c("S1", "S2"), treatment = c("Control", "Treated"))
tb <- tidy_microbiome(counts, sample_data)
tb <- calc_alpha_diversity(tb, metrics = "hill")
p <- plot_alpha(tb, metric = "hill_1", x = "treatment")
```
