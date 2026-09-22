# Calculate Microbial Co-Occurrence Network

Evaluates pairwise co-occurrence correlations across taxa and constructs
a network representation including edge lists, node degree, and
community structures.

## Usage

``` r
calc_network(
  tb,
  assay = "counts",
  method = "spearman",
  min_prevalence = 0.2,
  r_cutoff = 0.3,
  p_cutoff = 0.05,
  p_adj_method = "BH"
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- assay:

  Name of assay to compute correlations on (default: `"counts"`).

- method:

  Correlation method: `"spearman"` (default), `"pearson"`, or
  `"kendall"`.

- min_prevalence:

  Minimum taxon prevalence threshold (0-1) to filter before correlation.

- r_cutoff:

  Absolute correlation magnitude threshold for retaining edges (default:
  0.3).

- p_cutoff:

  Adjusted p-value threshold for edge significance (default: 0.05).

- p_adj_method:

  Multiple testing correction method (default: `"BH"`).

## Value

A `tidybiome_network` object containing:

- nodes:

  Tibble of taxa, taxonomy, and node degree.

- edges:

  Tibble of significant pairwise edges, correlation, and p-values.

- adjacency:

  Filtered adjacency matrix.

- params:

  List of parameter settings.

## Examples

``` r
data(gut_microbiome)
net <- calc_network(gut_microbiome, min_prevalence = 0.5, r_cutoff = 0.3)
print(net)
#> -- tidybiome_network (Microbial Co-Occurrence Network) --
#>   * Nodes (Taxa): 38
#>   * Significant edges: 491 (|r| >= 0.30, padj <= 0.05)
#>   * Positive: 404, Negative: 87
```
