# Calculate Microbial Co-Occurrence Network

Evaluates pairwise co-occurrence correlations across taxa and constructs
network nodes and edges filtered by correlation strength and statistical
significance.

## Usage

``` r
calc_network(
  tb,
  assay = "counts",
  method = c("spearman", "pearson"),
  min_prevalence = 0.2,
  r_cutoff = 0.4,
  p_cutoff = 0.05,
  p_adj_method = "BH"
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- assay:

  Assay to evaluate. Defaults to `"counts"`.

- method:

  Correlation method: `"spearman"` or `"pearson"`. Defaults to
  `"spearman"`.

- min_prevalence:

  Minimum proportion of samples where taxon must be detected (default:
  0.20).

- r_cutoff:

  Minimum absolute correlation coefficient threshold (default: 0.40).

- p_cutoff:

  Maximum multiple testing adjusted p-value threshold (default: 0.05).

- p_adj_method:

  Multiple testing adjustment method (default: `"BH"`).

## Value

An S3 object of class `tidybiome_network` containing:

- `nodes`: Tibble of network nodes with `taxon_id`, `degree`,
  `mean_abundance`, and taxonomy.

- `edges`: Tibble of network edges with `from`, `to`, `correlation`,
  `p_value`, `padj`, `weight`, and `direction`.

- `r_cutoff`: Filtering cutoff for r.

- `p_cutoff`: Filtering cutoff for p.
