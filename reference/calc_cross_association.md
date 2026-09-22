# Calculate Cross-Associations Between Taxa and Metadata Variables

Computes pairwise associations (Spearman or Pearson correlations and
significance) between microbial feature abundances and continuous
metadata covariates (e.g. host age, BMI, clinical metrics).

## Usage

``` r
calc_cross_association(
  tb,
  variables = NULL,
  assay = "counts",
  method = c("spearman", "pearson"),
  p_adj_method = "BH",
  sort = TRUE
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- variables:

  Character vector of numeric metadata columns in `tb` to correlate with
  taxa. If `NULL`, automatically selects all numeric columns in sample
  metadata. Defaults to `NULL`.

- assay:

  Name of assay to evaluate. Defaults to `"counts"`.

- method:

  Correlation method: `"spearman"` or `"pearson"`. Defaults to
  `"spearman"`.

- p_adj_method:

  Multiple testing correction method passed to
  [`stats::p.adjust`](https://rdrr.io/r/stats/p.adjust.html) (e.g.
  `"BH"`, `"bonferroni"`). Defaults to `"BH"`.

- sort:

  Logical; if `TRUE`, sorts output by adjusted p-value (`padj`).
  Defaults to `TRUE`.

## Value

A tidy tibble with columns:

- `taxon_id`: Taxon identifier.

- `variable`: Metadata covariate name.

- `correlation`: Estimated correlation coefficient.

- `p_value`: Raw association p-value from `cor.test`.

- `padj`: Multiple-testing adjusted p-value.

Joined with taxonomic lineages from `tax_table`.
