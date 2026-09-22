# Calculate Alpha Diversity and Simplex Compositional Metrics

Calculates alpha diversity indices including the unified Hill numbers
profile (\\q = 0, 1, 2\\), simplex compositional variation (bioRxiv
2026), and classic ecological indices (Shannon, Gini-Simpson, Inverse
Simpson, Observed, Chao1, and Faith's PD).

## Usage

``` r
calc_alpha_diversity(
  tb,
  metrics = c("hill", "simplex_variation", "shannon", "simpson", "inv_simpson",
    "observed", "chao1"),
  assay = "counts"
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- metrics:

  Character vector of metrics to calculate. Supported options:

  - `"hill"`: Calculates the Hill numbers profile: `hill_0` (species
    richness), `hill_1` (exponential of Shannon entropy), and `hill_2`
    (inverse Simpson index).

  - `"simplex_variation"`: Compositional total variance on the closed
    simplex (bioRxiv 2026), measuring true geometric dispersion among
    taxon ratios.

  - `"shannon"`: Shannon diversity index (\\H = -\sum p_i \ln p_i\\).

  - `"simpson"`: Gini-Simpson index (\\1 - \sum p_i^2\\).

  - `"inv_simpson"`: Inverse Simpson index (\\1 / \sum p_i^2\\).

  - `"observed"`: Number of observed taxa with counts \> 0.

  - `"chao1"`: Chao1 bias-corrected richness estimator.

  - `"faith_pd"`: Faith's phylogenetic diversity (requires `phy_tree`).

- assay:

  Character string naming the abundance assay to use. Defaults to
  `"counts"`.

## Value

An updated `tidy_microbiome` object with new alpha diversity columns
added.
