# Calculate Dimensionality Reduction and Ordination

Performs ordination on microbiome data, including modern compositional
Robust PCA (RPCA, generating feature loadings for true compositional
biplots without distortion) and classical Principal Coordinates Analysis
(PCoA with negative eigenvalue correction) or standard PCA.

## Usage

``` r
calc_ordination(
  tb,
  method = c("rpca", "pcoa", "cpca", "pca"),
  metric = "raitchison",
  contrast_group = NULL,
  contrast_alpha = 1,
  assay = "counts",
  n_components = 3
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- method:

  Character string specifying the ordination method:

  - `"rpca"`: Robust PCA (DEICODE-style SVD on the `rclr` matrix).
    Computes both sample coordinates and taxon loading vectors for clean
    biplots.

  - `"pcoa"`: Principal Coordinates Analysis with automatic Lingoes
    correction for negative eigenvalues.

  - `"cpca"`: Contrastive PCA (cPCA). Isolates axes of microbial
    variation that are enriched in a target condition (e.g.
    disease/treatment) relative to a background (control), eliminating
    uninformative common background variability.

  - `"pca"`: Standard Principal Component Analysis on centered
    abundance.

- metric:

  Distance metric to use if `method = "pcoa"`. Defaults to
  `"raitchison"`.

- contrast_group:

  Optional character string naming the binary metadata column if
  `method = "cpca"`.

- contrast_alpha:

  Numeric contrast trade-off parameter for cPCA (default: 1.0).

- assay:

  Character string naming the abundance assay to use. Defaults to
  `"counts"`.

- n_components:

  Number of ordination axes to retain. Defaults to 3.

## Value

An updated `tidy_microbiome` object with ordination results stored in
its metadata.
