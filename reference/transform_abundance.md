# Transform Abundance Data Using Compositional or Classical Normalization

Applies modern compositional transformations (such as Robust CLR) or
classic normalizations (TSS, pseudocount CLR, presence/absence) to the
abundance assay of a `tidy_microbiome` object.

## Usage

``` r
transform_abundance(
  tb,
  method = c("rclr", "clr", "relabundance", "coverage", "hellinger", "gmpr", "log10",
    "pa"),
  assay = "counts",
  name = NULL,
  pseudocount = 1,
  min_overlap = 2
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- method:

  Character string specifying the transformation method. Options:

  - `"rclr"`: Robust Centered Log-Ratio (Martino et al. 2019). Zero
    values remain zero; positive values are centered by the geometric
    mean of observed (non-zero) taxa.

  - `"clr"`: Classic Centered Log-Ratio with pseudocount.

  - `"relabundance"`: Total Sum Scaling (TSS), scaling each sample to
    sum to 1.

  - `"coverage"`: Coverage-based standardization using Good-Turing
    sample coverage.

  - `"hellinger"`: Hellinger transformation (square root of relative
    abundance).

  - `"gmpr"`: Geometric Mean of Pairwise Ratios (Chen et al. 2018),
    zero-tolerant size-factor normalization.

  - `"log10"`: Log10 transformation with pseudocount.

  - `"pa"`: Presence/Absence binary transformation (0 or 1).

- assay:

  Character string naming the source assay to transform. Defaults to
  `"counts"`.

- name:

  Optional character string naming the output assay. Defaults to
  `method`.

- pseudocount:

  Numeric pseudocount added for `"clr"` and `"log10"`. Defaults to 1.

- min_overlap:

  Minimum number of shared non-zero taxa required for GMPR pairwise
  ratios (default: 2).

## Value

An updated `tidy_microbiome` object with the new assay added.

## Examples

``` r
counts <- matrix(c(100, 0, 50, 0, 20, 10, 0, 5), nrow = 2, byrow = TRUE,
                 dimnames = list(c("ASV1", "ASV2"), c("S1", "S2", "S3", "S4")))
tb <- tidy_microbiome(counts)
tb <- transform_abundance(tb, method = "rclr")
tb <- transform_abundance(tb, method = "hellinger")
```
