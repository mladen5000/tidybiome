# Calculate Beta Diversity Dissimilarity and Distance Matrices

Computes beta diversity distance matrices between samples, including
modern compositionally coherent distances (Robust Aitchison,
Jensen-Shannon) and classic ecological metrics (Bray-Curtis, Jaccard,
classic Aitchison). The resulting distance matrix is stored in the
`tidy_microbiome` metadata and can be retrieved using
[`get_distance()`](https://tidybiome.org/reference/get_distance.md).

## Usage

``` r
calc_beta_diversity(
  tb,
  metric = c("raitchison", "unifrac", "wunifrac", "wasserstein", "aitchison", "bray",
    "jaccard", "jsd"),
  assay = "counts",
  pseudocount = 1
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- metric:

  Character string specifying the distance metric:

  - `"raitchison"`: Robust Aitchison distance (Euclidean distance on the
    `rclr` matrix). Recommended for zero-inflated microbiome data as it
    avoids pseudocount distortion.

  - `"wasserstein"`: Optimal Transport Tree-Wasserstein / Earth Mover's
    Distance across taxonomic hierarchy. Measures physical work needed
    to move microbial mass across lineages.

  - `"aitchison"`: Classic Aitchison distance (Euclidean distance on
    pseudocount CLR).

  - `"bray"`: Bray-Curtis dissimilarity.

  - `"jaccard"`: Jaccard binary dissimilarity.

  - `"jsd"`: Jensen-Shannon distance (\\\sqrt{\text{JSD}}\\).

- assay:

  Character string naming the abundance assay to use. Defaults to
  `"counts"`.

- pseudocount:

  Numeric pseudocount for `"aitchison"`. Defaults to 1.

## Value

An updated `tidy_microbiome` object with the computed distance stored in
its metadata.

## Examples

``` r
counts <- matrix(c(10, 0, 5, 20, 10, 0, 5, 5, 5), nrow = 3, ncol = 3,
                 dimnames = list(c("A", "B", "C"), c("S1", "S2", "S3")))
tb <- tidy_microbiome(counts)
tb <- calc_beta_diversity(tb, metric = "raitchison")
d <- get_distance(tb, "raitchison")
```
