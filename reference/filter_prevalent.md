# Filter Microbiome by Prevalence and Abundance Thresholds

Keeps only taxa exceeding the specified minimum prevalence (fraction of
samples) and mean relative abundance, keeping count matrices and
taxonomy in sync.

## Usage

``` r
filter_prevalent(
  tb,
  min_prevalence = 0.1,
  min_abundance = 0,
  assay = "counts",
  detection = 0
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- min_prevalence:

  Minimum proportion of samples where taxon must be detected (0 to 1).
  Defaults to `0.10`.

- min_abundance:

  Minimum mean abundance across samples. Defaults to `0`.

- assay:

  Assay to evaluate. Defaults to `"counts"`.

- detection:

  Abundance threshold for presence. Defaults to `0`.

## Value

A filtered `tidy_microbiome` object.
