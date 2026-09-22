# Extract Assay Matrix from tidy_microbiome

Extract Assay Matrix from tidy_microbiome

## Usage

``` r
assay(x, name = "counts", ...)
```

## Arguments

- x:

  A `tidy_microbiome` object.

- name:

  Character string naming the assay to extract (e.g. `"counts"`,
  `"rclr"`, `"tss"`). Defaults to `"counts"`.

- ...:

  Additional arguments (not used).

## Value

A numeric matrix with taxa as rows and samples as columns.
