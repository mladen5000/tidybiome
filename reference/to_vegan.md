# Convert tidy_microbiome to vegan Community and Environmental Format

Convert tidy_microbiome to vegan Community and Environmental Format

## Usage

``` r
to_vegan(tb, assay = "counts")
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- assay:

  Name of assay to extract. Defaults to `"counts"`.

## Value

An S3 object of class `tidybiome_vegan` containing:

- `comm`: Community matrix with samples in rows and taxa in columns.

- `env`: Sample metadata tibble.
