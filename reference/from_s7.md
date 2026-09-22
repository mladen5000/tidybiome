# Convert S7 Formal Object to tidy_microbiome S3 Tibble

Coerces an S7 `tidy_microbiome_s7` object back into a standard S3
`tidy_microbiome` tibble.

## Usage

``` r
from_s7(s7_obj)
```

## Arguments

- s7_obj:

  An object created by
  [`to_s7()`](https://tidybiome.org/reference/to_s7.md).

## Value

A `tidy_microbiome` S3 object.

## Examples

``` r
if (requireNamespace("S7", quietly = TRUE)) {
  data(gut_microbiome)
  tb_s7 <- to_s7(gut_microbiome)
  tb_restored <- from_s7(tb_s7)
  class(tb_restored)
}
#> [1] "tidy_microbiome" "tbl_df"          "tbl"             "data.frame"     
```
