# Extract Distance Matrix from tidy_microbiome

Extract Distance Matrix from tidy_microbiome

## Usage

``` r
get_distance(tb, metric = NULL)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- metric:

  Optional character string naming the distance metric to extract. If
  `NULL`, returns the most recently calculated distance matrix.

## Value

An object of class `dist`.
