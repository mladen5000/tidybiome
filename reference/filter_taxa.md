# Filter Taxa Based on Taxonomy Table Attributes

Filter Taxa Based on Taxonomy Table Attributes

## Usage

``` r
filter_taxa(tb, ...)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- ...:

  Logical predicates passed to
  [dplyr::filter](https://dplyr.tidyverse.org/reference/filter.html) on
  the taxonomy table.

## Value

A filtered `tidy_microbiome` object with synchronized assays, taxonomy,
and tree.
