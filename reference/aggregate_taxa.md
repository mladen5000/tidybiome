# Aggregate Taxa to a Higher Taxonomic Rank

Agglomerates/merges features (ASVs/OTUs) by a specified taxonomic rank
(e.g., Phylum, Family, Genus) by summing abundances across all assays.

## Usage

``` r
aggregate_taxa(tb, rank, na.rm = FALSE)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- rank:

  Character string specifying the taxonomic column to aggregate by.

- na.rm:

  Logical. If `TRUE`, removes unassigned taxa. If `FALSE` (default),
  groups them as `"Unclassified"`.

## Value

A new `tidy_microbiome` object aggregated at the specified rank.
