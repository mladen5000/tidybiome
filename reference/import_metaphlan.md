# Import MetaPhlAn Taxonomic Profiles

Parses MetaPhlAn (version 3 or 4) merged abundance profile tables
containing hierarchical clade strings (`k__...|p__...|s__...`) and
relative abundances.

## Usage

``` r
import_metaphlan(file, rank = "Species", sample_metadata = NULL)
```

## Arguments

- file:

  Path to MetaPhlAn merged abundance TSV file.

- rank:

  Taxonomic rank to extract (e.g. `"Species"`, `"Genus"`, or `"all"`).
  Defaults to `"Species"`.

- sample_metadata:

  Optional data frame or TSV file path containing sample metadata.

## Value

A `tidy_microbiome` object with relative abundance stored in assay
`"relabundance"` and `"counts"`.
