# Import QIIME 2 Artifacts or TSV Tables

Constructs a `tidy_microbiome` container from QIIME 2 exported TSV files
(feature table, taxonomy TSV with semicolon delimiters, sample metadata,
and optional tree).

## Usage

``` r
import_qiime2(
  feature_table,
  taxonomy = NULL,
  sample_metadata = NULL,
  tree = NULL
)
```

## Arguments

- feature_table:

  Matrix, data frame, or file path to QIIME 2 feature table (TSV).

- taxonomy:

  Optional file path or data frame to QIIME 2 taxonomy TSV
  (`Feature.ID`, `Taxon`, `Confidence`).

- sample_metadata:

  Optional file path or data frame containing sample metadata.

- tree:

  Optional file path to a Newick tree file (`tree.nwk`) or an
  [`ape::phylo`](https://rdrr.io/pkg/ape/man/read.tree.html) object.

## Value

A `tidy_microbiome` object.
