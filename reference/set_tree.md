# Attach and Validate Phylogenetic Tree to tidy_microbiome

Attach and Validate Phylogenetic Tree to tidy_microbiome

## Usage

``` r
set_tree(tb, tree, prune = TRUE)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- tree:

  An [`ape::phylo`](https://rdrr.io/pkg/ape/man/read.tree.html) object.

- prune:

  Logical; if `TRUE`, automatically prunes the tree to match taxon IDs
  present in `tb`. Defaults to `TRUE`.

## Value

A `tidy_microbiome` object with the attached tree.
