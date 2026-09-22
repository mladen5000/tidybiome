# Calculate GMPR Size Factors

Calculates sample size factors using the Geometric Mean of Pairwise
Ratios (GMPR) methodology (Chen et al. 2018), specifically tailored to
zero-inflated microbiome counts.

## Usage

``` r
calc_gmpr_size_factors(mat, min_overlap = 2)
```

## Arguments

- mat:

  A numeric matrix of counts (taxa as rows, samples as columns).

- min_overlap:

  Minimum number of shared non-zero taxa between sample pairs (default:
  2).

## Value

A named numeric vector of size factors.
