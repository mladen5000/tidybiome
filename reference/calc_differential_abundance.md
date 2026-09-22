# Multi-Engine Consensus Differential Abundance Analysis

Implements a state-of-the-art consensus differential abundance testing
framework inspired by recent bioRxiv preprints (*ConsensusMetaDA*,
*dar*, *LinDA*). Supports full formula specifications
(`formula = ~ treatment + age + batch`) to control for clinical
covariates, continuous biomarkers (e.g. `formula = ~ bmi`), and
multi-group experimental designs. Runs multiple complementary engines
(LinDA compositional regression, CLR linear models, two-part CAFT, and
rank/Wilcoxon tests) and computes consensus p-values (via the Cauchy
combination test), agreement scores, and effect sizes.

## Usage

``` r
calc_differential_abundance(
  tb,
  group = NULL,
  formula = NULL,
  covariates = NULL,
  contrast = NULL,
  methods = c("consensus", "caft", "linda", "clr_linear", "wilcoxon"),
  fdr_cutoff = 0.05,
  assay = "counts",
  pseudocount = 0.5
)
```

## Arguments

- tb:

  A `tidy_microbiome` object.

- group:

  Optional character string naming the primary predictor column in
  sample metadata.

- formula:

  Optional model formula specifying the target predictor and covariates
  (e.g. `~ treatment + age + batch`).

- covariates:

  Optional character vector of covariate names to adjust for if
  `formula` is omitted.

- contrast:

  Optional character string naming the specific predictor term of
  interest in `formula`. Defaults to the first term specified in
  `formula` or `group`.

- methods:

  Character vector of engines to run. Options:

  - `"consensus"`: Runs all available engines and aggregates them into a
    consensus verdict.

  - `"caft"`: Compositional Log-Linear Model with Zero Cells (bioRxiv
    Dec 2025), jointly modeling presence/absence probability and
    conditional log-abundance.

  - `"linda"`: Linear Models for Differential Abundance with
    compositional bias correction.

  - `"clr_linear"`: Linear regression on centered log-ratio abundances
    with covariate adjustment.

  - `"wilcoxon"`: Non-parametric two-sample Wilcoxon rank-sum or
    rank-regression test.

- fdr_cutoff:

  Numeric false discovery rate threshold (default: 0.05).

- assay:

  Character string naming the count assay to use. Defaults to
  `"counts"`.

- pseudocount:

  Numeric pseudocount for log transformations. Defaults to 0.5.

## Value

A
[tibble::tbl_df](https://tibble.tidyverse.org/reference/tbl_df-class.html)
with one row per taxon containing:

- `taxon_id`: Taxon identifier.

- `log2fc`: Average log2 fold change or regression slope across engines.

- `p_consensus`: Combined p-value using Cauchy combination test.

- `padj_consensus`: Benjamini-Hochberg adjusted consensus p-value.

- `agreement_score`: Number of engines agreeing on significance.

- `is_significant`: Logical flag whether taxon meets `fdr_cutoff`.

- Per-method statistics and taxonomic annotations.
