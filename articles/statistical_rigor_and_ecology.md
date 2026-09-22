# Statistical Rigor and Ecological Modeling in tidybiome

## 1. Overview

Microbiome data present profound statistical challenges: 1. **Extreme
Sparsity**: High percentages of true and structural zero counts. 2.
**Compositionality**: Sequencing yields relative proportions constrained
to the simplex, where shifts in one taxon artificially inflate or
deflate others. 3. **High Dimensionality**: Hundreds to thousands of
taxa with small sample cohorts ($`P \gg N`$). 4. **Confounding
Covariates**: Age, sex, antibiotic history, batch effects, and
sequencing depth bias naive comparisons.

This vignette explains how **`tidybiome`** addresses these challenges
using modern statistical methods, formula-based linear modeling,
phylogenetic distance vectorization, and ecological modeling.

------------------------------------------------------------------------

## 2. Zero-Safe Compositional Transformations (Robust CLR)

### 2.1 The Problem with Pseudocounts

The classical Centered Log-Ratio (CLR) transform is defined as:
``` math
\text{clr}(x_i) = \left[ \log\frac{x_{i1}}{g(x_i)}, \dots, \log\frac{x_{iP}}{g(x_i)} \right]
```
where $`g(x_i) = (\prod_{j=1}^P x_{ij})^{1/P}`$ is the geometric mean
across all taxa.

Because microbiome data contain many zeros, analysts traditionally add
an arbitrary pseudocount: $`\log(x_{ij} + c)`$ for $`c \in \{0.5, 1\}`$.
However, as demonstrated by Martino et al. (2019), the choice of $`c`$
introduces severe mathematical artifacts: - Artificially alters sample
geometric means. - Biases sample distance calculations and ordination
projections. - Generates spurious differential abundance biomarkers.

### 2.2 The Robust CLR (rCLR) Solution

`tidybiome` implements **rCLR**, which calculates the geometric mean
strictly over observed non-zero counts:
``` math
g_{>0}(x_i) = \left( \prod_{j: x_{ij} > 0} x_{ij} \right)^{1 / |\{j: x_{ij} > 0\}|}
```
For non-zero taxa:
``` math
\text{rclr}(x_{ij}) = \log\frac{x_{ij}}{g_{>0}(x_i)}
```
For zero counts ($`x_{ij} = 0`$), $`\text{rclr}(x_{ij})`$ remains 0
without any imputation or distortion:

``` r

library(tidybiome)
data("gut_microbiome", package = "tidybiome")

tb_rclr <- transform_abundance(gut_microbiome, method = "rclr")
# Zero values in raw counts remain exactly zero in rclr
raw_counts <- assay(tb_rclr, "counts")
rclr_counts <- assay(tb_rclr, "rclr")

all(rclr_counts[raw_counts == 0] == 0)
#> [1] TRUE
```

------------------------------------------------------------------------

## 3. Formula-Based Differential Abundance with Covariates

### 3.1 Confounder Control with Formula Syntax

In clinical and ecological surveys, naive two-group tests (like Wilcoxon
or two-sample $`t`$-tests) fail to account for confounding factors like
age, BMI, or sequencing batch.

`tidybiome`’s
[`calc_differential_abundance()`](https://tidybiome.org/reference/calc_differential_abundance.md)
accepts standard R formula syntax:

``` r

calc_differential_abundance(tb, formula = ~ treatment + age + batch, contrast = "treatment")
```

### 3.2 Continuous Predictors

In addition to categorical treatments, `tidybiome` supports continuous
response modeling (e.g. testing for taxa that correlate with continuous
age or biomarker concentration):

``` r

da_continuous <- calc_differential_abundance(
  gut_microbiome,
  formula = ~ age,
  methods = c("clr_linear", "linda")
)

da_continuous |>
  head(5)
#> # A tibble: 5 × 15
#>   taxon_id log2fc_linda  p_linda padj_linda log2fc_clr   p_clr padj_clr
#>   <chr>           <dbl>    <dbl>      <dbl>      <dbl>   <dbl>    <dbl>
#> 1 ASV16          0.127  0.00105      0.0210     0.121  0.00138   0.0509
#> 2 ASV03         -0.0542 0.000705     0.0210    -0.0618 0.00255   0.0509
#> 3 ASV01         -0.0899 0.0459       0.612     -0.0966 0.0299    0.399 
#> 4 ASV13          0.102  0.0616       0.616      0.0948 0.0715    0.661 
#> 5 ASV04          0.0952 0.104        0.682      0.0866 0.128     0.666 
#> # ℹ 8 more variables: p_consensus <dbl>, padj_consensus <dbl>, log2fc <dbl>,
#> #   agreement_score <int>, is_significant <lgl>, Kingdom <chr>, Phylum <chr>,
#> #   Genus <chr>
```

### 3.3 Engine Architectures

1.  **CAFT (bioRxiv Dec 2025)**: A two-part likelihood ratio model
    designed for zero-inflated microbiome data. Separates:
    - Zero-cell dropout probability via logistic regression:
      $`\text{logit}(P(Y=0)) \sim \mathbf{X}\boldsymbol{\beta}`$
    - Abundance shifts among present counts via log-linear regression:
      $`\log(Y \mid Y > 0) \sim \mathbf{X}\boldsymbol{\gamma}`$ Combines
      likelihood ratio test statistics while gracefully suppressing
      separation warnings.
2.  **LinDA (Zhou et al., 2022)**: Compositional linear regression
    selecting a robust reference taxon via mode/rank estimation, fully
    controlling for model matrix covariates.
3.  **CLR-Linear**: Fits multiple linear regression models on rCLR or
    CLR transformed features.
4.  **Cauchy Combination Aggregation**: Combines $`p`$-values across
    disparate statistical frameworks into a single unified consensus
    $`p`$-value (`padj_consensus`), weighted by engine agreement score.

------------------------------------------------------------------------

## 4. Vectorized Phylogenetic Ecology (UniFrac)

### 4.1 Vectorized Descendant-Edge Incidence

Calculating UniFrac traditionally requires deep recursive traversals
over the tree for every pair of samples ($`O(N^2 \cdot \text{nodes})`$),
leading many R packages to rely on external C/C++ libraries.

`tidybiome` implements a vectorized algorithm using an **edge-descendant
incidence matrix** $`\mathbf{M} \in \{0, 1\}^{E \times T}`$, where $`E`$
is the number of edges in the tree and $`T`$ is the number of tips. Each
tip is traced to the root in a single precomputation step
($`O(T \cdot \text{depth})`$). Multiplying:
``` math
\mathbf{W} = \mathbf{M} \times \mathbf{P}
```
where $`\mathbf{P} \in [0, 1]^{T \times N}`$ is the sample proportion
matrix, computes the total descendant branch weight for all edges across
all $`N`$ samples simultaneously in pure R matrix algebra.

### 4.2 Unweighted & Weighted UniFrac

- **Unweighted UniFrac**:
  ``` math
  d_{\text{unifrac}}(A, B) = \frac{\sum_e l_e |\mathbb{I}(W_{eA} > 0) - \mathbb{I}(W_{eB} > 0)|}{\sum_e l_e \mathbb{I}(W_{eA} + W_{eB} > 0)}
  ```
- **Weighted Normalized UniFrac**:
  ``` math
  d_{\text{wunifrac}}(A, B) = \frac{\sum_e l_e |W_{eA} - W_{eB}|}{\sum_e l_e (W_{eA} + W_{eB})}
  ```

``` r

taxa_ids <- attr(gut_microbiome, "tax_table")$taxon_id
set.seed(123)
tree <- ape::rtree(length(taxa_ids), tip.label = taxa_ids)
tb_phylo <- set_tree(gut_microbiome, tree)

# Fast vectorized calculations:
dist_unweighted <- calc_unifrac(tb_phylo, weighted = FALSE)
dist_weighted   <- calc_unifrac(tb_phylo, weighted = TRUE, normalized = TRUE)

as.matrix(dist_unweighted)[1:3, 1:3]
#>            Sample_01  Sample_02  Sample_03
#> Sample_01 0.00000000 0.09757923 0.03646537
#> Sample_02 0.09757923 0.00000000 0.06182808
#> Sample_03 0.03646537 0.06182808 0.00000000
as.matrix(dist_weighted)[1:3, 1:3]
#>            Sample_01  Sample_02  Sample_03
#> Sample_01 0.00000000 0.06479822 0.01380644
#> Sample_02 0.06479822 0.00000000 0.06485739
#> Sample_03 0.01380644 0.06485739 0.00000000
```

------------------------------------------------------------------------

## 5. Distance-Based Redundancy Analysis (db-RDA)

Distance-based Redundancy Analysis (Legendre & Anderson, 1999) is the
gold standard for constrained ordination of non-Euclidean ecological
distance matrices (Bray-Curtis, Jaccard, UniFrac).

[`tidybiome::run_dbrda()`](https://tidybiome.org/reference/run_dbrda.md)
fits the constrained linear ordination:
``` math
\mathbf{Y}_{\text{dist}} \sim \mathbf{X}_{\text{metadata}}
```

``` r

dbrda_fit <- run_dbrda(gut_microbiome, formula = ~ treatment + age, distance = "bray")
dbrda_fit
#> -- tidybiome_dbrda (Distance-Based Redundancy Analysis) --
#>   * Constrained axes: 2
#>   * Inertia explained: dbRDA1 = 77.84%, dbRDA2 = 1.01%
#>   * Samples: 24
```

The output returns: - `samples`: Ordination coordinates for all
samples. - `biplot`: Directional arrows for continuous and categorical
predictors. - `variance_explained`: Proportion of constrained inertia
explained by each canonical axis. - `anova`: Permutation ANOVA table
testing overall model significance.

------------------------------------------------------------------------

## 6. Microbial Co-Occurrence Networks

Pairwise co-occurrence networks identify putative cooperative guilds,
mutual exclusions, and ecological keystones:

``` r

net <- calc_network(
  gut_microbiome,
  method = "spearman",
  min_prevalence = 0.20,
  r_cutoff = 0.35,
  p_cutoff = 0.05
)

# Nodes table with degree centrality
head(net$nodes, 5)
#> # A tibble: 5 × 6
#>   taxon_id degree mean_abundance Kingdom  Phylum       Genus      
#>   <chr>     <int>          <dbl> <chr>    <chr>        <chr>      
#> 1 ASV01        29        0.162   Bacteria Bacteroidota Bacteroides
#> 2 ASV02        30        0.0347  Bacteria Bacteroidota Bacteroides
#> 3 ASV03        30        0.0916  Bacteria Bacteroidota Bacteroides
#> 4 ASV04        29        0.0172  Bacteria Bacteroidota Bacteroides
#> 5 ASV05        21        0.00361 Bacteria Bacteroidota Bacteroides

# Edges table with correlation strength and direction
head(net$edges, 5)
#> # A tibble: 5 × 7
#>   from  to    correlation       p_value        padj weight direction
#>   <chr> <chr>       <dbl>         <dbl>       <dbl>  <dbl> <chr>    
#> 1 ASV01 ASV02       0.876 0.0000000211  0.000000548  0.876 positive 
#> 2 ASV01 ASV03       0.887 0.00000000780 0.000000338  0.887 positive 
#> 3 ASV02 ASV03       0.897 0.00000000309 0.000000193  0.897 positive 
#> 4 ASV01 ASV04      -0.552 0.00516       0.0107       0.552 negative 
#> 5 ASV02 ASV04      -0.465 0.0221        0.0371       0.465 negative
```

The returned `tidybiome_network` object contains tidy data frames that
can be directly piped into `plot_network(net)` or external graph
packages like `igraph` or `tidygraph`.

------------------------------------------------------------------------

## 7. Container Integrity & Diagnostic Validation

To ensure analytical reproducibility, `tidybiome` provides
[`validate_tidy_microbiome()`](https://tidybiome.org/reference/validate_tidy_microbiome.md):

``` r

validate_tidy_microbiome(gut_microbiome, verbose = TRUE)
#> Container validation passed: all assays, taxonomy, and tree tips synchronized.
```

This verifies: 1. **Sample ID Synchronization**: Ensures row names of
metadata match column names of all assays. 2. **Taxon ID
Synchronization**: Ensures row names of taxonomy match row names of all
assays. 3. **Phylogenetic Tip Alignment**: If a tree is present,
verifies tree tip labels correspond 1:1 with taxon IDs.

------------------------------------------------------------------------

## References

1.  **Martino, C. et al. (2019)**. *A Novel Sparse Compositional Metric
    Alleviates Problems in Biomarker Discovery*. **mSystems**, 4(1),
    e00016-19.
2.  **Zhou, H. et al. (2022)**. *LinDA: linear models for differential
    abundance analysis of microbiome compositional data*. **Genome
    Biology**, 23, 95.
3.  **Lozupone, C. & Knight, R. (2005)**. *UniFrac: a new phylogenetic
    metric for comparing microbial communities*. **Applied and
    Environmental Microbiology**, 71(12), 8228-8235.
4.  **Legendre, P. & Anderson, M. J. (1999)**. *Distance-based
    redundancy analysis: testing multispecies responses in
    multifactorial ecological experiments*. **Ecological Monographs**,
    69(1), 1-24.
