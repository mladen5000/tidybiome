# Getting Started with tidybiome: A Modern, Tidyverse-Native Microbiome Suite

## 1. Introduction

Downstream microbiome analysis has historically been divided between
legacy packages like `phyloseq` and complex Bioconductor frameworks like
`mia` / `TreeSummarizedExperiment`. While powerful, these tools often
require navigating intricate S4 slots, rely on outdated heuristics (such
as arbitrary rarefaction or pseudocount-based CLR), and output visuals
with hard-to-read 50-color legends.

**`tidybiome`** brings microbiome data science directly into the modern
**tidyverse** ecosystem. It acts as a lightweight, aesthetic Swiss Army
knife built around the following core principles:

1.  **Native Tidyverse Ergonomics**: The `tidy_microbiome` container
    inherits from `tbl_df`. Standard `dplyr` verbs (`filter`, `mutate`,
    `select`, `arrange`, `slice`) work seamlessly on sample metadata
    while keeping count matrices, taxonomic lineages, distance caches,
    and phylogenetic trees synchronized under the hood.
2.  **State-of-the-Art & bioRxiv Methodologies**:
    - **Robust CLR (rCLR)**: Zero-safe compositional normalization
      without pseudocount distortion (Martino et al., 2019).
    - **Unified Hill Numbers Profile**: Richness ($`q=0`$), exponential
      Shannon ($`q=1`$), and inverse Simpson ($`q=2`$) on the intuitive
      scale of effective number of species.
    - **Phylogenetic UniFrac**: High-performance vectorized unweighted
      and weighted UniFrac distances via ancestor-descendant matrix
      algebra.
    - **Consensus Differential Abundance**: Robust multi-engine
      consensus (LinDA, CAFT zero-cell separation, CLR-linear models,
      Wilcoxon) with formula specification
      (`formula = ~ treatment + age + batch`), covariate adjustments,
      and Cauchy combination testing.
    - **Constrained Ordination & Networks**: Distance-based Redundancy
      Analysis (db-RDA) and microbial co-occurrence networks with false
      discovery rate (FDR) control and degree centrality.
3.  **Publication-Grade Visualizations**: Soft-gray `"Other"` pooling,
    colorblind-safe palettes (Okabe-Ito, Nature), confidence ellipses,
    and
    [`theme_tidybiome()`](https://tidybiome.org/reference/theme_tidybiome.md).

------------------------------------------------------------------------

## 2. Container Setup and Inspection

Let’s load `tidybiome`, `dplyr`, and `ggplot2`, and inspect the built-in
benchmark dataset `gut_microbiome` (24 samples across 40 taxa):

``` r

library(tidybiome)
library(dplyr)
library(ggplot2)

data("gut_microbiome", package = "tidybiome")
gut_microbiome
#> ── tidy_microbiome [24 samples × 40 taxa] ──────────── 
#>   • Assays:         counts 
#>   • Taxonomy ranks: Kingdom > Phylum > Genus 
#>   • Tree:           None 
#> ── Sample Metadata ─────────────────────────────────── 
#> # A tibble: 24 × 5
#>    sample_id treatment diet        age depth
#>    <chr>     <fct>     <fct>     <dbl> <dbl>
#>  1 Sample_01 Control   Standard     54 22212
#>  2 Sample_02 Control   HighFiber    34 17287
#>  3 Sample_03 Control   Standard     44 18716
#>  4 Sample_04 Control   HighFiber    46 18906
#>  5 Sample_05 Control   Standard     44 20561
#>  6 Sample_06 Control   HighFiber    39 14289
#>  7 Sample_07 Control   Standard     55 17019
#>  8 Sample_08 Control   HighFiber    39 18048
#>  9 Sample_09 Control   Standard     60 19046
#> 10 Sample_10 Control   HighFiber    39 20230
#> # ℹ 14 more rows
```

Notice the container printout: it displays the sample and feature
dimensions, registered assays (`counts`), taxonomy ranks
(`Kingdom > Phylum > Genus`), and prints the sample metadata table
directly.

------------------------------------------------------------------------

## 3. Tidy Wrangling and Data Integrity

### 3.1 Using dplyr Verbs

Because `gut_microbiome` is a tibble, you can filter, arrange, and
mutate sample metadata directly:

``` r

# Filter to samples with sequencing depth > 15,000 and calculate log10-depth
gut_filtered <- gut_microbiome |>
  filter(depth > 15000) |>
  mutate(log_depth = log10(depth))

gut_filtered
#> ── tidy_microbiome [21 samples × 40 taxa] ──────────── 
#>   • Assays:         counts 
#>   • Taxonomy ranks: Kingdom > Phylum > Genus 
#>   • Tree:           None 
#> ── Sample Metadata ─────────────────────────────────── 
#> # A tibble: 21 × 6
#>    sample_id treatment diet        age depth log_depth
#>    <chr>     <fct>     <fct>     <dbl> <dbl>     <dbl>
#>  1 Sample_01 Control   Standard     54 22212      4.35
#>  2 Sample_02 Control   HighFiber    34 17287      4.24
#>  3 Sample_03 Control   Standard     44 18716      4.27
#>  4 Sample_04 Control   HighFiber    46 18906      4.28
#>  5 Sample_05 Control   Standard     44 20561      4.31
#>  6 Sample_07 Control   Standard     55 17019      4.23
#>  7 Sample_08 Control   HighFiber    39 18048      4.26
#>  8 Sample_09 Control   Standard     60 19046      4.28
#>  9 Sample_10 Control   HighFiber    39 20230      4.31
#> 10 Sample_11 Control   Standard     53 18876      4.28
#> # ℹ 11 more rows
```

All underlying assay matrices (`counts`) are automatically sliced along
the column dimension to match the remaining samples!

### 3.2 S3 Bracket Subsetting (`[`, `[<-`)

Direct matrix-style bracket indexing is fully supported:

``` r

# Subset the first 6 samples and select metadata columns
gut_slice <- gut_microbiome[1:6, c("sample_id", "treatment", "age")]
gut_slice
#> ── tidy_microbiome [6 samples × 40 taxa] ───────────── 
#>   • Assays:         counts 
#>   • Taxonomy ranks: Kingdom > Phylum > Genus 
#>   • Tree:           None 
#> ── Sample Metadata ─────────────────────────────────── 
#> # A tibble: 6 × 3
#>   sample_id treatment   age
#>   <chr>     <fct>     <dbl>
#> 1 Sample_01 Control      54
#> 2 Sample_02 Control      34
#> 3 Sample_03 Control      44
#> 4 Sample_04 Control      46
#> 5 Sample_05 Control      44
#> 6 Sample_06 Control      39
```

### 3.3 Taxonomic Filtering (`filter_taxa`)

To filter along the feature/taxonomic dimension, use
[`filter_taxa()`](https://tidybiome.org/reference/filter_taxa.md). This
slices taxonomic rows, assay matrix rows, and automatically prunes
phylogenetic tree tips:

``` r

# Filter to Bacteroidota taxa
gut_bacteroidota <- filter_taxa(gut_microbiome, Phylum == "Bacteroidota")
gut_bacteroidota
#> ── tidy_microbiome [24 samples × 15 taxa] ──────────── 
#>   • Assays:         counts 
#>   • Taxonomy ranks: Kingdom > Phylum > Genus 
#>   • Tree:           None 
#> ── Sample Metadata ─────────────────────────────────── 
#> # A tibble: 24 × 5
#>    sample_id treatment diet        age depth
#>    <chr>     <fct>     <fct>     <dbl> <dbl>
#>  1 Sample_01 Control   Standard     54 22212
#>  2 Sample_02 Control   HighFiber    34 17287
#>  3 Sample_03 Control   Standard     44 18716
#>  4 Sample_04 Control   HighFiber    46 18906
#>  5 Sample_05 Control   Standard     44 20561
#>  6 Sample_06 Control   HighFiber    39 14289
#>  7 Sample_07 Control   Standard     55 17019
#>  8 Sample_08 Control   HighFiber    39 18048
#>  9 Sample_09 Control   Standard     60 19046
#> 10 Sample_10 Control   HighFiber    39 20230
#> # ℹ 14 more rows
```

### 3.4 Container Validation

You can run automated integrity checks at any point:

``` r

validate_tidy_microbiome(gut_bacteroidota, verbose = TRUE)
#> Container validation passed: all assays, taxonomy, and tree tips synchronized.
```

------------------------------------------------------------------------

## 4. Modern Compositional Transformations (rCLR)

Microbiome sequencing data are high-dimensional, sparse compositions.
Adding arbitrary pseudocounts (e.g. $`+1`$) before logarithmic
transformation creates severe artifacts for rare taxa.

`tidybiome` implements **Robust CLR (rCLR)** (Martino et al., 2019),
which calculates geometric means restricted strictly to observed
non-zero counts:

``` r

gut_microbiome <- transform_abundance(gut_microbiome, method = "rclr")
gut_microbiome
#> ── tidy_microbiome [24 samples × 40 taxa] ──────────── 
#>   • Assays:         counts, rclr 
#>   • Taxonomy ranks: Kingdom > Phylum > Genus 
#>   • Tree:           None 
#> ── Sample Metadata ─────────────────────────────────── 
#> # A tibble: 24 × 5
#>    sample_id treatment diet        age depth
#>    <chr>     <fct>     <fct>     <dbl> <dbl>
#>  1 Sample_01 Control   Standard     54 22212
#>  2 Sample_02 Control   HighFiber    34 17287
#>  3 Sample_03 Control   Standard     44 18716
#>  4 Sample_04 Control   HighFiber    46 18906
#>  5 Sample_05 Control   Standard     44 20561
#>  6 Sample_06 Control   HighFiber    39 14289
#>  7 Sample_07 Control   Standard     55 17019
#>  8 Sample_08 Control   HighFiber    39 18048
#>  9 Sample_09 Control   Standard     60 19046
#> 10 Sample_10 Control   HighFiber    39 20230
#> # ℹ 14 more rows
```

The transformed assay is registered under `rclr` without modifying
original sequence counts.

------------------------------------------------------------------------

## 5. Alpha Diversity: Unified Hill Profiles

Classic alpha diversity metrics (Chao1, Shannon entropy, Gini-Simpson)
operate on disparate mathematical scales. `tidybiome` computes the
unified **Hill numbers profile**: \* $`q = 0`$: Species Richness (count
of observed taxa) \* $`q = 1`$: Exponential Shannon entropy (effective
number of typical species) \* $`q = 2`$: Inverse Simpson index
(effective number of dominant species)

``` r

gut_microbiome <- calc_alpha_diversity(
  gut_microbiome,
  metrics = c("hill", "shannon", "simplex_variation")
)

gut_microbiome |>
  select(sample_id, treatment, hill_0, hill_1, hill_2, simplex_variation) |>
  head(6)
#> ── tidy_microbiome [6 samples × 40 taxa] ───────────── 
#>   • Assays:         counts, rclr 
#>   • Taxonomy ranks: Kingdom > Phylum > Genus 
#>   • Tree:           None 
#> ── Sample Metadata ─────────────────────────────────── 
#> # A tibble: 6 × 6
#>   sample_id treatment hill_0 hill_1 hill_2 simplex_variation
#>   <chr>     <fct>      <dbl>  <dbl>  <dbl>             <dbl>
#> 1 Sample_01 Control       37   32.1   28.8             0.375
#> 2 Sample_02 Control       33   29.0   26.2             0.338
#> 3 Sample_03 Control       37   32.3   29.2             0.359
#> 4 Sample_04 Control       35   30.3   27.3             0.430
#> 5 Sample_05 Control       35   30.1   27.0             0.423
#> 6 Sample_06 Control       33   29.1   26.5             0.346
```

------------------------------------------------------------------------

## 6. Beta Diversity & Phylogenetic UniFrac

### 6.1 Compositional & Ecological Distances

[`calc_beta_diversity()`](https://tidybiome.org/reference/calc_beta_diversity.md)
supports: \* `"raitchison"`: Robust Aitchison distance (Euclidean on
rCLR-transformed counts). \* `"bray"`: Fast vectorized Bray-Curtis
dissimilarity. \* `"jaccard"`: Vectorized binary Jaccard incidence
dissimilarity. \* `"tree_wasserstein"`: Earth Mover’s Distance across
the taxonomic tree hierarchy.

``` r

gut_microbiome <- calc_beta_diversity(gut_microbiome, metric = "raitchison")
```

### 6.2 Phylogenetic UniFrac

When an [`ape::phylo`](https://rdrr.io/pkg/ape/man/read.tree.html) tree
is attached (or simulated),
[`calc_unifrac()`](https://tidybiome.org/reference/calc_unifrac.md)
computes unweighted and weighted normalized UniFrac distances via
vectorized edge-descendant matrix multiplication:

``` r

taxa_ids <- attr(gut_microbiome, "tax_table")$taxon_id
set.seed(42)
mock_tree <- ape::rtree(length(taxa_ids), tip.label = taxa_ids)
gut_phylo <- set_tree(gut_microbiome, mock_tree)

# Compute unweighted UniFrac
u_dist <- calc_unifrac(gut_phylo, weighted = FALSE)
as.matrix(u_dist)[1:4, 1:4]
#>            Sample_01  Sample_02  Sample_03  Sample_04
#> Sample_01 0.00000000 0.05980237 0.03599487 0.11929143
#> Sample_02 0.05980237 0.00000000 0.02418161 0.10897549
#> Sample_03 0.03599487 0.02418161 0.00000000 0.08522762
#> Sample_04 0.11929143 0.10897549 0.08522762 0.00000000
```

------------------------------------------------------------------------

## 7. Dimensionality Reduction (RPCA Biplot & cPCA)

### 7.1 Robust PCA (RPCA) Biplots

RPCA decomposes the rCLR assay using singular value decomposition,
providing true sample coordinates and taxon loading vectors:

``` r

gut_microbiome <- calc_ordination(gut_microbiome, method = "rpca")
```

------------------------------------------------------------------------

## 8. Statistical Rigor: Differential Abundance with Covariates

Microbial differential abundance analysis is notoriously susceptible to
false discoveries. `tidybiome` provides a comprehensive **Multi-Engine
Consensus Architecture** supporting formula syntax:

``` r

da_results <- calc_differential_abundance(
  gut_microbiome,
  formula = ~ treatment + age,
  contrast = "treatment",
  methods = c("consensus", "caft", "linda", "clr_linear", "wilcoxon")
)

da_results |>
  filter(is_significant) |>
  select(taxon_id, Phylum, Genus, log2fc, padj_consensus, agreement_score) |>
  head(6)
#> # A tibble: 6 × 6
#>   taxon_id Phylum       Genus       log2fc padj_consensus agreement_score
#>   <chr>    <chr>        <chr>        <dbl>          <dbl>           <int>
#> 1 ASV05    Bacteroidota Bacteroides  -2.98       4.63e-15               4
#> 2 ASV04    Bacteroidota Bacteroides  -2.59       4.63e-15               4
#> 3 ASV10    Bacteroidota Prevotella   -2.37       4.63e-15               4
#> 4 ASV02    Bacteroidota Bacteroides   2.12       4.63e-15               4
#> 5 ASV01    Bacteroidota Bacteroides   2.03       4.63e-15               4
#> 6 ASV03    Bacteroidota Bacteroides   1.67       4.63e-15               4
```

- **CAFT** (bioRxiv 2025): Likelihood ratio testing separating zero-cell
  dropout from abundance shifts.
- **LinDA**: Linear models with rank-based reference taxon selection.
- **CLR-Linear**: Multiple linear regression on centered log-ratio
  abundances.
- **Cauchy Combination Test**: Aggregates $`p`$-values across distinct
  engines without assuming independence.

------------------------------------------------------------------------

## 9. Constrained Ordination (db-RDA) & Co-Occurrence Networks

### 9.1 Distance-Based Redundancy Analysis (db-RDA)

Test and visualize how environmental metadata variables constrain
microbial community variations:

``` r

dbrda_model <- run_dbrda(gut_microbiome, formula = ~ treatment + age, distance = "bray")
dbrda_model
#> -- tidybiome_dbrda (Distance-Based Redundancy Analysis) --
#>   * Constrained axes: 2
#>   * Inertia explained: dbRDA1 = 77.84%, dbRDA2 = 1.01%
#>   * Samples: 24
```

### 9.2 Microbial Co-Occurrence Networks

Build network graphs with correlation thresholds, FDR adjustment, and
node degree centrality:

``` r

net_model <- calc_network(
  gut_microbiome,
  method = "spearman",
  min_prevalence = 0.25,
  r_cutoff = 0.35,
  p_cutoff = 0.05
)

net_model
#> -- tidybiome_network (Microbial Co-Occurrence Network) --
#>   * Nodes (Taxa): 40
#>   * Significant edges: 483 (|r| >= 0.35, padj <= 0.05)
#>   * Positive: 400, Negative: 83
```

------------------------------------------------------------------------

## 10. Ecological Screening

Screen communities for prevalence, dominance, baseline divergence, and
host-microbe cross-associations:

``` r

# Taxon prevalence and mean abundance
calc_prevalence(gut_microbiome) |> head(4)
#> # A tibble: 4 × 8
#>   taxon_id prevalence prevalence_n mean_abundance median_abundance Kingdom 
#>   <chr>         <dbl>        <int>          <dbl>            <dbl> <chr>   
#> 1 ASV03             1           24          1712.            1447  Bacteria
#> 2 ASV08             1           24           675.             654. Bacteria
#> 3 ASV25             1           24           532.             508  Bacteria
#> 4 ASV20             1           24           526.             521  Bacteria
#> # ℹ 2 more variables: Phylum <chr>, Genus <chr>

# Dominant taxon per sample
calc_dominant(gut_microbiome, add_to_metadata = FALSE)
#> # A tibble: 2 × 3
#>   dominant_taxa n_samples frequency
#>   <chr>             <int>     <dbl>
#> 1 ASV01                23    0.958 
#> 2 ASV24                 1    0.0417

# Sample divergence relative to Control group
div_res <- calc_divergence(gut_microbiome, reference = list(treatment = "Control"))
mean(div_res$divergence[div_res$treatment == "Control"])
#> [1] 0.07758495
mean(div_res$divergence[div_res$treatment == "Treated"])
#> [1] 0.3550663

# Cross-associations between taxa and clinical variables
calc_cross_association(gut_microbiome, variables = c("age", "depth")) |> head(4)
#> # A tibble: 4 × 8
#>   taxon_id variable correlation p_value   padj Kingdom  Phylum       Genus      
#>   <chr>    <chr>          <dbl>   <dbl>  <dbl> <chr>    <chr>        <chr>      
#> 1 ASV02    age           -0.594 0.00221 0.0757 Bacteria Bacteroidota Bacteroides
#> 2 ASV30    depth          0.574 0.00333 0.0757 Bacteria Firmicutes   Blautia    
#> 3 ASV03    age           -0.569 0.00375 0.0757 Bacteria Bacteroidota Bacteroides
#> 4 ASV01    age           -0.557 0.00474 0.0757 Bacteria Bacteroidota Bacteroides
```

------------------------------------------------------------------------

## 11. Interoperability with Ecosystem Packages

Convert effortlessly between `tidybiome`, `vegan`, `phyloseq`, and
`mia`:

``` r

# Convert to vegan community matrix and metadata
veg <- to_vegan(gut_microbiome)

# Convert to phyloseq S4 object
ps <- to_phyloseq(gut_microbiome)

# Convert to Bioconductor TreeSummarizedExperiment
tse <- to_tse(gut_microbiome)

# Convert back to tidy_microbiome
tb_from_ps  <- as_tidybiome(ps)
tb_from_tse <- as_tidybiome(tse)
```

------------------------------------------------------------------------

## 12. Summary

`tidybiome` combines modern compositional statistical rigor, high
computational performance, and seamless tidyverse ergonomics into a
unified downstream microbiome analysis platform.
