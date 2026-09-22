# Changelog

## tidybiome 0.1.0

Initial production release of **`tidybiome`**, a modern,
tidyverse-native Swiss Army knife for downstream microbiome analysis.

### Major Features

#### Container Architecture & Data Integrity

- S3 `tidy_microbiome` container inheriting directly from `tbl_df`,
  enabling standard `dplyr` operations on sample metadata while
  synchronizing multi-assay matrices, taxonomic lineage tables, and
  phylogenetic trees.
- Full S3 matrix subsetting (`[.tidy_microbiome` and
  `[<-.tidy_microbiome`) maintaining sample ID matching across assays
  and tree tips.
- [`filter_taxa()`](https://tidybiome.org/reference/filter_taxa.md) for
  taxonomic expression filtering with automatic tree pruning via
  [`ape::keep.tip()`](https://rdrr.io/pkg/ape/man/drop.tip.html).
- [`validate_tidy_microbiome()`](https://tidybiome.org/reference/validate_tidy_microbiome.md)
  diagnostic validator verifying 1:1 synchronization between metadata,
  assays, taxonomy, and phylogenetic trees.
- Intuitive [`assay()`](https://tidybiome.org/reference/assay.md) and
  [`tax_table()`](https://tidybiome.org/reference/tax_table.md) S3
  accessors for seamless interoperability with legacy Bioconductor and
  phyloseq workflows.

#### Modern Normalization & Compositional Transformations

- Robust Centered Log-Ratio (`transform_abundance(method = "rclr")`,
  Martino et al. 2019) eliminating zero-count pseudocount distortion.
- Geometric Mean of Pairwise Ratios
  ([`calc_gmpr_size_factors()`](https://tidybiome.org/reference/calc_gmpr_size_factors.md),
  `transform_abundance(method = "gmpr")`, Chen et al. 2018) for
  zero-tolerant size factor normalization in sparse cohorts.
- Hellinger transformation (`transform_abundance(method = "hellinger")`)
  for linear ordination of compositional community data.
- Total Sum Scaling (`"relabundance"`), Coverage-based standardization
  (`"coverage"`), and classic CLR (`"clr"`).

#### Diversity & Vectorized Distance Metrics

- Unified **Hill Numbers Profile**
  (`calc_alpha_diversity(metrics = "hill")`) calculating $`q=0`$
  (Richness), $`q=1`$ (Exponential Shannon), and $`q=2`$ (Inverse
  Simpson) on the intuitive scale of effective species.
- Simplex Compositional Variation (Aitchison total variance on the
  closed simplex, bioRxiv 2026).
- High-performance vectorized **Unweighted & Weighted UniFrac**
  ([`calc_unifrac()`](https://tidybiome.org/reference/calc_unifrac.md))
  via ancestor-descendant matrix multiplication.
- Vectorized Bray-Curtis, Jaccard, and Jensen-Shannon distance routines.
- Optimal transport Tree-Wasserstein Earth Mover’s Distance and
  Distance-based Test for Homogeneity
  ([`dth_test()`](https://tidybiome.org/reference/dth_test.md), bioRxiv
  2025).

#### Statistical Modeling & Differential Abundance

- Full formula interface in
  `calc_differential_abundance(formula = ~ treatment + age + batch)`
  supporting multi-covariate regression adjustments and continuous
  response predictors.
- CAFT (Compositional Analysis with Fractional Thresholding, bioRxiv
  Dec 2025) two-part likelihood ratio engine separating zero-cell
  dropout from abundance shifts.
- LinDA (Linear Models for Differential Abundance, Zhou et al. 2022)
  with rank-based reference taxon selection.
- CLR-linear multiple regression models and stratified Wilcoxon rank-sum
  testing.
- Multi-engine consensus architecture aggregating $`p`$-values via the
  Cauchy combination test with agreement scoring.

#### Ecological Modeling & Screening

- Distance-based Redundancy Analysis
  ([`run_dbrda()`](https://tidybiome.org/reference/run_dbrda.md),
  [`plot_dbrda()`](https://tidybiome.org/reference/plot_dbrda.md)) for
  constrained ordination of non-Euclidean ecological distance matrices
  with biplot environmental constraint vectors.
- Microbial Co-Occurrence Networks
  ([`calc_network()`](https://tidybiome.org/reference/calc_network.md),
  [`plot_network()`](https://tidybiome.org/reference/plot_network.md))
  with correlation thresholds, FDR control, and degree centrality.
- Ecological screening utilities:
  [`calc_prevalence()`](https://tidybiome.org/reference/calc_prevalence.md),
  [`filter_prevalent()`](https://tidybiome.org/reference/filter_prevalent.md),
  [`calc_dominant()`](https://tidybiome.org/reference/calc_dominant.md),
  [`calc_divergence()`](https://tidybiome.org/reference/calc_divergence.md),
  and
  [`calc_cross_association()`](https://tidybiome.org/reference/calc_cross_association.md).

#### Pipeline Ingestion Suite

- [`import_dada2()`](https://tidybiome.org/reference/import_dada2.md):
  Ingests DADA2 sequence tables and taxonomy matrices, automatically
  cleaning long nucleotide names into concise `ASV001` identifiers while
  preserving sequence strings in taxonomy annotations.
- [`import_qiime2()`](https://tidybiome.org/reference/import_qiime2.md):
  Imports QIIME 2 TSV feature tables, parses semicolon taxonomy strings,
  and attaches Newick phylogenetic trees.
- [`import_metaphlan()`](https://tidybiome.org/reference/import_metaphlan.md):
  Parses MetaPhlAn 3/4 merged abundance tables, splitting hierarchical
  clade pipe strings into standard taxonomic rank columns.

#### Publication Aesthetics & Standalone Dashboards

- [`plot_composition()`](https://tidybiome.org/reference/plot_composition.md):
  Stacked bar charts with intelligent top-$`N`$ pooling and soft-gray
  `"Other"` baseline bars.
- [`plot_heatmap()`](https://tidybiome.org/reference/plot_heatmap.md):
  Compositional heatmaps with hierarchical clustering of samples and
  taxa, multi-scale transforms, and faceted metadata annotations.
- [`plot_ordination()`](https://tidybiome.org/reference/plot_ordination.md):
  RPCA, cPCA, and PCoA scatter plots with 95% confidence ellipses and
  driving taxon biplot vectors.
- [`report_tidybiome()`](https://tidybiome.org/reference/report_tidybiome.md):
  Zero-dependency standalone HTML quality control dashboard summarizing
  cohort dimensions, sequencing depth distributions, zero sparsity, top
  taxa, and metadata dictionaries.
- Minimalist publication theme:
  [`theme_tidybiome()`](https://tidybiome.org/reference/theme_tidybiome.md).

#### Sample Quality Control & Diagnostic Suite

- `calc_qc_metrics(tb, augment = TRUE)`: Calculates per-sample library
  depth (`qc_total_reads`), observed richness (`qc_n_features`),
  sample-level sparsity (`qc_sparsity`), dominance of the top taxon
  (`qc_top_taxon_share`), and Shannon entropy (`qc_shannon`). Seamlessly
  augments sample metadata for downstream dplyr filtering pipelines.
- [`plot_qc()`](https://tidybiome.org/reference/plot_qc.md):
  Publication-ready QC diagnostic scatter plot displaying sequencing
  depth versus observed feature richness with log scaling, metadata
  color mappings, and cutoff thresholds.
- [`summary.tidy_microbiome()`](https://tidybiome.org/reference/summary.tidy_microbiome.md):
  Comprehensive S3 ecological summary method computing library size
  quantiles, global count matrix sparsity, taxonomic rank inventories,
  phylogenetic tree topology, and metadata dictionaries.

#### Ecosystem Interoperability

- Bidirectional conversion bridges:
  - [`to_vegan()`](https://tidybiome.org/reference/to_vegan.md) and
    [`from_vegan()`](https://tidybiome.org/reference/from_vegan.md) with
    wrappers for
    [`run_permanova()`](https://tidybiome.org/reference/run_permanova.md)
    and
    [`run_betadisper()`](https://tidybiome.org/reference/run_betadisper.md).
  - [`to_phyloseq()`](https://tidybiome.org/reference/to_phyloseq.md),
    [`as_phyloseq()`](https://tidybiome.org/reference/to_phyloseq.md),
    and
    [`as_tidybiome()`](https://tidybiome.org/reference/as_tidybiome.md).
  - [`to_tse()`](https://tidybiome.org/reference/to_tse.md),
    [`to_mia()`](https://tidybiome.org/reference/to_tse.md),
    [`as_mia()`](https://tidybiome.org/reference/to_tse.md), and
    [`as_tidybiome()`](https://tidybiome.org/reference/as_tidybiome.md).
  - [`to_s7()`](https://tidybiome.org/reference/to_s7.md),
    [`from_s7()`](https://tidybiome.org/reference/from_s7.md), and S7
    generic dispatch (`as_tidybiome`) with formal `tidy_microbiome_s7`
    class definition supporting next-generation R OOP.

### Performance & Scalability Optimizations

- **Vectorized Cross-Associations
  ([`calc_cross_association()`](https://tidybiome.org/reference/calc_cross_association.md))**:
  - Replaced per-taxon and per-variable nested loops and repetitive
    [`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html)
    invocations with BLAS matrix operations
    ([`stats::cor()`](https://rdrr.io/r/stats/cor.html)) and closed-form
    analytic Student’s $`t`$-statistic transformation.
  - Implemented single-pass matrix rank transforms for Spearman
    correlations, accelerating association testing across thousands of
    features by over 50x.
- **Candidate-Pruned Co-Occurrence Networks
  ([`calc_network()`](https://tidybiome.org/reference/calc_network.md))**:
  - Filtered candidate edges at the upper-triangular index level prior
    to candidate edge tibble materialization, avoiding $`O(T^2)`$ memory
    explosion on large OTU/ASV tables.
  - Retained conservative whole-matrix FDR corrections on candidate
    sets.
- **Skew-Symmetric GMPR Size Factors
  ([`calc_gmpr_size_factors()`](https://tidybiome.org/reference/calc_gmpr_size_factors.md))**:
  - Halved pairwise median ratio iterations by exploiting skew-symmetry
    ($`r_{ji} = 1 / r_{ij}`$) and precomputed positive-count incidence
    masks.
- **C-Compiled Manhattan UniFrac Distances
  ([`calc_unifrac()`](https://tidybiome.org/reference/calc_unifrac.md))**:
  - Formulated both unweighted and weighted UniFrac numerators as
    C-compiled Manhattan distances
    (`stats::dist(..., method = "manhattan")`) over branch-weighted
    occurrence and abundance matrices, eliminating nested R sample
    loops.
- **Vectorized QR Decomposition for Differential Abundance
  ([`calc_differential_abundance()`](https://tidybiome.org/reference/calc_differential_abundance.md))**:
  - Precomputes a single QR decomposition of the design matrix in LinDA
    and CLR-linear regression engines.
  - Simultaneously solves coefficients, residuals, standard errors, and
    multi-level $`F`$-tests for all taxa via BLAS operations
    ([`qr.coef()`](https://rdrr.io/r/base/qr.html),
    [`qr.resid()`](https://rdrr.io/r/base/qr.html)), eliminating
    per-taxon [`stats::lm()`](https://rdrr.io/r/stats/lm.html) fitting
    loops.
- **Compiled C Agglomeration
  ([`aggregate_taxa()`](https://tidybiome.org/reference/aggregate_taxa.md))**:
  - Utilizes `base::rowsum(..., reorder = FALSE)` for taxonomic feature
    aggregation, executing in compiled C rather than R loops.
- **Memory-Efficient Long-Format Extraction
  ([`tidy_abundance()`](https://tidybiome.org/reference/tidy_abundance.md))**:
  - Added targeted `taxa` and `samples` subsetting, and optional
    metadata/taxonomy joins (`include_metadata`, `include_taxonomy`) to
    prevent unintended multi-million-row tibble materializations.
- **Ancestor-Incidence Matrix Faith’s Phylogenetic Diversity
  (`calc_faith_pd()`)**:
  - Replaced repeated $`O(N)`$ tree-pruning via
    [`ape::keep.tip()`](https://rdrr.io/pkg/ape/man/drop.tip.html) with
    a single ancestor-descendant tip incidence matrix multiplication
    against sample presence masks, computing Faith’s PD across all
    samples in a single vectorized step.
