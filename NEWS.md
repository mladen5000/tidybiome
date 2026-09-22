# tidybiome 0.1.0

Initial production release of **`tidybiome`**, a modern, tidyverse-native Swiss Army knife for downstream microbiome analysis.

## Major Features

### Container Architecture & Data Integrity
* S3 `tidy_microbiome` container inheriting directly from `tbl_df`, enabling standard `dplyr` operations on sample metadata while synchronizing multi-assay matrices, taxonomic lineage tables, and phylogenetic trees.
* Full S3 matrix subsetting (`[.tidy_microbiome` and `[<-.tidy_microbiome`) maintaining sample ID matching across assays and tree tips.
* `filter_taxa()` for taxonomic expression filtering with automatic tree pruning via `ape::keep.tip()`.
* `validate_tidy_microbiome()` diagnostic validator verifying 1:1 synchronization between metadata, assays, taxonomy, and phylogenetic trees.
* Intuitive `assay()` and `tax_table()` S3 accessors for seamless interoperability with legacy Bioconductor and phyloseq workflows.

### Modern Normalization & Compositional Transformations
* Robust Centered Log-Ratio (`transform_abundance(method = "rclr")`, Martino et al. 2019) eliminating zero-count pseudocount distortion.
* Geometric Mean of Pairwise Ratios (`calc_gmpr_size_factors()`, `transform_abundance(method = "gmpr")`, Chen et al. 2018) for zero-tolerant size factor normalization in sparse cohorts.
* Hellinger transformation (`transform_abundance(method = "hellinger")`) for linear ordination of compositional community data.
* Total Sum Scaling (`"relabundance"`), Coverage-based standardization (`"coverage"`), and classic CLR (`"clr"`).

### Diversity & Vectorized Distance Metrics
* Unified **Hill Numbers Profile** (`calc_alpha_diversity(metrics = "hill")`) calculating $q=0$ (Richness), $q=1$ (Exponential Shannon), and $q=2$ (Inverse Simpson) on the intuitive scale of effective species.
* Simplex Compositional Variation (Aitchison total variance on the closed simplex, bioRxiv 2026).
* High-performance vectorized **Unweighted & Weighted UniFrac** (`calc_unifrac()`) via ancestor-descendant matrix multiplication.
* Vectorized Bray-Curtis, Jaccard, and Jensen-Shannon distance routines.
* Optimal transport Tree-Wasserstein Earth Mover's Distance and Distance-based Test for Homogeneity (`dth_test()`, bioRxiv 2025).

### Statistical Modeling & Differential Abundance
* Full formula interface in `calc_differential_abundance(formula = ~ treatment + age + batch)` supporting multi-covariate regression adjustments and continuous response predictors.
* CAFT (Compositional Analysis with Fractional Thresholding, bioRxiv Dec 2025) two-part likelihood ratio engine separating zero-cell dropout from abundance shifts.
* LinDA (Linear Models for Differential Abundance, Zhou et al. 2022) with rank-based reference taxon selection.
* CLR-linear multiple regression models and stratified Wilcoxon rank-sum testing.
* Multi-engine consensus architecture aggregating $p$-values via the Cauchy combination test with agreement scoring.

### Ecological Modeling & Screening
* Distance-based Redundancy Analysis (`run_dbrda()`, `plot_dbrda()`) for constrained ordination of non-Euclidean ecological distance matrices with biplot environmental constraint vectors.
* Microbial Co-Occurrence Networks (`calc_network()`, `plot_network()`) with correlation thresholds, FDR control, and degree centrality.
* Ecological screening utilities: `calc_prevalence()`, `filter_prevalent()`, `calc_dominant()`, `calc_divergence()`, and `calc_cross_association()`.

### Pipeline Ingestion Suite
* `import_dada2()`: Ingests DADA2 sequence tables and taxonomy matrices, automatically cleaning long nucleotide names into concise `ASV001` identifiers while preserving sequence strings in taxonomy annotations.
* `import_qiime2()`: Imports QIIME 2 TSV feature tables, parses semicolon taxonomy strings, and attaches Newick phylogenetic trees.
* `import_metaphlan()`: Parses MetaPhlAn 3/4 merged abundance tables, splitting hierarchical clade pipe strings into standard taxonomic rank columns.

### Publication Aesthetics & Standalone Dashboards
* `plot_composition()`: Stacked bar charts with intelligent top-$N$ pooling and soft-gray `"Other"` baseline bars.
* `plot_heatmap()`: Compositional heatmaps with hierarchical clustering of samples and taxa, multi-scale transforms, and faceted metadata annotations.
* `plot_ordination()`: RPCA, cPCA, and PCoA scatter plots with 95% confidence ellipses and driving taxon biplot vectors.
* `report_tidybiome()`: Zero-dependency standalone HTML quality control dashboard summarizing cohort dimensions, sequencing depth distributions, zero sparsity, top taxa, and metadata dictionaries.
* Minimalist publication theme: `theme_tidybiome()`.

### Sample Quality Control & Diagnostic Suite
* `calc_qc_metrics(tb, augment = TRUE)`: Calculates per-sample library depth (`qc_total_reads`), observed richness (`qc_n_features`), sample-level sparsity (`qc_sparsity`), dominance of the top taxon (`qc_top_taxon_share`), and Shannon entropy (`qc_shannon`). Seamlessly augments sample metadata for downstream dplyr filtering pipelines.
* `plot_qc()`: Publication-ready QC diagnostic scatter plot displaying sequencing depth versus observed feature richness with log scaling, metadata color mappings, and cutoff thresholds.
* `summary.tidy_microbiome()`: Comprehensive S3 ecological summary method computing library size quantiles, global count matrix sparsity, taxonomic rank inventories, phylogenetic tree topology, and metadata dictionaries.

### Ecosystem Interoperability
* Bidirectional conversion bridges:
  - `to_vegan()` and `from_vegan()` with wrappers for `run_permanova()` and `run_betadisper()`.
  - `to_phyloseq()`, `as_phyloseq()`, and `as_tidybiome()`.
  - `to_tse()`, `to_mia()`, `as_mia()`, and `as_tidybiome()`.
  - `to_s7()`, `from_s7()`, and S7 generic dispatch (`as_tidybiome`) with formal `tidy_microbiome_s7` class definition supporting next-generation R OOP.

