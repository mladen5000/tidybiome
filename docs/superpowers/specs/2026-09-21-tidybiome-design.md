# Specification: `tidybiome` R Package

## 1. Overview & Vision
`tidybiome` is a tidyverse-native, high-aesthetic R package designed as a modern "Swiss Army knife" for downstream microbiome analysis. While drawing structural and functional inspiration from Bioconductor's `mia` and legacy `phyloseq`, `tidybiome` is built to be:

1. **Simpler and Tidy-first**: Microbiome objects behave directly like tibbles in user workflows (`dplyr` verbs work natively), avoiding cumbersome S4 slot manipulation.
2. **Up-to-Date with bioRxiv & Modern Literature**: Emphasizes modern, compositional, and consensus methodologies (Robust CLR, Robust Aitchison, Unified Hill Numbers profile, RPCA biplots, and Consensus Differential Abundance) while maintaining full backward-compatibility with classic methods (Bray-Curtis, Shannon, classic CLR).
3. **Publication-Grade Visuals**: Replaces clunky, crowded plots with modern, highly aesthetic ggplot2 visualizations featuring intelligent top-$N$ taxonomic pooling, colorblind-safe palettes, confidence ellipses, and minimalist styling via `theme_tidybiome()`.

---

## 2. Core Architecture & Data Container

### 2.1 The `tidy_microbiome` S3 Class
The central data structure is `tidy_microbiome`, an S3 class inheriting from `tbl_df`, `tbl`, and `data.frame`.

* **Outer Appearance**: When printed in the R console, it displays as an enriched sample metadata tibble with informative summary header information (number of samples, number of taxa, total read depth range, active assays, and taxonomic ranks available).
* **Inner Attributes**:
  * `assays`: A named list of numeric matrices (taxa in rows, samples in columns). The primary raw count assay is `"counts"`. Transformed assays (e.g., `"relabundance"`, `"rclr"`, `"clr"`) are added dynamically.
  * `tax_table`: A tibble or data frame mapping taxon identifiers (`taxon_id`) to taxonomic ranks (`Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species`).
  * `phy_tree`: An optional phylogenetic tree (object of class `phylo` from `ape`).
  * `metadata`: A list for storing experiment-level provenance, distance matrices, ordination objects, and differential abundance results.

### 2.2 Tidyverse Method Dispatch
`tidybiome` implements S3 methods for core `dplyr` generics:
* `filter.tidy_microbiome()`: Filters sample rows. Automatically slices corresponding sample columns in all matrices in `assays`.
* `mutate.tidy_microbiome()`: Adds or modifies sample-level metadata variables.
* `select.tidy_microbiome()`: Selects sample metadata columns, always retaining the sample identifier.
* `arrange.tidy_microbiome()`: Re-orders samples, synchronizing order across sample metadata and matrix columns.
* `slice.tidy_microbiome()`: Slices samples by integer index.

### 2.3 Extraction Helpers
* `tidy_taxa(tb)`: Returns the taxonomy table as a clean, standardized tibble.
* `tidy_abundance(tb, assay = "counts", rank = NULL, long = TRUE)`: Returns abundance data optionally rolled up to a specific taxonomic rank, in either wide format or tidy long format (`sample_id`, `taxon_id`, `abundance`, plus taxonomy and sample metadata).

### 2.4 Interoperability Bridges
* `tidy_microbiome(counts, sample_data = NULL, tax_table = NULL, phy_tree = NULL)`: Direct constructor from matrices/data frames.
* `as_tidybiome(x)`: S3 generic with methods for:
  * `phyloseq` (from `phyloseq`)
  * `TreeSummarizedExperiment` / `SummarizedExperiment` (from Bioconductor)
  * `matrix` / `data.frame`
* `to_phyloseq(tb)`: Converts a `tidy_microbiome` object into a `phyloseq::phyloseq` object.
* `to_tse(tb)`: Converts a `tidy_microbiome` object into a Bioconductor `TreeSummarizedExperiment::TreeSummarizedExperiment`.

---

## 3. Analytical Methods: Modern & Classic

### 3.1 Normalization & Compositional Transformations (`transform_abundance`)
Function: `transform_abundance(tb, method = c("rclr", "clr", "relabundance", "log10", "pa", "coverage"), ...)`

* **Modern (bioRxiv / CoDA standards)**:
  * **Robust CLR (`rclr`)**: Computes the centered log-ratio using geometric means calculated strictly over non-zero values (Martino et al. / DEICODE), avoiding spurious correlations caused by pseudocounts.
  * **Coverage-based Standardization (`coverage`)**: Standardizes sequencing depth based on sample sample-coverage estimation rather than destructive random thinning.
* **Classic**:
  * **Total Sum Scaling / Relative Abundance (`relabundance`)**: Proportion of total sample reads ($x_{ij} / \sum_k x_{kj}$).
  * **Classic CLR (`clr`)**: Standard centered log-ratio with user-configurable pseudocount ($\log(x_{ij} + \delta) - \text{mean}(\log(x_{\cdot j} + \delta))$).
  * **Log10 (`log10`)**: $\log_{10}(x + 1)$.
  * **Presence/Absence (`pa`)**: Binary indicator matrix ($x > 0$).

### 3.2 Alpha Diversity (`calc_alpha_diversity`)
Function: `calc_alpha_diversity(tb, metrics = c("hill", "shannon", "simpson", "inv_simpson", "observed", "chao1", "faith_pd"), ...)`

Adds calculated indices directly as new columns to the sample metadata of `tb`.

* **Modern**:
  * **Unified Hill Numbers Profile (`hill`)**: Calculates effective species numbers for orders $q \in \{0, 1, 2\}$:
    * $q=0$: Species richness.
    * $q=1$: Exponential of Shannon entropy ($\exp(H)$).
    * $q=2$: Inverse Simpson index ($1 / \sum p_i^2$).
    All three indices exist on the identical intuitive scale ("effective number of species").
* **Classic**:
  * Shannon index, Gini-Simpson ($1 - D$), Inverse Simpson ($1/D$), Observed richness, Chao1 bias-corrected richness estimator.
  * **Faith's Phylogenetic Diversity (`faith_pd`)**: Sum of branch lengths on the phylogenetic tree connecting observed taxa (when `phy_tree` is present).

### 3.3 Beta Diversity & Distances (`calc_beta_diversity`)
Function: `calc_beta_diversity(tb, metric = c("raitchison", "aitchison", "bray", "jaccard", "unifrac", "wunifrac", "jsd"), assay = "counts")`

Returns an updated `tidy_microbiome` with the computed distance matrix stored in `tb@metadata$distances[[metric]]` or returns the `dist` object directly.

* **Modern**:
  * **Robust Aitchison Distance (`raitchison`)**: Euclidean distance computed on the `rclr` assay. Compositionally coherent, zero-robust, and avoids pseudocount distortions.
  * **Jensen-Shannon Divergence (`jsd`)**: Symmetrical information-theoretic divergence between probability distributions.
  * **Generalized / Weighted UniFrac (`wunifrac`, `unifrac`)**: Distance accounting for phylogenetic relationships between lineages.
* **Classic**:
  * **Bray-Curtis (`bray`)**: Standard ecological dissimilarity.
  * **Jaccard (`jaccard`)**: Presence/absence dissimilarity.
  * **Classic Aitchison (`aitchison`)**: Euclidean distance on standard pseudocount-CLR matrix.

### 3.4 Ordination (`calc_ordination`)
Function: `calc_ordination(tb, method = c("rpca", "pcoa", "nmds", "pca"), metric = "raitchison", ...)`

* **Modern**:
  * **Robust PCA (`rpca`)**: Matrix decomposition (singular value decomposition) directly on the zero-filtered `rclr` matrix. Generates both sample coordinates and taxon loading vectors for true compositional biplots without distortion.
  * **PCoA with Lingoes / Cailliez Correction (`pcoa`)**: Principal coordinates analysis with automatic correction for non-Euclidean negative eigenvalues.
* **Classic**:
  * Non-metric Multidimensional Scaling (`nmds`) via `vegan::metaMDS`.
  * Standard Principal Component Analysis (`pca`).

### 3.5 Differential Abundance Testing (`calc_differential_abundance`)
Function: `calc_differential_abundance(tb, group, formula = NULL, methods = c("consensus", "linda", "clr_linear", "wilcoxon"), ...)`

Returns a tidy tibble with one row per taxon, containing effect sizes (log2 fold change), per-method p-values, combined p-values, and consensus metrics.

* **Modern (bioRxiv Consensus Engine)**:
  * **Consensus Engine (`consensus`)**: Inspired by recent 2024–2025 bioRxiv preprints (*ConsensusMetaDA*, *dar*, *LinDA*). Instead of relying on a single method prone to false discoveries:
    1. Runs **LinDA** (Linear Models for Differential Abundance, robust to compositionality and zero inflation).
    2. Runs **CLR-Linear / Moderated Regression** (linear model on CLR/rCLR abundance).
    3. Runs **Non-parametric Wilcoxon / Rank-sum** test.
    4. Combines p-values using Cauchy or Fisher combination tests.
    5. Calculates an **Agreement Score** (number of engines agreeing on direction and significance at FDR $\le 0.05$) and ranks taxa by consensus evidence.
* **Classic**:
  * Standard two-sample Wilcoxon rank-sum test or Welch's t-test on relative abundance or log-transformed counts.

### 3.6 Core Microbiome Analysis (`calc_core_microbiome`)
Function: `calc_core_microbiome(tb, method = c("inflection", "threshold"), min_prevalence = 0.8, min_abundance = 0.001)`

* **Modern**:
  * **Data-driven Inflection Curve (`inflection`)**: Evaluates taxon prevalence across a continuous spectrum of abundance thresholds and identifies the natural mathematical inflection points defining the true core community without arbitrary cutoffs.
* **Classic**:
  * **Static Threshold (`threshold`)**: Filters taxa meeting predefined minimum prevalence (e.g. 80% of samples) and relative abundance (e.g. $> 0.1\%$).

---

## 4. Visualizations & Aesthetics

All plotting functions return clean, un-rendered `ggplot` objects, allowing users to add additional ggplot2 layers (`+ labs(...)`, `+ facet_wrap(...)`, etc.).

### 4.1 Visual Functions
1. **`plot_composition(tb, rank = "Genus", top_n = 8, group_by = NULL, palette = "tidybiome")`**:
   * Intelligent top-$N$ aggregation: Ranks taxa by mean relative abundance, displays the top $N$, and automatically groups all remaining taxa into an aesthetically styled `"Other"` group at the bottom.
   * Eliminates border clutter, ensures clean vertical/horizontal bars, and applies colorblind-safe palettes.
2. **`plot_ordination(tb, ordination = NULL, color = NULL, shape = NULL, ellipse = TRUE, biplot = FALSE, top_taxa = 5)`**:
   * Plots sample ordination coordinates with customizable 95% normal/t confidence ellipses.
   * If `biplot = TRUE`, overlays sleek arrows representing the top driving taxa from RPCA loadings.
3. **`plot_alpha(tb, metric = "hill_1", x = NULL, color = NULL, test = TRUE)`**:
   * Publication-grade hybrid plot: half-violin + boxplot + jittered sample points.
   * Automatically computes and annotates pairwise statistical brackets (Wilcoxon or ANOVA/Tukey) if `test = TRUE`.
4. **`plot_da_volcano(da_res, fdr_cutoff = 0.05, fc_cutoff = 1.0, label_top = 8)`**:
   * Volcano plot of effect size (log2FC) versus $-\log_{10}(\text{p-value})$.
   * Highlights consensus biomarkers with automatic label repulsion using clean point annotations.
5. **`plot_core(tb, rank = "Genus")`**:
   * Plots prevalence versus relative abundance curves across taxa, visually demarcating the detected core community.
6. **`theme_tidybiome()`**:
   * Minimalist, modern ggplot2 theme: crisp typography, light gray gridlines, transparent panel backgrounds, and neat facet labels.

---

## 5. Directory Structure of the R Package

```
tidybiome/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── tidy_microbiome.R       # S3 class constructor, validators, print, dplyr methods
│   ├── extractors.R            # tidy_taxa, tidy_abundance, agglomeration
│   ├── bridges.R               # as_tidybiome, to_phyloseq, to_tse
│   ├── transformations.R       # transform_abundance (rclr, clr, relabundance, coverage)
│   ├── alpha_diversity.R       # calc_alpha_diversity (hill numbers, shannon, simpson, faith_pd)
│   ├── beta_diversity.R        # calc_beta_diversity (raitchison, bray, jaccard, unifrac, jsd)
│   ├── ordination.R            # calc_ordination (rpca, pcoa, nmds)
│   ├── differential_abundance.R# calc_differential_abundance (consensus engine, linda, clr_lm)
│   ├── core_microbiome.R       # calc_core_microbiome (inflection, threshold)
│   ├── plots.R                 # plot_composition, plot_ordination, plot_alpha, plot_da_volcano, plot_core
│   ├── theme.R                 # theme_tidybiome, curated color palettes
│   └── data_demo.R             # built-in synthetic & benchmark microbiome datasets
├── data/                       # Built-in demo dataset (gut_microbiome)
├── man/                        # Rd documentation files (generated via roxygen2)
├── tests/
│   ├── testthat.R
│   └── testthat/
│       ├── test-container.R
│       ├── test-transformations.R
│       ├── test-diversity.R
│       ├── test-differential-abundance.R
│       └── test-plots.R
└── vignettes/
    └── introduction_to_tidybiome.Rmd
```

---

## 6. Verification & Validation Plan
1. **Automated Unit Tests (`testthat`)**:
   * Container integrity: verifies that `filter()`, `mutate()`, and `select()` maintain strict alignment between sample metadata and count matrix columns.
   * Mathematical verification:
     * Hill numbers ($q=0, 1, 2$) reproduce expected values on known synthetic communities.
     * rCLR handles zero-heavy matrices without returning `NaN` or `Inf`.
     * Robust Aitchison distance produces symmetric, valid positive semi-definite distances.
     * RPCA generates valid eigenvalues and orthogonal eigenvector axes.
     * Consensus DA accurately detects spiked differential taxa and combines p-values.
2. **End-to-End Demo Script**:
   * A full walkthrough script demonstrating import -> filtering -> rCLR transformation -> Hill diversity -> RPCA biplot -> Consensus DA -> Volcano plot.
