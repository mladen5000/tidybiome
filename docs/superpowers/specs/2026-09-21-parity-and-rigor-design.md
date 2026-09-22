# Technical Specification: Parity, Rigor, and Scalability for tidybiome

**Date**: 2026-09-21  
**Author**: Antigravity  
**Status**: Approved  
**Target Package**: `tidybiome` (v0.2.0)

---

## 1. Executive Summary & Goals

This specification defines architectural and algorithmic enhancements to bring `tidybiome` to feature and rigor parity with Bioconductor's `mia` and CRAN's `phyloseq`. It systematically resolves:
1. **Data Integrity Hazards**: Ensuring S3 tibble subsetting (`[`, `filter`, `slice`, `distinct`) automatically and atomically synchronizes assay matrices and phylogenetic trees.
2. **Computational Bottlenecks**: Replacing $O(N^2 \cdot P)$ nested R loops in distance calculations with compiled C or vectorized matrix algebra; optimizing `tidy_abundance(long = TRUE)` memory consumption.
3. **Statistical Modeling Deficits**: Providing full formula-based covariate adjustment (`~ treatment + age + batch`), continuous predictor modeling, and multi-group contrasts in differential abundance.
4. **Phylogenetic Ecology Gap**: Implementing native **Unweighted and Weighted UniFrac** distance metrics and preserving phylogenetic trees across community transformations.
5. **Advanced Ecological Modeling**: Introducing distance-based redundancy analysis (db-RDA) constrained ordinations and microbial co-occurrence network analysis.

---

## 2. Component Architecture & Detailed Design

### 2.1 Component 1: Data Integrity, Synchronization, and Indexing

#### S3 Method Overloads
* Implement `[.tidy_microbiome`:
  * Detect row subsetting (`i`):
    * Slices sample metadata tibble.
    * Slices column dimension of all matrices in `attr(tb, "assays")` to match remaining `sample_id`s.
    * Re-validates that `sample_id` exists. If dropped, restores it.
  * Detect column subsetting (`j`):
    * Slices metadata columns while preserving `sample_id` and all attributes (`assays`, `tax_table`, `phy_tree`, `metadata`).
* Implement `[<-.tidy_microbiome` to preserve attributes on column assignment.
* Implement `filter_taxa(tb, ...)`:
  * Accepts standard tidyverse expressions evaluated in the context of `attr(tb, "tax_table")`.
  * Computes matching `taxon_id`s.
  * Slices all assay rows, tax_table rows, and prunes `attr(tb, "phy_tree")` tips via `ape::keep.tip`.
* Implement `validate_tidy_microbiome(tb)`:
  * Enforces exact synchronization between `sample_id` in metadata and assay column names.
  * Enforces exact synchronization between `taxon_id` in `tax_table` and assay row names.
  * Validates tip labels in `phy_tree` against taxon IDs when a tree is present.

---

### 2.2 Component 2: High-Performance Vectorized Distances & Memory Optimization

#### Vectorized Distance Metrics in `R/beta_diversity.R`
* **Bray-Curtis (`calc_bray_curtis`)**:
  * Check if `vegan` is available $\rightarrow$ delegate to `vegan::vegdist(t(mat), method = "bray")`.
  * Pure-R fallback: Compute total-sum scaling `prop = sweep(mat, 2, colSums(mat), "/")`, then compute `stats::dist(t(prop), method = "manhattan") * 0.5`. This transforms $O(N^2 \cdot P)$ loop operations into an optimized C-level BLAS/LAPACK Manhattan distance.
* **Jaccard (`calc_jaccard`)**:
  * Check if `vegan` is available $\rightarrow$ delegate to `vegan::vegdist(t(mat), method = "jaccard", binary = TRUE)`.
  * Pure-R fallback: Binarize matrix `pa = (mat > 0) * 1`. Compute intersection matrix via `C = crossprod(pa)`. Union is $U_{jk} = n_j + n_k - C_{jk}$. Distance is $1 - (C / U)$.
* **Jensen-Shannon Divergence (`calc_jsd`)**:
  * Vectorize calculation across columns using matrix multiplication and log proportions, avoiding sample-pair R iteration.
* **Lean `tidy_abundance()` in `R/extractors.R`**:
  * Eliminate `expand.grid()` and multi-pass `merge()`.
  * Pre-allocate vectors:
    * `taxon_id = rep(rownames(mat), times = ncol(mat))`
    * `sample_id = rep(colnames(mat), each = nrow(mat))`
    * `abundance = as.vector(mat)`
  * Direct vector indexing for metadata and taxonomy merging to eliminate memory explosion.

---

### 2.3 Component 3: Rigorous Covariate-Adjusted Differential Abundance

#### Formula-Driven Modeling in `R/differential_abundance.R`
* Upgrade `calc_differential_abundance()`:
  * Signature:
    ```r
    calc_differential_abundance(
      tb,
      formula = NULL,
      group = NULL,
      covariates = NULL,
      contrast = NULL,
      methods = c("clr_linear", "linda", "wilcoxon", "consensus"),
      fdr_cutoff = 0.05,
      assay = "counts",
      pseudocount = 0.5
    )
    ```
  * Formula Parsing:
    * Supports formulas like `~ treatment + age + batch` or `~ bmi + sex`.
    * Automatically extracts target predictor (first term or user-specified `contrast`) and adjusts for all remaining covariates using `stats::model.matrix()`.
  * Support for Continuous Predictors:
    * If target variable is numeric/continuous (e.g. BMI, age), fits linear models across taxa: $\text{CLR}(Y_i) \sim X_{\text{continuous}} + Z_{\text{covariates}}$, returning slope $\beta$, standard error, t-statistic, and adjusted p-value.
  * Support for Multi-Level Categorical Factors:
    * Fits omnibus ANOVA / F-test across factor levels plus pairwise contrast coefficients against reference level.
  * Rigorous Effect Size & Multiple Testing:
    * Benjamini-Hochberg FDR adjustment across all models.
    * Real LinDA-style weighted least squares and shrinkage for compositional correction when covariates are present.

---

### 2.4 Component 4: Phylogenetic Ecology Suite (UniFrac & Tree Integration)

#### UniFrac Distances in `R/beta_diversity.R`
* Implement `calc_unifrac(tb, weighted = FALSE, normalized = TRUE, assay = "counts")`:
  * Requires an attached `ape::phylo` tree.
  * Computes the tip-to-root ancestor-edge presence/abundance matrix:
    * For each branch $e$ in the tree with length $l_e$, calculates the total abundance of descendant tips in sample $j$: $p_{e, j}$.
    * **Unweighted UniFrac**:
      $$U = \frac{\sum_e l_e \cdot |I(p_{e, j} > 0) - I(p_{e, k} > 0)|}{\sum_e l_e \cdot I(p_{e, j} + p_{e, k} > 0)}$$
    * **Weighted UniFrac**:
      $$W = \frac{\sum_e l_e \cdot |p_{e, j} - p_{e, k}|}{\sum_e l_e \cdot (p_{e, j} + p_{e, k})}$$
  * Integrate `"unifrac"` and `"wunifrac"` directly into `calc_beta_diversity(tb, metric = "unifrac")` and `calc_beta_diversity(tb, metric = "wunifrac")`.
* **Tree Preservation during Agglomeration**:
  * In `aggregate_taxa()`:
    * When collapsing by rank, map the taxonomic rank to representative subtree tips or prune the tree to keep one tip per aggregated group, preserving phylogenetic topology rather than blindly dropping the tree.
* **Tip Glomming**:
  * Introduce `glom_taxa_tree(tb, h = 0.1)` to cluster tips by phylogenetic cophenetic distance.

---

### 2.5 Component 5: Constrained Ordination & Co-Occurrence Networks

#### Constrained Ordination in `R/connectors_vegan.R` & `R/ordination.R`
* Implement `run_dbrda(tb, formula, assay = "counts", distance = "raitchison", ...)`:
  * Runs Distance-based Redundancy Analysis (`vegan::dbrda` or `vegan::capscale`).
  * Returns tidy coordinates:
    * `samples`: Sample coordinates in constrained space (`dbRDA1`, `dbRDA2`).
    * `biplot`: Environmental variable vectors / arrows.
    * `variance_explained`: Proportion of inertia explained by constrained vs unconstrained axes.
* Provide `plot_dbrda(dbrda_res, color = ...)` for biplot visualization with environmental constraint vectors.

#### Microbial Co-Occurrence Networks in `R/associations.R` & `R/plots.R`
* Implement `calc_network(tb, assay = "counts", method = "spearman", min_prevalence = 0.2, r_cutoff = 0.5, p_cutoff = 0.05)`:
  * Computes taxon-taxon correlation matrix.
  * Filters edges by correlation magnitude $|r| \ge r_{\text{cutoff}}$ and $p_{\text{adj}} \le p_{\text{cutoff}}$.
  * Returns an S3 object `tidybiome_network` with tidy `nodes` (degrees, mean abundance, taxonomy) and `edges` (from, to, weight, direction).
* Implement `plot_network(network_res)`:
  * Visualizes co-occurrence network with node sizes proportional to abundance, colored by taxonomy (e.g. Phylum), and edges colored by positive/negative correlation.

---

## 3. Verification & Test Plan

1. **Unit Tests (`tests/testthat/`)**:
   * `test-container.R`: Verify that `tb[1:3, ]`, `tb[, c("sample_id", "group")]`, `filter_taxa()`, and `validate_tidy_microbiome()` maintain 100% matrix synchronization.
   * `test-beta-diversity.R`: Verify that vectorized Bray-Curtis and Jaccard exactly match vegan output; verify unweighted and weighted UniFrac distances on trees.
   * `test-differential-abundance.R`: Verify formula interface with continuous variables, multiple covariates, and multi-group factors.
   * `test-ordination.R`: Verify db-RDA coordinates, biplot vectors, and inertia breakdown.
   * `test-associations.R`: Verify network nodes, edges, and degree calculations.
2. **End-to-End Validation**:
   * Execute full showcase script (`demo/tidybiome_showcase.R`).
   * Run full test suite with 0 failures, 0 errors, 0 warnings.
