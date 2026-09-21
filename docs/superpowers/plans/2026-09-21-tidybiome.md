# `tidybiome` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a complete, production-grade, tidyverse-native R package named `tidybiome` that acts as an up-to-date, aesthetic "Swiss Army knife" for downstream microbiome analysis, incorporating state-of-the-art bioRxiv methodologies (rCLR, rAitchison, Hill diversity profile, RPCA biplots, Consensus Differential Abundance) alongside classic methods.

**Architecture:** An S3 container `tidy_microbiome` that inherits from `tbl_df`, `tbl`, and `data.frame` so it natively supports `dplyr` verbs while keeping assays (matrices) and taxonomic hierarchies tightly synchronized. Clean functional verbs (`transform_abundance`, `calc_alpha_diversity`, `calc_beta_diversity`, `calc_ordination`, `calc_differential_abundance`, `calc_core_microbiome`) operate on the container and return enriched objects or tidy tibbles, paired with high-aesthetic `ggplot2` plotting functions.

**Tech Stack:** R (>= 4.1), `tibble`, `dplyr`, `rlang`, `ggplot2`, `vegan` (or base matrix linear algebra equivalents), `testthat` (>= 3.0).

## Global Constraints

- Must strictly adhere to tidy data principles: functions take `tidy_microbiome` or tibble and return `tidy_microbiome` or tibble.
- Keep dependency footprint lean: base R linear algebra (SVD, eigen, dist) where possible; avoid heavy S4 Bioconductor dependencies for core usage, while providing seamless interoperability wrappers when `phyloseq` or `TreeSummarizedExperiment` are present.
- Visuals must be publication-ready: clean typography, automatic top-N grouping with "Other" pool, colorblind-safe palettes, no cramped 50-color legends.
- Every function must include comprehensive roxygen2 documentation and `testthat` unit tests.

---

### Task 1: Package Scaffolding & Core S3 Container

**Files:**
- Create: `DESCRIPTION`
- Create: `NAMESPACE`
- Create: `R/tidy_microbiome.R`
- Create: `tests/testthat.R`
- Create: `tests/testthat/test-container.R`

**Interfaces:**
- Produces: `tidy_microbiome(counts, sample_data = NULL, tax_table = NULL, phy_tree = NULL)`
- Produces: S3 methods `filter.tidy_microbiome`, `mutate.tidy_microbiome`, `select.tidy_microbiome`, `arrange.tidy_microbiome`, `print.tidy_microbiome`

- [ ] **Step 1: Write the failing tests for container creation and dplyr verbs**

```r
# tests/testthat/test-container.R
test_that("tidy_microbiome initializes and displays correctly", {
  counts <- matrix(c(10, 0, 5, 20, 15, 2, 0, 8), nrow = 2, ncol = 4,
                   dimnames = list(c("ASV1", "ASV2"), c("S1", "S2", "S3", "S4")))
  sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"),
                            group = c("A", "A", "B", "B"))
  tax_table <- data.frame(taxon_id = c("ASV1", "ASV2"),
                          Phylum = c("Bacteroidota", "Firmicutes"),
                          Genus = c("Bacteroides", "Lactobacillus"))
  
  tb <- tidy_microbiome(counts, sample_data, tax_table)
  expect_s3_class(tb, "tidy_microbiome")
  expect_s3_class(tb, "tbl_df")
  expect_equal(nrow(tb), 4)
  expect_equal(attr(tb, "assays")$counts, counts)
})

test_that("dplyr verbs synchronize sample metadata and assays", {
  counts <- matrix(c(10, 0, 5, 20, 15, 2, 0, 8), nrow = 2, ncol = 4,
                   dimnames = list(c("ASV1", "ASV2"), c("S1", "S2", "S3", "S4")))
  sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"),
                            group = c("A", "A", "B", "B"),
                            depth = c(15, 22, 15, 10))
  tb <- tidy_microbiome(counts, sample_data)
  
  # Filter
  tb_filt <- dplyr::filter(tb, group == "A")
  expect_equal(nrow(tb_filt), 2)
  expect_equal(colnames(attr(tb_filt, "assays")$counts), c("S1", "S2"))
  
  # Mutate
  tb_mut <- dplyr::mutate(tb, log_depth = log(depth))
  expect_true("log_depth" %in% colnames(tb_mut))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-container.R")'`
Expected: FAIL ("could not find function 'tidy_microbiome'")

- [ ] **Step 3: Implement `DESCRIPTION`, `NAMESPACE`, and `R/tidy_microbiome.R`**

Create `DESCRIPTION` with package metadata, imports: `tibble`, `dplyr`, `rlang`.
Implement `tidy_microbiome()` constructor with validation (checking sample names match column names of counts, and taxon names match rownames of counts). Implement `filter.tidy_microbiome`, `mutate.tidy_microbiome`, `select.tidy_microbiome`, `arrange.tidy_microbiome`, and formatted `print.tidy_microbiome`.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-container.R")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add DESCRIPTION NAMESPACE R/tidy_microbiome.R tests/
git commit -m "feat: implement tidy_microbiome S3 container with dplyr dispatch"
```

---

### Task 2: Extractors, Agglomeration, and Bioconductor Bridges

**Files:**
- Create: `R/extractors.R`
- Create: `R/bridges.R`
- Test: `tests/testthat/test-extractors-bridges.R`

**Interfaces:**
- Consumes: `tidy_microbiome`
- Produces: `tidy_taxa(tb)`, `tidy_abundance(tb, assay = "counts", rank = NULL, long = TRUE)`
- Produces: `aggregate_taxa(tb, rank = "Genus")`
- Produces: `as_tidybiome(x)`, `to_phyloseq(tb)`, `to_tse(tb)`

- [ ] **Step 1: Write the failing tests**

```r
# tests/testthat/test-extractors-bridges.R
test_that("tidy_abundance extracts long format and rolls up taxonomy", {
  counts <- matrix(c(10, 5, 20, 8), nrow = 2, ncol = 2,
                   dimnames = list(c("ASV1", "ASV2"), c("S1", "S2")))
  sample_data <- data.frame(sample_id = c("S1", "S2"), group = c("A", "B"))
  tax_table <- data.frame(taxon_id = c("ASV1", "ASV2"),
                          Phylum = c("Bacteroidota", "Bacteroidota"),
                          Genus = c("Bacteroides", "Prevotella"))
  tb <- tidy_microbiome(counts, sample_data, tax_table)
  
  long_df <- tidy_abundance(tb, assay = "counts", long = TRUE)
  expect_true(all(c("sample_id", "taxon_id", "abundance", "group", "Genus") %in% colnames(long_df)))
  
  # Agglomerate by Phylum
  tb_phylum <- aggregate_taxa(tb, rank = "Phylum")
  expect_equal(nrow(attr(tb_phylum, "assays")$counts), 1) # merged into 1 phylum
  expect_equal(attr(tb_phylum, "assays")$counts["Bacteroidota", "S1"], 15)
})
```

- [ ] **Step 2: Run test to verify failure**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-extractors-bridges.R")'`
Expected: FAIL

- [ ] **Step 3: Implement extractors, agglomeration, and bridges**

Implement `tidy_taxa()`, `tidy_abundance()`, `aggregate_taxa()` in `R/extractors.R`.
Implement `as_tidybiome()`, `to_phyloseq()`, and `to_tse()` in `R/bridges.R` with graceful degradation if optional packages are not loaded.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-extractors-bridges.R")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add R/extractors.R R/bridges.R tests/testthat/test-extractors-bridges.R
git commit -m "feat: add extractors, taxonomic agglomeration, and bioconductor bridges"
```

---

### Task 3: Normalization & Compositional Transformations (rCLR, CLR, Coverage)

**Files:**
- Create: `R/transformations.R`
- Test: `tests/testthat/test-transformations.R`

**Interfaces:**
- Consumes: `tidy_microbiome`
- Produces: `transform_abundance(tb, method = c("rclr", "clr", "relabundance", "log10", "pa", "coverage"), pseudocount = 1, name = NULL)`

- [ ] **Step 1: Write failing tests for modern rCLR and classic transformations**

```r
# tests/testthat/test-transformations.R
test_that("rclr handles zero-inflated matrices without NaNs or pseudocounts", {
  counts <- matrix(c(100, 0, 50, 0,
                     20,  10, 0, 5), nrow = 2, byrow = TRUE,
                   dimnames = list(c("ASV1", "ASV2"), c("S1", "S2", "S3", "S4")))
  tb <- tidy_microbiome(counts)
  
  tb_rclr <- transform_abundance(tb, method = "rclr")
  rclr_mat <- attr(tb_rclr, "assays")$rclr
  expect_false(any(is.nan(rclr_mat)))
  expect_false(any(is.infinite(rclr_mat)))
  # Zeros remain zero in robust CLR
  expect_equal(rclr_mat[1, 2], 0)
  expect_equal(rclr_mat[2, 3], 0)
  
  # Non-zeros are centered by geometric mean of positive values
  pos_vals <- counts[counts[, 1] > 0, 1]
  geom_mean <- exp(mean(log(pos_vals)))
  expect_equal(rclr_mat[1, 1], log(100 / geom_mean))
})

test_that("relabundance scales each sample to sum to 1", {
  counts <- matrix(c(10, 90, 50, 50), nrow = 2, dimnames = list(c("A", "B"), c("S1", "S2")))
  tb <- tidy_microbiome(counts)
  tb_rel <- transform_abundance(tb, method = "relabundance")
  expect_equal(colSums(attr(tb_rel, "assays")$relabundance), c(S1 = 1, S2 = 1))
})
```

- [ ] **Step 2: Run test to verify failure**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-transformations.R")'`
Expected: FAIL

- [ ] **Step 3: Implement `transform_abundance`**

Implement:
- Robust CLR (`rclr`): $\text{rclr}(x_i) = \log(x_i / g(x_{>0}))$ where $g(x_{>0}) = \exp\left(\frac{1}{|x_{>0}|} \sum_{j \in x_{>0}} \log x_j\right)$ for non-zero entries; zero entries mapped to 0.
- Classic CLR (`clr`): centered log-ratio with pseudocount.
- Relative abundance (`relabundance`): column proportions.
- Coverage-based standardization (`coverage`): scaling counts by estimated sample coverage $C = 1 - f_1 / n$.
- Log10 and Presence/Absence.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-transformations.R")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add R/transformations.R tests/testthat/test-transformations.R
git commit -m "feat: add robust CLR, coverage standardization, and classic transformations"
```

---

### Task 4: Alpha Diversity (Unified Hill Numbers Profile & Classic Indices)

**Files:**
- Create: `R/alpha_diversity.R`
- Test: `tests/testthat/test-alpha-diversity.R`

**Interfaces:**
- Consumes: `tidy_microbiome`
- Produces: `calc_alpha_diversity(tb, metrics = c("hill", "shannon", "simpson", "inv_simpson", "observed", "chao1"), assay = "counts")`

- [ ] **Step 1: Write failing tests for Hill numbers and classic alpha diversity**

```r
# tests/testthat/test-alpha-diversity.R
test_that("Hill numbers profile calculates q=0, 1, 2 on effective species scale", {
  # Even community with 4 species: effective species should all be 4
  counts <- matrix(c(25, 25, 25, 25), nrow = 4, ncol = 1,
                   dimnames = list(c("T1", "T2", "T3", "T4"), "S1"))
  tb <- tidy_microbiome(counts)
  tb_alpha <- calc_alpha_diversity(tb, metrics = c("hill", "shannon", "simpson"))
  
  expect_equal(tb_alpha$hill_0, 4) # Richness
  expect_equal(round(tb_alpha$hill_1, 5), 4) # Exp(Shannon)
  expect_equal(round(tb_alpha$hill_2, 5), 4) # Inverse Simpson
  expect_equal(round(tb_alpha$shannon, 5), round(log(4), 5))
})
```

- [ ] **Step 2: Run test to verify failure**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-alpha-diversity.R")'`
Expected: FAIL

- [ ] **Step 3: Implement `calc_alpha_diversity` in `R/alpha_diversity.R`**

Implement vectorized calculations for:
- Hill numbers:
  - $q=0$: $\sum [x_i > 0]$
  - $q=1$: $\exp(-\sum p_i \log p_i)$ where $p_i = x_i / \sum x$
  - $q=2$: $1 / \sum p_i^2$
- Shannon: $-\sum p_i \log p_i$
- Gini-Simpson: $1 - \sum p_i^2$
- Inverse Simpson: $1 / \sum p_i^2$
- Chao1: $S_{obs} + \frac{f_1(f_1 - 1)}{2(f_2 + 1)}$
Append results directly into the `tidy_microbiome` tibble.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-alpha-diversity.R")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add R/alpha_diversity.R tests/testthat/test-alpha-diversity.R
git commit -m "feat: add unified Hill numbers profile and classic alpha diversity metrics"
```

---

### Task 5: Beta Diversity & Distance Matrices (rAitchison, JSD, Bray-Curtis)

**Files:**
- Create: `R/beta_diversity.R`
- Test: `tests/testthat/test-beta-diversity.R`

**Interfaces:**
- Consumes: `tidy_microbiome`
- Produces: `calc_beta_diversity(tb, metric = c("raitchison", "aitchison", "bray", "jaccard", "jsd"), assay = "counts")`
- Produces: `get_distance(tb, metric = NULL)`

- [ ] **Step 1: Write failing tests for robust Aitchison and ecological distances**

```r
# tests/testthat/test-beta-diversity.R
test_that("Robust Aitchison distance computes Euclidean distance on rclr", {
  counts <- matrix(c(10, 0, 5,
                     20, 10, 0,
                     5,  5,  5), nrow = 3, ncol = 3,
                   dimnames = list(c("A", "B", "C"), c("S1", "S2", "S3")))
  tb <- tidy_microbiome(counts)
  tb_dist <- calc_beta_diversity(tb, metric = "raitchison")
  
  d <- get_distance(tb_dist, "raitchison")
  expect_s3_class(d, "dist")
  expect_equal(attr(d, "Size"), 3)
  expect_true(all(d >= 0))
})
```

- [ ] **Step 2: Run test to verify failure**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-beta-diversity.R")'`
Expected: FAIL

- [ ] **Step 3: Implement `calc_beta_diversity` and distance algorithms**

Implement native, fast distance calculations:
- `raitchison`: automated computation of rCLR followed by Euclidean distance.
- `aitchison`: standard CLR with pseudocount followed by Euclidean distance.
- `bray`: Bray-Curtis dissimilarity $d_{jk} = \frac{\sum |x_{ij} - x_{ik}|}{\sum (x_{ij} + x_{ik})}$.
- `jaccard`: Jaccard binary dissimilarity.
- `jsd`: Jensen-Shannon divergence between probability vectors.
Store distance objects in `attr(tb, "metadata")$distances[[metric]]`.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-beta-diversity.R")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add R/beta_diversity.R tests/testthat/test-beta-diversity.R
git commit -m "feat: implement robust Aitchison, JSD, and beta diversity distance suite"
```

---

### Task 6: Ordination & Biplots (RPCA & PCoA)

**Files:**
- Create: `R/ordination.R`
- Test: `tests/testthat/test-ordination.R`

**Interfaces:**
- Consumes: `tidy_microbiome`
- Produces: `calc_ordination(tb, method = c("rpca", "pcoa", "pca", "nmds"), metric = "raitchison", n_components = 3)`
- Produces: `get_ordination(tb, method = NULL)`

- [ ] **Step 1: Write failing tests for RPCA matrix decomposition and biplot coordinates**

```r
# tests/testthat/test-ordination.R
test_that("RPCA generates sample coordinates, taxon loadings, and explained variance", {
  counts <- matrix(c(100, 20, 5, 2,
                     5, 10, 80, 70,
                     1, 1, 50, 40), nrow = 3, byrow = TRUE,
                   dimnames = list(c("Tax1", "Tax2", "Tax3"), c("S1", "S2", "S3", "S4")))
  tb <- tidy_microbiome(counts)
  tb_ord <- calc_ordination(tb, method = "rpca")
  
  ord <- get_ordination(tb_ord, "rpca")
  expect_true(all(c("samples", "taxa", "variance_explained") %in% names(ord)))
  expect_equal(nrow(ord$samples), 4)
  expect_equal(nrow(ord$taxa), 3)
  expect_true(sum(ord$variance_explained) <= 1.0001)
})
```

- [ ] **Step 2: Run test to verify failure**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-ordination.R")'`
Expected: FAIL

- [ ] **Step 3: Implement `calc_ordination` in `R/ordination.R`**

- `rpca`: Performs singular value decomposition (`svd`) on the mean-centered rCLR matrix. Extracts sample factor scores ($U D$) and feature loading vectors ($V$), plus eigenvalues / explained variance ratio.
- `pcoa`: Performs spectral decomposition on Gower-centered distance matrices with automatic Lingoes or Cailliez correction for negative eigenvalues.
- `pca`: Standard PCA via SVD on normalized matrix.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-ordination.R")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add R/ordination.R tests/testthat/test-ordination.R
git commit -m "feat: implement RPCA compositional biplot ordination and corrected PCoA"
```

---

### Task 7: SOTA Differential Abundance & Consensus Engine

**Files:**
- Create: `R/differential_abundance.R`
- Test: `tests/testthat/test-differential-abundance.R`

**Interfaces:**
- Consumes: `tidy_microbiome`
- Produces: `calc_differential_abundance(tb, group, formula = NULL, methods = c("consensus", "linda", "clr_linear", "wilcoxon"), fdr_cutoff = 0.05)`

- [ ] **Step 1: Write failing tests for LinDA, CLR-linear, and Consensus DA**

```r
# tests/testthat/test-differential-abundance.R
test_that("calc_differential_abundance runs consensus engines and flags biomarkers", {
  # Taxon 1 strongly differential between groups A and B
  counts <- matrix(c(100, 120, 110,  5,   8,   6,
                     50,  45,  55,  48,  52,  50), nrow = 2, byrow = TRUE,
                   dimnames = list(c("DiffTaxon", "NullTaxon"), paste0("S", 1:6)))
  sample_data <- data.frame(sample_id = paste0("S", 1:6),
                            group = c("A", "A", "A", "B", "B", "B"))
  tb <- tidy_microbiome(counts, sample_data)
  
  da_res <- calc_differential_abundance(tb, group = "group", methods = c("consensus", "linda", "clr_linear", "wilcoxon"))
  expect_s3_class(da_res, "tbl_df")
  expect_true(all(c("taxon_id", "log2fc", "p_consensus", "padj_consensus", "agreement_score", "is_significant") %in% colnames(da_res)))
  
  # DiffTaxon should be detected as significant with high agreement
  diff_row <- da_res[da_res$taxon_id == "DiffTaxon", ]
  expect_true(diff_row$is_significant)
  expect_true(diff_row$agreement_score >= 2)
})
```

- [ ] **Step 2: Run test to verify failure**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-differential-abundance.R")'`
Expected: FAIL

- [ ] **Step 3: Implement modern DA engines and consensus combination**

Implement in `R/differential_abundance.R`:
1. **LinDA** (Linear Models for Differential Abundance): Runs linear models on log-transformed relative abundances with sample-specific compositional bias correction (centered mean shift adjustment).
2. **CLR-Linear**: Ordinary/moderated linear regression on CLR or rCLR transformed matrix.
3. **Wilcoxon / Rank-Sum**: Non-parametric two-sample test.
4. **Consensus Aggregator**:
   - Combines p-values across engines using Cauchy combination test (robust to correlated statistics) or Fisher's method.
   - Calculates Benjamini-Hochberg FDR adjustments.
   - Computes `agreement_score` ($\sum [p_{adj, m} \le \alpha \text{ and } \text{sign}(\beta_m) = \text{sign}(\beta)]$).

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-differential-abundance.R")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add R/differential_abundance.R tests/testthat/test-differential-abundance.R
git commit -m "feat: implement bioRxiv-inspired consensus differential abundance engine"
```

---

### Task 8: Core Microbiome Analysis (Data-Driven Inflection & Thresholds)

**Files:**
- Create: `R/core_microbiome.R`
- Test: `tests/testthat/test-core-microbiome.R`

**Interfaces:**
- Consumes: `tidy_microbiome`
- Produces: `calc_core_microbiome(tb, method = c("inflection", "threshold"), min_prevalence = 0.8, min_abundance = 0.001)`

- [ ] **Step 1: Write failing tests for data-driven core microbiome detection**

```r
# tests/testthat/test-core-microbiome.R
test_that("calc_core_microbiome detects core members accurately", {
  # Core taxon present in all 5 samples, transient taxon only in 1
  counts <- matrix(c(50, 40, 60, 55, 45,
                     0,  0,  10, 0,  0), nrow = 2, byrow = TRUE,
                   dimnames = list(c("CoreTaxon", "RareTaxon"), paste0("S", 1:5)))
  tb <- tidy_microbiome(counts)
  
  core_res <- calc_core_microbiome(tb, method = "threshold", min_prevalence = 0.8)
  expect_s3_class(core_res, "tbl_df")
  expect_true(all(c("taxon_id", "prevalence", "mean_abundance", "is_core") %in% colnames(core_res)))
  expect_true(core_res$is_core[core_res$taxon_id == "CoreTaxon"])
  expect_false(core_res$is_core[core_res$taxon_id == "RareTaxon"])
})
```

- [ ] **Step 2: Run test to verify failure**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-core-microbiome.R")'`
Expected: FAIL

- [ ] **Step 3: Implement `calc_core_microbiome`**

Implement:
- `method = "threshold"`: Static filter on prevalence ($\text{prevalence} \ge \text{min\_prevalence}$) and relative abundance.
- `method = "inflection"`: Evaluates the continuous cumulative distribution curve of prevalence vs abundance across all taxa and computes second-derivative inflection points to discover natural community boundaries without arbitrary thresholds.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-core-microbiome.R")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add R/core_microbiome.R tests/testthat/test-core-microbiome.R
git commit -m "feat: implement data-driven inflection and threshold core microbiome detection"
```

---

### Task 9: Publication-Grade Visuals & Aesthetic Palette Suite

**Files:**
- Create: `R/theme.R`
- Create: `R/plots.R`
- Test: `tests/testthat/test-plots.R`

**Interfaces:**
- Consumes: `tidy_microbiome`, `ggplot2`
- Produces: `theme_tidybiome()`, `scale_color_tidybiome()`, `scale_fill_tidybiome()`
- Produces: `plot_composition(tb, rank = "Genus", top_n = 8, group_by = NULL)`
- Produces: `plot_ordination(tb, method = "rpca", color = NULL, ellipse = TRUE, biplot = FALSE)`
- Produces: `plot_alpha(tb, metric = "hill_1", x = NULL, color = NULL, test = TRUE)`
- Produces: `plot_da_volcano(da_res, fdr_cutoff = 0.05, fc_cutoff = 1.0, label_top = 8)`
- Produces: `plot_core(core_res)`

- [ ] **Step 1: Write failing tests verifying all plotting functions return valid ggplot objects**

```r
# tests/testthat/test-plots.R
test_that("all plotting functions return valid ggplot2 objects", {
  counts <- matrix(c(100, 20, 5, 2,
                     5, 10, 80, 70), nrow = 2, byrow = TRUE,
                   dimnames = list(c("Tax1", "Tax2"), c("S1", "S2", "S3", "S4")))
  sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"),
                            group = c("A", "A", "B", "B"))
  tax_table <- data.frame(taxon_id = c("Tax1", "Tax2"),
                          Phylum = c("Bacteroidota", "Firmicutes"),
                          Genus = c("Bacteroides", "Faecalibacterium"))
  tb <- tidy_microbiome(counts, sample_data, tax_table)
  
  # plot_composition
  p_comp <- plot_composition(tb, rank = "Genus", top_n = 1)
  expect_s3_class(p_comp, "ggplot")
  
  # plot_alpha
  tb <- calc_alpha_diversity(tb, metrics = "hill")
  p_alpha <- plot_alpha(tb, metric = "hill_1", x = "group")
  expect_s3_class(p_alpha, "ggplot")
  
  # plot_ordination
  tb <- calc_ordination(tb, method = "rpca")
  p_ord <- plot_ordination(tb, method = "rpca", color = "group")
  expect_s3_class(p_ord, "ggplot")
})
```

- [ ] **Step 2: Run test to verify failure**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-plots.R")'`
Expected: FAIL

- [ ] **Step 3: Implement modern theme, color palettes, and plotting suite**

Implement in `R/theme.R`:
- `theme_tidybiome()`: Minimalist, clean publication theme (high-contrast sans font, gentle grey gridlines, no heavy black outer box, clean strip background).
- Curated color palettes: Nature-inspired, Okabe-Ito colorblind-friendly.

Implement in `R/plots.R`:
- `plot_composition()`: Automatically aggregates all taxa beyond `top_n` into a soft-gray `"Other"` category placed at the bottom, stacked relative abundance bars with clean spacing.
- `plot_ordination()`: PCoA or RPCA biplot with optional 95% confidence ellipses, sample centroids, and sleek loading arrows with repelled labels for top taxa.
- `plot_alpha()`: Raincloud / half-violin + boxplot + jittered sample points with optional automated pairwise statistical significance brackets.
- `plot_da_volcano()`: High-impact biomarker volcano plot with automatic text repulsion for top significant taxa.
- `plot_core()`: Continuous prevalence-abundance landscape plot.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-plots.R")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add R/theme.R R/plots.R tests/testthat/test-plots.R
git commit -m "feat: implement publication-grade aesthetic plots and theme_tidybiome"
```

---

### Task 10: Demo Dataset, Vignette, and End-to-End Verification Pipeline

**Files:**
- Create: `R/data_demo.R`
- Create: `vignettes/introduction_to_tidybiome.Rmd`
- Create: `demo/tidybiome_showcase.R`
- Modify: `DESCRIPTION`
- Test: `tests/testthat/test-e2e.R`

**Interfaces:**
- Consumes: All `tidybiome` modules
- Produces: Built-in dataset `data(gut_microbiome)`
- Produces: Complete end-to-end runnable showcase script demonstrating full workflow

- [ ] **Step 1: Write end-to-end integration test**

```r
# tests/testthat/test-e2e.R
test_that("full end-to-end tidybiome workflow executes seamlessly", {
  data("gut_microbiome", package = "tidybiome")
  expect_s3_class(gut_microbiome, "tidy_microbiome")
  
  # Tidy pipeline: filter -> transform -> alpha -> beta -> ordination -> DA
  res <- gut_microbiome |>
    dplyr::filter(depth > 500) |>
    transform_abundance(method = "rclr") |>
    calc_alpha_diversity(metrics = "hill") |>
    calc_beta_diversity(metric = "raitchison") |>
    calc_ordination(method = "rpca")
  
  expect_true("hill_1" %in% colnames(res))
  expect_true("raitchison" %in% names(attr(res, "metadata")$distances))
  expect_true("rpca" %in% names(attr(res, "metadata")$ordinations))
})
```

- [ ] **Step 2: Implement demo dataset generator and showcase vignette**

Create a realistic synthetic human gut microbiome benchmark dataset (`gut_microbiome`) comprising 30 samples (e.g. 15 Control vs. 15 Treatment), 50 ASVs across 4 phyla, with spiked differential abundance and realistic sparsity. Save in `data/` and document in `R/data_demo.R`.
Create `demo/tidybiome_showcase.R` and `vignettes/introduction_to_tidybiome.Rmd` documenting each updated bioRxiv method step-by-step.

- [ ] **Step 3: Run full package test suite**

Run: `Rscript -e 'devtools::test()' || Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: All tests PASS.

- [ ] **Step 4: Commit**

```bash
git add data/ R/data_demo.R vignettes/ demo/ tests/testthat/test-e2e.R
git commit -m "feat: add demo dataset, showcase script, vignette, and end-to-end tests"
```
