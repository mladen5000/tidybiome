# Parity, Rigor, and Scalability Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Elevate `tidybiome` to scientific, mathematical, and computational parity with `mia` and `phyloseq` by adding robust subsetting synchronization, high-speed vectorized distance engines, formula-based covariate differential abundance, UniFrac phylogenetic metrics, db-RDA constrained ordinations, and co-occurrence networks.

**Architecture:** Maintain the lightweight S3 tibble inheritance while strengthening data integrity with S3 subsetting methods (`[.tidy_microbiome`), replacing nested R loops with matrix algebra / vegan acceleration, implementing formula parsing (`stats::model.matrix`) for differential abundance, and computing branch-ancestor UniFrac matrices via `ape`.

**Tech Stack:** R (>= 4.1.0), tibble, dplyr, rlang, ggplot2, stats, utils, vegan, ape, testthat (>= 3.0.0).

## Global Constraints
- Must pass `R CMD check` cleanly without errors or warnings.
- Backward compatibility: Existing 161 test assertions must continue to pass without regression.
- S3 class `tidy_microbiome` must maintain strict sample and feature synchronization across all operations.

---

### Task 1: Data Integrity & Robust Indexing

**Files:**
- Modify: `R/tidy_microbiome.R`
- Modify: `NAMESPACE`
- Test: `tests/testthat/test-container.R`

**Interfaces:**
- Produces: `[.tidy_microbiome(x, i, j, ...)`, `[<-.tidy_microbiome(x, i, j, value)`, `filter_taxa(tb, ...)`, `validate_tidy_microbiome(tb)`

- [ ] **Step 1: Write failing tests for S3 bracket subsetting and `filter_taxa`**

Add tests in `tests/testthat/test-container.R`:
```r
test_that("bracket subsetting slices sample metadata and assay matrices synchronously", {
  counts <- matrix(1:6, nrow = 2, ncol = 3, dimnames = list(c("T1", "T2"), c("S1", "S2", "S3")))
  tb <- tidy_microbiome(counts)
  sub_tb <- tb[1:2, ]
  expect_s3_class(sub_tb, "tidy_microbiome")
  expect_equal(nrow(sub_tb), 2)
  expect_equal(colnames(attr(sub_tb, "assays")$counts), c("S1", "S2"))
})

test_that("filter_taxa slices features across assays, taxonomy, and tree", {
  counts <- matrix(1:6, nrow = 2, ncol = 3, dimnames = list(c("T1", "T2"), c("S1", "S2", "S3")))
  tax <- tibble::tibble(taxon_id = c("T1", "T2"), Phylum = c("Bacteroidota", "Firmicutes"))
  tb <- tidy_microbiome(counts, tax_table = tax)
  filt_tb <- filter_taxa(tb, Phylum == "Bacteroidota")
  expect_equal(nrow(attr(filt_tb, "assays")$counts), 1)
  expect_equal(rownames(attr(filt_tb, "assays")$counts), "T1")
  expect_equal(attr(filt_tb, "tax_table")$taxon_id, "T1")
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-container.R")'`
Expected: FAIL with "object 'filter_taxa' not found" or assay matrix mismatch.

- [ ] **Step 3: Implement S3 bracket subsetting, `filter_taxa`, and `validate_tidy_microbiome`**

In `R/tidy_microbiome.R`:
Implement `[.tidy_microbiome`, `[<-.tidy_microbiome`, `filter_taxa`, and `validate_tidy_microbiome`.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-container.R")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/tidy_microbiome.R tests/testthat/test-container.R NAMESPACE
git commit -m "feat(container): add bracket subsetting, filter_taxa, and container validation"
```

---

### Task 2: High-Performance Vectorized Distances & Memory Optimization

**Files:**
- Modify: `R/beta_diversity.R`
- Modify: `R/extractors.R`
- Test: `tests/testthat/test-beta-diversity.R`
- Test: `tests/testthat/test-extractors-bridges.R`

**Interfaces:**
- Modifies: `calc_bray_curtis()`, `calc_jaccard()`, `calc_jsd()` to use vectorized matrix operations or `vegan::vegdist`.
- Modifies: `tidy_abundance(long = TRUE)` to use direct vector indexing instead of `expand.grid` and `merge`.

- [ ] **Step 1: Write test verifying distance calculation consistency and speed**

In `tests/testthat/test-beta-diversity.R`:
```r
test_that("vectorized bray and jaccard match vegan vegdist exactly", {
  counts <- matrix(c(10, 0, 5, 20, 10, 0, 5, 5, 5), nrow = 3, ncol = 3,
                   dimnames = list(c("A", "B", "C"), c("S1", "S2", "S3")))
  tb <- tidy_microbiome(counts)
  tb_bray <- calc_beta_diversity(tb, metric = "bray")
  d_bray <- get_distance(tb_bray, "bray")
  expect_equal(as.matrix(d_bray), as.matrix(vegan::vegdist(t(counts), method = "bray")), tolerance = 1e-6)
})
```

- [ ] **Step 2: Implement vectorized distance calculation and lean `tidy_abundance`**

Update `calc_bray_curtis`, `calc_jaccard`, `calc_jsd` in `R/beta_diversity.R`.
Update `tidy_abundance` in `R/extractors.R` using preallocated vectors.

- [ ] **Step 3: Run tests to verify correctness**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-beta-diversity.R")'`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add R/beta_diversity.R R/extractors.R tests/testthat/test-beta-diversity.R
git commit -m "perf(beta): vectorize bray/jaccard/jsd distance calculations and optimize tidy_abundance"
```

---

### Task 3: Rigorous Covariate-Adjusted Differential Abundance

**Files:**
- Modify: `R/differential_abundance.R`
- Modify: `NAMESPACE`
- Test: `tests/testthat/test-differential-abundance.R`

**Interfaces:**
- Produces: `calc_differential_abundance(tb, formula = NULL, group = NULL, covariates = NULL, contrast = NULL, ...)`

- [ ] **Step 1: Write tests for formula-based differential abundance with covariates and continuous predictors**

In `tests/testthat/test-differential-abundance.R`:
```r
test_that("calc_differential_abundance adjusts for covariates using formula", {
  set.seed(42)
  counts <- matrix(rpois(100, lambda = 50), nrow = 10, ncol = 10,
                   dimnames = list(paste0("ASV", 1:10), paste0("S", 1:10)))
  sample_data <- data.frame(
    sample_id = paste0("S", 1:10),
    treatment = rep(c("A", "B"), each = 5),
    age = rnorm(10, mean = 40, sd = 5)
  )
  tb <- tidy_microbiome(counts, sample_data = sample_data)
  res <- calc_differential_abundance(tb, formula = ~ treatment + age, contrast = "treatment")
  expect_s3_class(res, "tbl_df")
  expect_true("p_consensus" %in% colnames(res))
  expect_true("padj_consensus" %in% colnames(res))
})

test_that("calc_differential_abundance supports continuous target variable", {
  set.seed(42)
  counts <- matrix(rpois(100, lambda = 50), nrow = 10, ncol = 10,
                   dimnames = list(paste0("ASV", 1:10), paste0("S", 1:10)))
  sample_data <- data.frame(
    sample_id = paste0("S", 1:10),
    bmi = rnorm(10, mean = 25, sd = 4)
  )
  tb <- tidy_microbiome(counts, sample_data = sample_data)
  res <- calc_differential_abundance(tb, formula = ~ bmi)
  expect_s3_class(res, "tbl_df")
  expect_true("log2fc" %in% colnames(res))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-differential-abundance.R")'`
Expected: FAIL.

- [ ] **Step 3: Implement formula parser and covariate regression in `calc_differential_abundance`**

Update `R/differential_abundance.R` to parse `formula`, extract predictor and covariates, fit multiple linear regression on CLR and relative abundances, and compute consensus statistics.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-differential-abundance.R")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/differential_abundance.R tests/testthat/test-differential-abundance.R NAMESPACE
git commit -m "feat(da): implement formula interface with covariate adjustment and continuous predictor support"
```

---

### Task 4: Phylogenetic Ecology Suite (UniFrac & Tree Integration)

**Files:**
- Modify: `R/beta_diversity.R`
- Modify: `R/extractors.R`
- Modify: `NAMESPACE`
- Test: `tests/testthat/test-beta-diversity.R`
- Test: `tests/testthat/test-connectors-tree-aliases.R`

**Interfaces:**
- Produces: `calc_unifrac(tb, weighted = FALSE, normalized = TRUE, assay = "counts")`
- Enhances: `calc_beta_diversity(tb, metric = "unifrac")` and `calc_beta_diversity(tb, metric = "wunifrac")`
- Enhances: `aggregate_taxa()` to retain/prune trees

- [ ] **Step 1: Write test for unweighted and weighted UniFrac**

In `tests/testthat/test-beta-diversity.R`:
```r
test_that("calc_unifrac computes unweighted and weighted UniFrac distances on trees", {
  counts <- matrix(c(10, 0, 5, 20), nrow = 2, ncol = 2,
                   dimnames = list(c("T1", "T2"), c("S1", "S2")))
  tree <- ape::read.tree(text = "(T1:0.5,T2:0.5);")
  tb <- tidy_microbiome(counts, phy_tree = tree)
  
  tb_unifrac <- calc_beta_diversity(tb, metric = "unifrac")
  d_u <- get_distance(tb_unifrac, "unifrac")
  expect_s3_class(d_u, "dist")
  
  tb_wunifrac <- calc_beta_diversity(tb, metric = "wunifrac")
  d_w <- get_distance(tb_wunifrac, "wunifrac")
  expect_s3_class(d_w, "dist")
})
```

- [ ] **Step 2: Implement UniFrac algorithms and tree retention**

In `R/beta_diversity.R`, implement `calc_unifrac()` using branch ancestor projection.
Update `calc_beta_diversity()` to support `"unifrac"` and `"wunifrac"`.
In `R/extractors.R`, update `aggregate_taxa()` to prune `phy_tree` to group exemplars when possible.

- [ ] **Step 3: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-beta-diversity.R")'`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add R/beta_diversity.R R/extractors.R tests/testthat/test-beta-diversity.R NAMESPACE
git commit -m "feat(phylo): implement unweighted and weighted UniFrac and preserve trees during agglomeration"
```

---

### Task 5: Advanced Ecological Modeling (Constrained Ordination & Networks)

**Files:**
- Modify: `R/connectors_vegan.R`
- Modify: `R/associations.R`
- Modify: `R/plots.R`
- Modify: `NAMESPACE`
- Test: `tests/testthat/test-connectors-vegan-engines.R`
- Test: `tests/testthat/test-associations.R`

**Interfaces:**
- Produces: `run_dbrda(tb, formula, distance = "bray", ...)`, `plot_dbrda(dbrda_res, color = ...)`
- Produces: `calc_network(tb, method = "spearman", r_cutoff = 0.4, p_cutoff = 0.05, ...)`, `plot_network(network_res)`

- [ ] **Step 1: Write tests for `run_dbrda` and `calc_network`**

In `tests/testthat/test-connectors-vegan-engines.R` and `tests/testthat/test-associations.R`.

- [ ] **Step 2: Implement `run_dbrda` and `calc_network`**

In `R/connectors_vegan.R`: wrap `vegan::dbrda` / `vegan::capscale` returning tidy sample coordinates and biplot constraint vectors.
In `R/associations.R`: implement `calc_network()` to compute taxon co-occurrence correlations, degree centrality, and edge weights.
In `R/plots.R`: implement `plot_dbrda()` and `plot_network()`.

- [ ] **Step 3: Run tests to verify they pass**

Run test suites for connectors and associations.
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add R/connectors_vegan.R R/associations.R R/plots.R tests/testthat NAMESPACE
git commit -m "feat(ecology): add db-RDA constrained ordination and microbial co-occurrence networks"
```

---

### Task 6: End-to-End Showcase & Full Suite Verification

**Files:**
- Modify: `demo/tidybiome_showcase.R`
- Modify: `README.md`
- Test: All test files in `tests/testthat/`

- [ ] **Step 1: Update demo showcase with UniFrac, formula DA, db-RDA, and network analysis**
- [ ] **Step 2: Execute entire test suite**

Run: `Rscript -e 'testthat::test_dir("tests/testthat")'`
Expected: 100% passing tests, 0 failures, 0 errors.

- [ ] **Step 3: Run full showcase script**

Run: `Rscript demo/tidybiome_showcase.R`
Expected: Clean exit 0 with all generated figures and tables.

- [ ] **Step 4: Commit and finalize**

```bash
git add demo/tidybiome_showcase.R README.md
git commit -m "docs(showcase): update end-to-end showcase with parity and rigor features"
```
