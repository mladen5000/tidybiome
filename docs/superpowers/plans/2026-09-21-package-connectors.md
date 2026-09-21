# Ecosystem Connectors & Essential mia Functions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement full ecosystem connectors (`phyloseq`, `mia` / `TreeSummarizedExperiment`, `vegan`, `ape`) and the most important missing `mia` analytical functions (`calc_prevalence`, `calc_dominant`, `filter_prevalent`, `calc_divergence`, `calc_cross_association`) with strict TDD.

**Architecture:**
- `R/bridges.R`: Bidirectional data converters (`to_vegan`, `from_vegan`, `to_phyloseq`, `as_phyloseq`, `to_mia`, `as_mia`, `set_tree`, `get_tree`).
- `R/connectors_vegan.R`: Tidy analytical wrappers for `vegan` (`run_permanova`, `run_betadisper`, `run_nmds`).
- `R/prevalence.R`: Prevalence & dominance screening (`calc_prevalence`, `calc_dominant`, `filter_prevalent`).
- `R/divergence.R`: Ecological divergence against reference/control state (`calc_divergence`).
- `R/associations.R`: Taxon $\times$ covariate correlation analysis (`calc_cross_association`).

**Tech Stack:** R (>= 4.0), `tibble`, `dplyr`, `rlang`, `vegan`, `ape`, `testthat`, `roxygen2`.

## Global Constraints

- Pure tidyverse ergonomics: outputs are clean tibbles or updated `tidy_microbiome` objects.
- Strict TDD: Write failing test -> verify failure -> write minimal code -> verify pass -> commit.
- Zero test failures, warnings, or skips.

---

### Task 1: Vegan Bidirectional Converters & S3 Methods (`R/bridges.R`)

- [ ] **Step 1: Write failing test** (`tests/testthat/test-connectors-vegan.R`)
- [ ] **Step 2: Verify test fails**
- [ ] **Step 3: Implement minimal code** (`to_vegan`, `from_vegan`, `print.tidybiome_vegan`, `as_tidybiome.tidybiome_vegan`)
- [ ] **Step 4: Verify test passes**
- [ ] **Step 5: Commit**

---

### Task 2: Tidy PERMANOVA, Betadisper & NMDS Wrappers (`R/connectors_vegan.R`)

- [ ] **Step 1: Write failing test** (`tests/testthat/test-connectors-vegan-engines.R`)
- [ ] **Step 2: Verify test fails**
- [ ] **Step 3: Implement minimal code** (`run_permanova`, `run_betadisper`, `run_nmds`)
- [ ] **Step 4: Verify test passes**
- [ ] **Step 5: Commit**

---

### Task 3: Phylogenetic Tree Integration & Bridge Aliases (`R/bridges.R`)

- [ ] **Step 1: Write failing test** (`tests/testthat/test-connectors-tree-aliases.R`)
- [ ] **Step 2: Verify test fails**
- [ ] **Step 3: Implement minimal code** (`set_tree`, `get_tree`, `as_phyloseq`, `as_mia`, `to_mia`)
- [ ] **Step 4: Verify test passes**
- [ ] **Step 5: Commit**

---

### Task 4: Prevalence, Dominance & Filtering Functions (`R/prevalence.R`)

- [ ] **Step 1: Write failing test** (`tests/testthat/test-prevalence-dominance.R`)
- [ ] **Step 2: Verify test fails**
- [ ] **Step 3: Implement minimal code** (`calc_prevalence`, `calc_dominant`, `filter_prevalent`)
- [ ] **Step 4: Verify test passes**
- [ ] **Step 5: Commit**

---

### Task 5: Sample Divergence from Reference Cohort (`R/divergence.R`)

- [ ] **Step 1: Write failing test** (`tests/testthat/test-divergence.R`)
- [ ] **Step 2: Verify test fails**
- [ ] **Step 3: Implement minimal code** (`calc_divergence`)
- [ ] **Step 4: Verify test passes**
- [ ] **Step 5: Commit**

---

### Task 6: Cross-Association Analysis (`R/associations.R`)

- [ ] **Step 1: Write failing test** (`tests/testthat/test-associations.R`)
- [ ] **Step 2: Verify test fails**
- [ ] **Step 3: Implement minimal code** (`calc_cross_association`)
- [ ] **Step 4: Verify test passes**
- [ ] **Step 5: Commit**

---

### Task 7: Full Verification, Documentation & Showcase Update

- [ ] **Step 1: Update documentation and export in NAMESPACE**
- [ ] **Step 2: Run complete unit test suite across all modules**
- [ ] **Step 3: Update `demo/tidybiome_showcase.R` and execute end-to-end**
- [ ] **Step 4: Commit**
