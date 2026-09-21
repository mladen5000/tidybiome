# Design Specification: Ecosystem Package Connectors for `tidybiome`

**Date:** 2026-09-21  
**Author:** Antigravity  
**Status:** Approved  
**Scope:** Bidirectional data converters and tidy analysis wrappers connecting `tidybiome` with `phyloseq`, `mia` / `TreeSummarizedExperiment`, `vegan`, and `ape`.

---

## 1. Objectives & Motivation

While `tidybiome` provides cutting-edge statistical methods and tidy S3 ergonomics, microbiome researchers often work in mixed ecosystems where datasets originate in or must be shared with established frameworks:
- **`phyloseq`**: The ubiquitous legacy standard across thousands of published studies.
- **`mia` / `TreeSummarizedExperiment`**: The Bioconductor standard for tree-aware multi-assay experiments.
- **`vegan`**: The golden standard for community ecology multivariate statistics (PERMANOVA, dispersion tests, NMDS, CCA).
- **`ape`**: The core R package for phylogenetic tree representations.

This specification details:
1. Bidirectional data conversion bridges (`to_*`, `from_*`, `as_tidybiome`).
2. High-level tidy wrappers executing `vegan` routines (`adonis2`, `betadisper`, `metaMDS`) directly on `tidy_microbiome` objects and returning clean tibbles.
3. Phylogenetic tree attachment and synchronization with `ape::phylo`.

---

## 2. Architecture & File Structure

```
tidybiome/
├── R/
│   ├── bridges.R            # S3 generics & bidirectional data converters (phyloseq, TreeSE/mia, vegan, ape)
│   └── connectors_vegan.R   # Tidy analysis wrappers (run_permanova, run_betadisper, run_nmds)
├── tests/testthat/
│   └── test-connectors.R    # Comprehensive unit tests for all connectors and analysis wrappers
└── man/                     # Auto-generated .Rd documentation
```

### Dependency Management
- **`vegan`** and **`ape`**: Listed in `Imports` or `Suggests` with graceful handling; both are available CRAN packages.
- **`phyloseq`** and **`TreeSummarizedExperiment`**: Listed in `Suggests`. When a user invokes `to_phyloseq()` or `to_tse()`, `requireNamespace()` verifies availability. If missing, clear instructions are returned on how to install via BiocManager (`BiocManager::install("phyloseq")`).

---

## 3. Data Conversion Layer (`R/bridges.R`)

### 3.1 Vegan Community Data Interoperability
- **`to_vegan(tb, assay = "counts")`**:
  - Extracts the requested assay matrix ($P$ taxa $\times$ $N$ samples) and transposes it to the `vegan` convention ($N$ samples in rows $\times$ $P$ taxa in columns).
  - Returns an S3 object of class `tidybiome_vegan` containing:
    - `$comm`: $N \times P$ community matrix with sample IDs as row names and taxon IDs as column names.
    - `$env`: Sample metadata tibble with `sample_id` as row names.
  - S3 methods: `print.tidybiome_vegan` providing a concise summary of sample/taxa counts and metadata variables.
- **`from_vegan(comm, env = NULL, tax_table = NULL, phy_tree = NULL)`**:
  - Accepts a community matrix (samples in rows, taxa in columns), transposes it into an assay matrix, validates dimensions against `env` and `tax_table`, and returns a valid `tidy_microbiome` object.
- **`as_tidybiome.vegan` / `as_tidybiome.tidybiome_vegan`**:
  - S3 method dispatch allowing seamless coercion.

### 3.2 `phyloseq` Interoperability
- **`to_phyloseq(tb)` & `as_phyloseq(tb)`**:
  - Converts `tidy_microbiome` to `phyloseq::phyloseq`.
  - Maps `assays$counts` to `phyloseq::otu_table(..., taxa_are_rows = TRUE)`.
  - Maps sample metadata to `phyloseq::sample_data()`.
  - Maps `tax_table` to `phyloseq::tax_table()`.
  - Maps `phy_tree` to `phyloseq::phy_tree()` if present.
- **`from_phyloseq(ps)` & `as_tidybiome.phyloseq(ps)`**:
  - Extracts `otu_table` (handling both `taxa_are_rows = TRUE` and `FALSE`), `sample_data`, `tax_table`, and `phy_tree`.
  - Reconstructs a clean `tidy_microbiome` S3 object.

### 3.3 `mia` / `TreeSummarizedExperiment` Interoperability
- **`to_tse(tb)` / `to_mia(tb)` / `as_mia(tb)`**:
  - Converts `tidy_microbiome` to `TreeSummarizedExperiment::TreeSummarizedExperiment`.
  - Transfers all assays in `attr(tb, "assays")` (`counts`, `rclr`, `relabundance`, etc.).
  - Transfers sample metadata to `colData`.
  - Transfers `tax_table` to `rowData`.
  - Transfers `phy_tree` to `rowTree` if present.
- **`from_tse(tse)` & `as_tidybiome.TreeSummarizedExperiment(tse)` / `as_tidybiome.SummarizedExperiment(se)`**:
  - Extracts primary assay (`counts` or first assay), metadata (`colData`), taxonomy (`rowData`), and phylogenetic tree (`rowTree`).
  - Transfers any additional secondary assays into `attr(tb, "assays")`.

### 3.4 Phylogenetic Tree Integration (`ape`)
- **`set_tree(tb, tree, prune = TRUE)`**:
  - Attaches an `ape::phylo` object to `attr(tb, "phy_tree")`.
  - Validates that tip labels match the taxon IDs. If `prune = TRUE`, prunes tree tips to match available taxa (or drops taxa missing from tree).
- **`get_tree(tb)`**:
  - Extracts `attr(tb, "phy_tree")`. Returns `NULL` if not set.

---

## 4. Tidy Analysis Wrappers (`R/connectors_vegan.R`)

### 4.1 Permutational Multivariate Analysis of Variance (`run_permanova`)
- **Signature**:
  ```r
  run_permanova(
    tb,
    formula,
    assay = "counts",
    method = "bray",
    permutations = 999,
    by = "terms",
    ...
  )
  ```
- **Behavior**:
  - Validates `formula` against variables in `tb`.
  - Converts `tb` to `vegan` format via `to_vegan(tb, assay = assay)`.
  - Executes `vegan::adonis2(formula = formula, data = env, method = method, permutations = permutations, by = by, ...)`.
  - Formats output into a clean, tidy tibble:
    `term`, `df`, `sum_sq`, `r2`, `f_stat`, `p_value`.
  - Attaches metadata attributes: `method`, `permutations`, `by`.

### 4.2 Multivariate Dispersion Testing (`run_betadisper`)
- **Signature**:
  ```r
  run_betadisper(
    tb,
    group,
    assay = "counts",
    method = "bray",
    permutations = 999,
    ...
  )
  ```
- **Behavior**:
  - Computes distance matrix $D$ on the requested assay using `vegan::vegdist()` or internal distance metric.
  - Calls `vegan::betadisper(d = D, group = group_vec)`.
  - Runs permutation test `vegan::permutest()` on the dispersion object.
  - Returns an S3 object of class `tidybiome_betadisper` with:
    - `$sample_distances`: Tibble containing `sample_id`, `group`, and `distance_to_centroid`.
    - `$group_summary`: Tibble containing mean distance and standard error per group.
    - `$test`: Tidy tibble containing F-statistic, degrees of freedom, and permutation $p$-value.
  - Pluggable into `ggplot2` for raincloud or boxplot visualization of dispersion across cohorts.

### 4.3 Non-Metric Multidimensional Scaling (`run_nmds`)
- **Signature**:
  ```r
  run_nmds(
    tb,
    assay = "counts",
    method = "bray",
    k = 2,
    trymax = 50,
    trace = 0,
    ...
  )
  ```
- **Behavior**:
  - Executes `vegan::metaMDS(comm, distance = method, k = k, trymax = trymax, trace = trace, ...)`.
  - Extracts ordination sample scores and attaches sample metadata.
  - Returns a tidy ordination structure compatible with `plot_ordination()`, containing stress value, 2D coordinates (`NMDS1`, `NMDS2`), and sample covariates.

---

## 5. Error Handling & Edge Cases
1. **Missing Third-Party Packages**:
   - For `vegan`: If not installed, error cleanly with message `"Package 'vegan' is required. Please run install.packages('vegan')."`.
   - For `phyloseq` / `TreeSummarizedExperiment`: Error with `"Package 'phyloseq' / 'TreeSummarizedExperiment' is required. Please install via BiocManager::install('...')."`.
2. **Dimension & Name Mismatches**:
   - `from_vegan()` automatically aligns row names with sample IDs and column names with taxon IDs. If row names or column names are missing, informative synthetic IDs are created (`Sample_1`, `Taxon_1`).
3. **Tree Tip Label Alignment**:
   - `set_tree()` checks whether `tree$tip.label` overlaps with `attr(tb, "tax_table")$taxon_id`. If overlap is partial, emits informative warning and prunes cleanly when `prune = TRUE`.

---

## 6. Verification Plan

### 6.1 Unit Tests (`tests/testthat/test-connectors.R`)
1. **Vegan conversion**: Test `to_vegan()` and `from_vegan()` round-trip consistency of count matrices and metadata.
2. **PERMANOVA wrapper**: Test `run_permanova(gut_microbiome, ~ treatment)` against direct `vegan::adonis2` output; confirm tidy tibble columns and $p$-value.
3. **Betadisper wrapper**: Test `run_betadisper(gut_microbiome, "treatment")`; confirm sample distance tibble and permutation test statistics.
4. **NMDS wrapper**: Test `run_nmds(gut_microbiome)`; confirm NMDS coordinates and stress value extraction.
5. **Tree integration**: Test `set_tree()` and `get_tree()` using an `ape::rtree(40)` random tree; verify tip label synchronization and pruning.
6. **Graceful fallbacks**: Test behavior of `to_phyloseq()` and `to_tse()` when packages are absent, and when mocked S4 objects are passed to `from_phyloseq()` / `from_tse()`.

### 6.2 Showcase Integration (`demo/tidybiome_showcase.R`)
- Add a dedicated showcase section demonstrating:
  - Export to `vegan` format via `to_vegan()`.
  - Execution of `run_permanova()` and `run_betadisper()`.
  - Execution of `run_nmds()` and visualization with `plot_ordination()`.
  - Attachment of an `ape` phylogenetic tree via `set_tree()`.
