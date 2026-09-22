# tidybiome: Tidyverse-Native Modern Microbiome Analysis Suite

[![R-CMD-check](https://img.shields.io/badge/R-4.1+-blue.svg)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests: 239 passing](https://img.shields.io/badge/tests-239%20passed-success.svg)](https://github.com/tidybiome/tidybiome)

**`tidybiome`** is a tidyverse-native, publication-aesthetic R package designed as an up-to-date "Swiss Army knife" for downstream microbiome analysis. It adheres strictly to modern `tidyR` and `tidyverse` principles while incorporating cutting-edge methodologies from recent literature and **bioRxiv preprints (2024–2026)**.

---

## Why `tidybiome`?

Traditional microbiome workflows (e.g. `phyloseq` or Bioconductor's `mia` / `TreeSummarizedExperiment`) often require complex S4 slot accessors, produce dated visuals (e.g., 50-color rainbow legends, harsh black borders), and rely on legacy heuristics (such as arbitrary rarefaction or pseudocount-based CLR).

`tidybiome` delivers:
1. **True Tidy Ergonomics**: The core `tidy_microbiome` container inherits from `tbl_df`. Standard `dplyr` verbs (`filter`, `mutate`, `select`, `arrange`, `slice`) work natively on sample metadata while keeping count matrices, taxonomy tables, and phylogenetic trees synchronized.
2. **Updated SOTA & bioRxiv Methods**:
   * **Normalization**: Robust CLR (`rclr`, Martino et al.) without pseudocount distortion, alongside Coverage-based standardization and classic TSS.
   * **Alpha Diversity**: Unified **Hill Numbers Profile** ($q = 0, 1, 2$) placing Richness, Exponential Shannon, and Inverse Simpson on the same intuitive scale of *effective number of species*, plus **Simplex Compositional Diversity** (bioRxiv 2026).
   * **Beta Diversity**: **Robust Aitchison Distance** (`raitchison`), **Tree-Wasserstein Distance** (Earth Mover's Distance across the taxonomic hierarchy), **Unweighted & Weighted UniFrac** (`calc_unifrac()`), and the **DTH Test for Homogeneity** (bioRxiv 2025).
   * **Dimensionality Reduction**: **Robust PCA (RPCA)** singular value decomposition with true taxon loading vectors for compositional biplots, plus **Contrastive PCA (cPCA)** isolating condition-specific dysbiosis.
   * **Differential Abundance**: **Multi-Engine Consensus Engine** (bioRxiv 2024–2025: *ConsensusMetaDA*, *dar*, *LinDA*) with formula specification (`formula = ~ treatment + age + batch`), covariate adjustments, continuous predictors, and **CAFT** zero-cell separation modeling.
   * **Ecological Modeling**: Distance-based redundancy analysis (`run_dbrda()`, `plot_dbrda()`) and microbial co-occurrence networks (`calc_network()`, `plot_network()`).
   * **Core Microbiome**: Data-driven inflection curve analysis to discover natural community boundaries without arbitrary cutoffs.
3. **Publication-Grade Visuals**: Intelligent top-$N$ taxonomic pooling with soft-gray `"Other"` base bars, colorblind-safe palettes (Okabe-Ito, Nature), confidence ellipses, and `theme_tidybiome()`.

---

## Architecture & Container Design

A `tidy_microbiome` object is a tibble on the outside (holding sample metadata) with internal multi-assay matrices, taxonomic lineage tables, distance caches, ordination embeddings, and phylogenetic trees kept strictly in lockstep:

```
┌─────────────────────────────────────────────────────────────────┐
│                     tidy_microbiome container                   │
│                                                                 │
│   Sample Metadata (tibble)                                      │
│   ├── sample_id (primary key)                                   │
│   └── host / clinical / environmental covariates                │
│                                                                 │
│   Assay Matrices (taxa x samples)                               │
│   ├── counts (raw sequence reads)                               │
│   ├── rclr (robust centered log-ratio, zero-safe)               │
│   └── tss (relative proportions)                                │
│                                                                 │
│   Taxonomy Table (taxa x ranks)                                 │
│   └── taxon_id, Kingdom, Phylum, Class, Order, Family, Genus    │
│                                                                 │
│   Phylogenetic Tree (ape::phylo)                                │
│   └── Tips mapped 1:1 to taxon_id                               │
│                                                                 │
│   Distance Cache & Ordination Embeddings                        │
│   └── raitchison, bray, unifrac, rpca, cpca, dbrda              │
└─────────────────────────────────────────────────────────────────┘
```

---

## Comparison: `tidybiome` vs `phyloseq`, `mia`, and `vegan`

| Analysis Area | Standard in `phyloseq` / `mia` / `vegan` | `tidybiome` Capabilities |
| :--- | :--- | :--- |
| **Object System** | S4 (`phyloseq`, `TreeSummarizedExperiment`) | S3 `tidy_microbiome` inheriting from `tbl_df` (native `dplyr` pipelines) |
| **Data Integrity** | Manual subsetting across slots | Automatic S3 bracket subsetting (`tb[i, j]`), `filter_taxa()`, and `validate_tidy_microbiome()` |
| **Differential Abundance** | DESeq2 (`phyloseq`), ALDEx2/ANCOMBC (`mia`), None (`vegan`) | **CAFT (bioRxiv Dec 2025)** zero-cell separation model; **Multi-Engine Consensus DA** (*ConsensusMetaDA* / *dar* 2024–2025) with full formula covariate regression (`~ treatment + age + batch`) |
| **Beta Diversity** | Bray-Curtis, Jaccard (`vegan`/`phyloseq`), standard Aitchison | **Vectorized Bray/Jaccard/JSD**; **Unweighted & Weighted UniFrac**; **Tree-Wasserstein Distance**; **DTH Test for Homogeneity (bioRxiv 2025)** |
| **Alpha Diversity** | Shannon, Simpson, Chao1 on raw/rarefied counts | **Simplex Compositional Diversity (bioRxiv 2026)**; Unified **Hill Numbers Profile ($q=0, 1, 2$)** on the effective species scale |
| **Ordination & Projection** | PCoA, NMDS, CCA, RDA (`vegan`/`phyloseq`) | **Contrastive PCA (cPCA)** isolating dysbiosis against healthy background; **Robust PCA (RPCA)** biplots; **db-RDA** constrained ordination |
| **Microbial Networks** | External packages (`SpiecEasi`, `igraph`) | Built-in `calc_network()` and `plot_network()` with correlation filtering, FDR correction, and degree centrality |
| **Core Microbiome** | Static cutoffs (e.g. 80% prev, 0.1% abund) | **Data-driven inflection curve analysis**: Discovers natural mathematical breakpoints in the prevalence landscape |
| **Visual Aesthetics** | Cluttered 50-color bars, harsh borders | **Intelligent top-N pooling** into soft-gray `"Other"`, colorblind-safe palettes, `theme_tidybiome` |

---

## Installation

```r
# Install from local source repository:
devtools::install_local(".")

# Or install dependencies first if needed:
install.packages(c("tibble", "dplyr", "rlang", "ggplot2", "vegan", "ape"))
```

---

## Quickstart

```r
library(tidybiome)
library(dplyr)
library(ggplot2)

# Load built-in benchmark dataset (24 samples x 40 taxa)
data("gut_microbiome", package = "tidybiome")

# 1. Tidy manipulation and filtering
clean_tb <- gut_microbiome |>
  filter(depth > 10000) |>
  mutate(log_depth = log10(depth))

# 2. Robust CLR transformation (zero-safe, no pseudocount artifacts)
clean_tb <- transform_abundance(clean_tb, method = "rclr")

# 3. Unified Hill diversity profile (q=0 richness, q=1 exponential Shannon, q=2 inverse Simpson)
clean_tb <- calc_alpha_diversity(clean_tb, metrics = "hill")

# 4. Beta diversity & RPCA biplot
clean_tb <- calc_beta_diversity(clean_tb, metric = "raitchison")
clean_tb <- calc_ordination(clean_tb, method = "rpca")

# 5. Formula-based Consensus Differential Abundance with Covariate Adjustment
da_res <- calc_differential_abundance(
  clean_tb,
  formula = ~ treatment + age,
  contrast = "treatment",
  methods = c("consensus", "linda", "clr_linear", "wilcoxon")
)

# Inspect top significant biomarkers
da_res |>
  filter(is_significant) |>
  select(taxon_id, Phylum, Genus, log2fc, padj_consensus, agreement_score) |>
  head(6)

# 6. Aesthetic Visualizations
plot_composition(clean_tb, rank = "Phylum", top_n = 4, group_by = "treatment")
plot_ordination(clean_tb, method = "rpca", color = "treatment", ellipse = TRUE, biplot = TRUE)
plot_alpha(clean_tb, metric = "hill_1", x = "treatment", test = TRUE)
plot_da_volcano(da_res)

# 7. Distance-Based Redundancy Analysis (db-RDA)
dbrda_res <- run_dbrda(clean_tb, ~ treatment + age, distance = "bray")
plot_dbrda(dbrda_res, color = "treatment")

# 8. Microbial Co-Occurrence Networks
net_res <- calc_network(clean_tb, method = "spearman", min_prevalence = 0.3, r_cutoff = 0.3, p_cutoff = 0.05)
plot_network(net_res)

# 9. Compositional Heatmap with Hierarchical Clustering
plot_heatmap(clean_tb, rank = "Genus", top_n = 15, annotation_col = "treatment")

# 10. Phylogenetic Integration & UniFrac Distances
tree <- ape::rtree(40, tip.label = attr(clean_tb, "tax_table")$taxon_id)
clean_tb <- set_tree(clean_tb, tree)
clean_tb <- calc_beta_diversity(clean_tb, metric = "unifrac")
clean_tb <- calc_beta_diversity(clean_tb, metric = "wunifrac")

# 11. Container Integrity Validation & HTML Report Dashboard
validate_tidy_microbiome(clean_tb, verbose = TRUE)
report_tidybiome(clean_tb, output = "report.html", browse = FALSE)
```

---

## Complete API Reference

### Pipeline Ingestion & Import
| Function | Description | Key Arguments |
| :--- | :--- | :--- |
| `import_dada2()` | Imports DADA2 sequence table matrix/RDS and taxonomy. | `seqtab, taxa, sample_metadata, clean_names = TRUE` |
| `import_qiime2()` | Imports QIIME 2 exported feature table, taxonomy TSV, and Newick tree. | `feature_table, taxonomy, sample_metadata, tree` |
| `import_metaphlan()` | Parses MetaPhlAn merged abundance profile TSV with hierarchical clades. | `file, rank = "Species", sample_metadata` |

### Container & Data Wrangling
| Function | Description | Key Arguments |
| :--- | :--- | :--- |
| `tidy_microbiome()` | Constructs a new synchronized `tidy_microbiome` container. | `counts, sample_data, tax_table, phy_tree` |
| `[.tidy_microbiome` | S3 bracket subsetting (`tb[i, j]`) preserving assay column and tree tip alignment. | `x[i, j]` |
| `filter_taxa()` | Filters taxa based on taxonomic expressions, pruning assays and tree tips. | `tb, ...` (e.g. `Phylum == "Bacteroidota"`) |
| `validate_tidy_microbiome()` | Verifies 1:1 integrity between sample metadata, assays, taxonomy, and tree tips. | `tb, verbose = TRUE` |
| `tidy_abundance()` | Extracts long or wide abundance tables joined with metadata and taxonomy. | `tb, assay = "counts", long = TRUE` |
| `assay()` | Extracts or sets specific assay matrices (`counts`, `rclr`, `tss`, `gmpr`). | `tb, name = "counts"` |
| `tax_table()` | Extracts the taxonomic lineage table. | `tb` |
| `aggregate_taxa()` | Agglomerates abundances to a specific taxonomic rank, preserving trees. | `tb, rank = "Genus"` |

### Normalization & Transformations
| Function | Description | Key Arguments |
| :--- | :--- | :--- |
| `transform_abundance()` | Modern compositional and coverage transformations. | `tb, method = c("rclr", "gmpr", "hellinger", "tss", "coverage")` |
| `calc_gmpr_size_factors()` | Calculates sample size factors via Geometric Mean of Pairwise Ratios. | `mat, min_overlap = 2` |

### Diversity & Distance Metrics
| Function | Description | Key Arguments |
| :--- | :--- | :--- |
| `calc_alpha_diversity()` | Calculates Hill profiles ($q=0, 1, 2$), Shannon, Simpson, and Simplex variation. | `tb, metrics = c("hill", "shannon", "simplex_variation")` |
| `calc_beta_diversity()` | Calculates distance matrices across samples. | `tb, metric = c("raitchison", "bray", "jaccard", "jsd", "unifrac", "wunifrac", "tree_wasserstein")` |
| `calc_unifrac()` | Computes vectorized Unweighted and Weighted normalized UniFrac distances. | `tb, weighted = FALSE, normalized = TRUE` |
| `dth_test()` | Distance-based Test for Homogeneity permutation test (bioRxiv 2025). | `tb, group, metric = "tree_wasserstein", n_perm = 999` |

### Dimensionality Reduction & Ordination
| Function | Description | Key Arguments |
| :--- | :--- | :--- |
| `calc_ordination()` | Dimensionality reduction via RPCA, cPCA, or classical PCoA. | `tb, method = c("rpca", "cpca", "pcoa")` |
| `run_nmds()` | Non-metric multidimensional scaling via `vegan::metaMDS`. | `tb, distance = "bray", k = 2` |
| `run_dbrda()` | Distance-based redundancy analysis with biplot environmental arrows. | `tb, formula = ~ treatment + age, distance = "bray"` |

### Differential Abundance & Biomarker Discovery
| Function | Description | Key Arguments |
| :--- | :--- | :--- |
| `calc_differential_abundance()` | Consensus and multi-engine differential abundance with formula & covariates. | `tb, formula = ~ treatment + age, methods = c("consensus", "caft", "linda", "clr_linear", "wilcoxon")` |

### Ecological & Community Screening
| Function | Description | Key Arguments |
| :--- | :--- | :--- |
| `calc_prevalence()` | Calculates taxon detection prevalence and mean abundance. | `tb, detection_limit = 0` |
| `filter_prevalent()` | Filters community to taxa meeting detection prevalence thresholds. | `tb, min_prevalence = 0.1` |
| `calc_dominant()` | Identifies dominant taxon per sample. | `tb, add_to_metadata = TRUE` |
| `calc_divergence()` | Computes sample divergence from a reference group or baseline. | `tb, reference = list(treatment = "Control")` |
| `calc_cross_association()` | Computes pairwise correlations between taxa and continuous metadata variables. | `tb, variables = c("age", "depth"), method = "spearman"` |
| `calc_network()` | Constructs microbial co-occurrence networks with FDR filtering and degree centrality. | `tb, min_prevalence = 0.2, r_cutoff = 0.3, p_cutoff = 0.05` |

### Ecosystem Connectors & Interoperability
| Tool | Function | Description |
| :--- | :--- | :--- |
| **`vegan`** | `to_vegan(tb)`, `from_vegan(comm, env)` | Bidirectional bridge with vegan community matrices. |
| **`vegan`** | `run_permanova(tb, formula)` | Direct PERMANOVA (`vegan::adonis2`) returning broom-style tibbles. |
| **`vegan`** | `run_betadisper(tb, group)` | Multivariate dispersion permutation test. |
| **`phyloseq`** | `to_phyloseq(tb)`, `as_phyloseq(tb)`, `as_tidybiome(ps)` | Bidirectional bridge with `phyloseq::phyloseq` S4 objects. |
| **`mia` / `TreeSE`** | `to_tse(tb)`, `to_mia(tb)`, `as_mia(tb)`, `as_tidybiome(tse)` | Bidirectional bridge with `TreeSummarizedExperiment` S4 objects. |
| **`ape`** | `set_tree(tb, tree)`, `get_tree(tb)` | Attaches and retrieves phylogenetic trees. |

### Aesthetic Visualizations & Dashboards
| Function | Description | Key Features |
| :--- | :--- | :--- |
| `plot_composition()` | Relative abundance stacked bar chart. | Intelligent top-$N$ pooling, soft-gray `"Other"` base bar, Okabe-Ito palette. |
| `plot_heatmap()` | Compositional microbiome heatmap. | Hierarchical sample & taxon clustering, data scaling (log10, rclr, relative), metadata grouping. |
| `plot_ordination()` | Ordination scatter plot. | 95% confidence ellipses, top taxon biplot loading vectors. |
| `plot_alpha()` | Alpha diversity box/violin plot. | Jittered points, automated Wilcoxon or ANOVA significance brackets. |
| `plot_da_volcano()` | Volcano plot for differential abundance. | Dual thresholds, top biomarker labels, agreement score sizing. |
| `plot_dbrda()` | db-RDA constrained ordination plot. | Sample centroids and environmental constraint vector arrows. |
| `plot_network()` | Microbial co-occurrence network plot. | Positive (blue) vs negative (coral) edges, node sizes by degree. |
| `report_tidybiome()` | Standalone HTML diagnostic cohort dashboard. | Quality control cards, depth summary, metadata dictionary, zero external dependencies. |
| `theme_tidybiome()` | Publication-ready ggplot2 theme. | Clean minimalist typography, light borders, subtle gridlines. |

---

## Testing & Quality Assurance

The package includes a comprehensive unit test suite covering container integrity, distance calculation vectorization, formula modeling, and ecosystem bridges:

```bash
# Run the test suite:
Rscript -e 'testthat::test_dir("tests/testthat")'
# [ FAIL 0 | WARN 0 | SKIP 0 | PASS 239 ]

# Run the 12-step end-to-end showcase:
Rscript demo/tidybiome_showcase.R
```

In **RStudio**:
- Press `Cmd + Shift + T` (Mac) or `Ctrl + Shift + T` (Windows/Linux) to run the full test suite in the Build pane.
- Press `Cmd + Shift + L` (`devtools::load_all()`) for instant in-memory reloading during development.
- Click the editor gutter on any line to set visual breakpoints.

---

## References

1. **Martino, C. et al. (2019)**. *A Novel Sparse Compositional Metric Alleviates Problems in Biomarker Discovery*. **mSystems**, 4(1), e00016-19.
2. **Chao, A. et al. (2014)**. *Rarefaction and extrapolation with Hill numbers: a framework for analyzing species diversity*. **Ecological Monographs**, 84(1), 45-67.
3. **Zhou, H. et al. (2022)**. *LinDA: linear models for differential abundance analysis of microbiome compositional data*. **Genome Biology**, 23, 95.
4. **Lozupone, C. & Knight, R. (2005)**. *UniFrac: a new phylogenetic metric for comparing microbial communities*. **Applied and Environmental Microbiology**, 71(12), 8228-8235.
5. **Legendre, P. & Anderson, M. J. (1999)**. *Distance-based redundancy analysis: testing multispecies responses in multifactorial ecological experiments*. **Ecological Monographs**, 69(1), 1-24.
6. **ConsensusMetaDA & CAFT (bioRxiv 2024–2025)**. *Consensus differential abundance testing and zero-cell likelihood ratio modeling for high-dimensional compositional microbiome surveys*.

---

## License

MIT © tidybiome authors.
