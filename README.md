# tidybiome: Tidyverse-Native Modern Microbiome Analysis Suite

[![R-CMD-check](https://img.shields.io/badge/R-4.1+-blue.svg)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**`tidybiome`** is a tidyverse-native, high-aesthetic R package designed as an up-to-date, simpler "Swiss Army knife" for downstream microbiome analysis. It adheres strictly to `tidyR` and `tidyverse` principles while incorporating cutting-edge methodologies from recent literature and **bioRxiv preprints (2024–2026)**.

---

## Why `tidybiome`?

Traditional microbiome workflows (e.g. `phyloseq` or Bioconductor's `mia` / `TreeSummarizedExperiment`) often require complex S4 slot gymnastics, produce dated visuals (e.g., 50-color rainbow legends, harsh black borders), and rely on legacy heuristics (such as arbitrary rarefaction or pseudocount-based CLR).

`tidybiome` delivers:
1. **True Tidy Ergonomics**: The core `tidy_microbiome` container inherits from `tbl_df`. Standard `dplyr` verbs (`filter`, `mutate`, `select`, `arrange`, `slice`) work natively on sample metadata while keeping count matrices and taxonomy tables synchronized.
2. **Updated SOTA & bioRxiv Methods**:
   * **Normalization**: Robust CLR (`rclr`, Martino et al.) without pseudocount distortion, alongside Coverage-based standardization and classic TSS.
   * **Alpha Diversity**: Unified **Hill Numbers Profile** ($q = 0, 1, 2$) placing Richness, Exponential Shannon, and Inverse Simpson on the same intuitive scale of *effective number of species*.
   * **Beta Diversity**: **Robust Aitchison Distance** (`raitchison`) and Jensen-Shannon Divergence (`jsd`).
   * **Dimensionality Reduction**: **Robust PCA (RPCA)** singular value decomposition with true taxon loading vectors for compositional biplots.
   * **Differential Abundance**: **Multi-Engine Consensus Engine** (bioRxiv 2024–2025: *ConsensusMetaDA*, *dar*, *LinDA*) combining LinDA, CLR-linear models, and Wilcoxon tests via the Cauchy combination test and agreement scoring.
   * **Core Microbiome**: Data-driven inflection curve analysis to discover natural community boundaries without arbitrary cutoffs.
3. **Publication-Grade Visuals**: Intelligent top-$N$ taxonomic pooling with soft-gray `"Other"` base bars, colorblind-safe palettes (Okabe-Ito, Nature), confidence ellipses, and `theme_tidybiome()`.

---

## Comparison: `tidybiome` vs Classic Tools

| Analysis Area | Classic `phyloseq` / `mia` | `tidybiome` Modern SOTA Method |
| :--- | :--- | :--- |
| **Object System** | S4 (`phyloseq`, `TreeSummarizedExperiment`) | S3 `tidy_microbiome` inheriting from `tbl_df` (native `dplyr`) |
| **Normalization** | Rarefaction, pseudocount CLR | **Robust CLR (rCLR)** & **Coverage Standardization** |
| **Alpha Diversity** | Disconnected Shannon, Simpson, Chao1 | **Unified Hill Numbers Profile ($q=0, 1, 2$)** |
| **Beta Diversity** | Bray-Curtis, Jaccard | **Robust Aitchison Distance** & **JSD** |
| **Ordination** | Basic PCoA / NMDS | **Robust PCA (RPCA / DEICODE-style biplot)** & corrected PCoA |
| **Differential Abundance** | Single-tool (often high false discovery rate) | **Multi-Engine Consensus DA** (LinDA, CLR-linear, Wilcoxon + Cauchy test) |
| **Visual Aesthetics** | Cluttered 50-color bars, no auto-pooling | **Intelligent top-N pooling**, colorblind-safe palettes, `theme_tidybiome` |

---

## Installation

```r
# Install from source repository
# install.packages("devtools")
devtools::install_local(".")
```

---

## Quickstart

```r
library(tidybiome)
library(dplyr)
library(ggplot2)

# Load built-in benchmark dataset
data("gut_microbiome", package = "tidybiome")
gut_microbiome

# 1. Tidy manipulation
clean_tb <- gut_microbiome |>
  filter(depth > 10000) |>
  mutate(log_depth = log10(depth))

# 2. Robust CLR transformation (no pseudocount artifacts)
clean_tb <- transform_abundance(clean_tb, method = "rclr")

# 3. Unified Hill diversity profile (q=0, 1, 2)
clean_tb <- calc_alpha_diversity(clean_tb, metrics = "hill")

# 4. Robust Aitchison distance & RPCA biplot
clean_tb <- calc_beta_diversity(clean_tb, metric = "raitchison")
clean_tb <- calc_ordination(clean_tb, method = "rpca")

# 5. Multi-engine consensus differential abundance
da_res <- calc_differential_abundance(
  clean_tb,
  group = "treatment",
  methods = c("consensus", "linda", "clr_linear", "wilcoxon")
)

# 6. Aesthetic visualizations
plot_composition(clean_tb, rank = "Phylum", top_n = 4, group_by = "treatment")
plot_ordination(clean_tb, method = "rpca", color = "treatment", ellipse = TRUE, biplot = TRUE)
plot_alpha(clean_tb, metric = "hill_1", x = "treatment", test = TRUE)
plot_da_volcano(da_res)
```

---

## Running Demo & Tests

```bash
# Run the complete end-to-end showcase
Rscript demo/tidybiome_showcase.R

# Run the testthat test suite (81 assertions)
Rscript -e 'testthat::test_dir("tests/testthat")'
```

---

## License

MIT License (c) 2026 tidybiome authors.
