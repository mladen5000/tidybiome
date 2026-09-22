# Package index

## Pipeline Ingestion & Import

Ingest outputs from DADA2, QIIME 2, and MetaPhlAn pipelines.

- [`import_dada2()`](https://tidybiome.org/reference/import_dada2.md) :
  Import DADA2 Pipeline Outputs into tidy_microbiome
- [`import_qiime2()`](https://tidybiome.org/reference/import_qiime2.md)
  : Import QIIME 2 Artifacts or TSV Tables
- [`import_metaphlan()`](https://tidybiome.org/reference/import_metaphlan.md)
  : Import MetaPhlAn Taxonomic Profiles

## Container & Data Wrangling

Creating, subsetting, and validating tidy_microbiome objects.

- [`tidy_microbiome()`](https://tidybiome.org/reference/tidy_microbiome.md)
  : Create a tidy_microbiome Object
- [`summary(`*`<tidy_microbiome>`*`)`](https://tidybiome.org/reference/summary.tidy_microbiome.md)
  : Summary of a tidy_microbiome Object
- [`filter_taxa()`](https://tidybiome.org/reference/filter_taxa.md) :
  Filter Taxa Based on Taxonomy Table Attributes
- [`validate_tidy_microbiome()`](https://tidybiome.org/reference/validate_tidy_microbiome.md)
  : Validate Integrity of a tidy_microbiome Object
- [`tidy_abundance()`](https://tidybiome.org/reference/tidy_abundance.md)
  : Extract Abundance Data in Tidy Format
- [`tidy_taxa()`](https://tidybiome.org/reference/tidy_taxa.md) :
  Extract Taxonomy Annotations as a Tidy Tibble
- [`assay()`](https://tidybiome.org/reference/assay.md) : Extract Assay
  Matrix from tidy_microbiome
- [`tax_table()`](https://tidybiome.org/reference/tax_table.md) :
  Extract Taxonomy Table from tidy_microbiome
- [`aggregate_taxa()`](https://tidybiome.org/reference/aggregate_taxa.md)
  : Aggregate Taxa to a Higher Taxonomic Rank

## Transformations & Normalization

Compositional and sequencing depth standardizations.

- [`transform_abundance()`](https://tidybiome.org/reference/transform_abundance.md)
  : Transform Abundance Data Using Compositional or Classical
  Normalization
- [`calc_gmpr_size_factors()`](https://tidybiome.org/reference/calc_gmpr_size_factors.md)
  : Calculate GMPR Size Factors

## Sample Quality Control & Diagnostics

Per-sample library depth, feature richness, sparsity, and diagnostic
plots.

- [`calc_qc_metrics()`](https://tidybiome.org/reference/calc_qc_metrics.md)
  : Calculate Sample-Level Quality Control (QC) Metrics
- [`plot_qc()`](https://tidybiome.org/reference/plot_qc.md) : Plot
  Sample Quality Control Diagnostics

## Diversity & Distance Metrics

Unified alpha diversity profiles and beta diversity dissimilarities.

- [`calc_alpha_diversity()`](https://tidybiome.org/reference/calc_alpha_diversity.md)
  : Calculate Alpha Diversity and Simplex Compositional Metrics
- [`calc_beta_diversity()`](https://tidybiome.org/reference/calc_beta_diversity.md)
  : Calculate Beta Diversity Dissimilarity and Distance Matrices
- [`calc_unifrac()`](https://tidybiome.org/reference/calc_unifrac.md) :
  Calculate Unweighted or Weighted UniFrac Distance
- [`get_distance()`](https://tidybiome.org/reference/get_distance.md) :
  Extract Distance Matrix from tidy_microbiome
- [`dth_test()`](https://tidybiome.org/reference/dth_test.md) :
  Distance-based Test for Homogeneity (DTH)

## Dimensionality Reduction & Ordination

RPCA biplots, contrastive PCA, NMDS, and db-RDA.

- [`calc_ordination()`](https://tidybiome.org/reference/calc_ordination.md)
  : Calculate Dimensionality Reduction and Ordination
- [`get_ordination()`](https://tidybiome.org/reference/get_ordination.md)
  : Extract Ordination Results from tidy_microbiome
- [`run_dbrda()`](https://tidybiome.org/reference/run_dbrda.md) :
  Distance-Based Redundancy Analysis via vegan::dbrda
- [`run_nmds()`](https://tidybiome.org/reference/run_nmds.md) : Run
  Non-Metric Multidimensional Scaling via vegan::metaMDS

## Differential Abundance

Consensus and multi-engine biomarker discovery with formulas and
covariates.

- [`calc_differential_abundance()`](https://tidybiome.org/reference/calc_differential_abundance.md)
  : Multi-Engine Consensus Differential Abundance Analysis

## Ecological Screening, Networks & Core Microbiome

Prevalence, dominance, divergence, cross-associations, and co-occurrence
graphs.

- [`calc_prevalence()`](https://tidybiome.org/reference/calc_prevalence.md)
  : Calculate Prevalence and Abundance Metrics Across Taxa
- [`filter_prevalent()`](https://tidybiome.org/reference/filter_prevalent.md)
  : Filter Microbiome by Prevalence and Abundance Thresholds
- [`calc_dominant()`](https://tidybiome.org/reference/calc_dominant.md)
  : Identify Dominant Taxa per Sample and Cohort-Wide
- [`calc_divergence()`](https://tidybiome.org/reference/calc_divergence.md)
  : Calculate Ecological Divergence Relative to a Reference State
- [`calc_cross_association()`](https://tidybiome.org/reference/calc_cross_association.md)
  : Calculate Cross-Associations Between Taxa and Metadata Variables
- [`calc_network()`](https://tidybiome.org/reference/calc_network.md) :
  Calculate Microbial Co-Occurrence Network
- [`calc_core_microbiome()`](https://tidybiome.org/reference/calc_core_microbiome.md)
  : Identify Core Microbiome Members

## Ecosystem Interoperability Bridges

Connectors to vegan, phyloseq, mia, S7, and ape.

- [`to_vegan()`](https://tidybiome.org/reference/to_vegan.md) : Convert
  tidy_microbiome to vegan Community and Environmental Format
- [`from_vegan()`](https://tidybiome.org/reference/from_vegan.md) :
  Construct tidy_microbiome from vegan Community Format
- [`run_permanova()`](https://tidybiome.org/reference/run_permanova.md)
  : Run PERMANOVA via vegan::adonis2
- [`run_betadisper()`](https://tidybiome.org/reference/run_betadisper.md)
  : Test for Homogeneity of Multivariate Dispersions via
  vegan::betadisper
- [`to_phyloseq()`](https://tidybiome.org/reference/to_phyloseq.md)
  [`as_phyloseq()`](https://tidybiome.org/reference/to_phyloseq.md) :
  Convert tidy_microbiome to phyloseq Object
- [`to_tse()`](https://tidybiome.org/reference/to_tse.md)
  [`as_mia()`](https://tidybiome.org/reference/to_tse.md)
  [`to_mia()`](https://tidybiome.org/reference/to_tse.md) : Convert
  tidy_microbiome to TreeSummarizedExperiment Object
- [`to_s7()`](https://tidybiome.org/reference/to_s7.md) : Convert
  tidy_microbiome to S7 Formal Object
- [`from_s7()`](https://tidybiome.org/reference/from_s7.md) : Convert S7
  Formal Object to tidy_microbiome S3 Tibble
- [`as_tidybiome()`](https://tidybiome.org/reference/as_tidybiome.md) :
  Coerce Objects to tidy_microbiome
- [`set_tree()`](https://tidybiome.org/reference/set_tree.md) : Attach
  and Validate Phylogenetic Tree to tidy_microbiome
- [`get_tree()`](https://tidybiome.org/reference/get_tree.md) : Get
  Phylogenetic Tree from tidy_microbiome

## Publication-Grade Visualizations & Dashboards

High-aesthetic ggplot2 extensions and standalone HTML dashboards.

- [`plot_composition()`](https://tidybiome.org/reference/plot_composition.md)
  : Plot Taxonomic Composition with Intelligent Top-N Grouping
- [`plot_heatmap()`](https://tidybiome.org/reference/plot_heatmap.md) :
  Plot Compositional Microbiome Heatmap
- [`plot_ordination()`](https://tidybiome.org/reference/plot_ordination.md)
  : Plot Microbiome Ordination and Compositional Biplots
- [`plot_alpha()`](https://tidybiome.org/reference/plot_alpha.md) : Plot
  Alpha Diversity with Statistical Comparisons
- [`plot_da_volcano()`](https://tidybiome.org/reference/plot_da_volcano.md)
  : Plot Consensus Differential Abundance Volcano Plot
- [`plot_dbrda()`](https://tidybiome.org/reference/plot_dbrda.md) : Plot
  Distance-Based Redundancy Analysis (db-RDA) Biplot
- [`plot_network()`](https://tidybiome.org/reference/plot_network.md) :
  Plot Microbial Co-Occurrence Network
- [`plot_core()`](https://tidybiome.org/reference/plot_core.md) : Plot
  Core Microbiome Landscape
- [`plot_qc()`](https://tidybiome.org/reference/plot_qc.md) : Plot
  Sample Quality Control Diagnostics
- [`report_tidybiome()`](https://tidybiome.org/reference/report_tidybiome.md)
  : Generate Standalone Diagnostic HTML Summary Report
- [`theme_tidybiome()`](https://tidybiome.org/reference/theme_tidybiome.md)
  : Modern Minimalist Publication Theme for tidybiome
- [`tidybiome_pal()`](https://tidybiome.org/reference/tidybiome_pal.md)
  [`scale_color_tidybiome()`](https://tidybiome.org/reference/tidybiome_pal.md)
  [`scale_fill_tidybiome()`](https://tidybiome.org/reference/tidybiome_pal.md)
  : Curated Color Palettes for tidybiome

## Example Datasets

Built-in demonstration datasets.

- [`gut_microbiome`](https://tidybiome.org/reference/gut_microbiome.md)
  : Synthetic Benchmark Human Gut Microbiome Dataset
