#' Synthetic Benchmark Human Gut Microbiome Dataset
#'
#' @description
#' A realistic synthetic human gut microbiome benchmark dataset designed for testing
#' and demonstrating modern microbiome analysis workflows.
#' Contains 24 samples (12 Control vs. 12 Treated) and 40 amplicon sequence variants (ASVs)
#' across 4 major bacterial phyla (Bacteroidota, Firmicutes, Actinobacteriota, and Proteobacteria).
#'
#' @format A `tidy_microbiome` object with 24 samples and 40 taxa:
#' \describe{
#'   \item{sample_id}{Unique sample identifier (e.g. `Sample_01`).}
#'   \item{treatment}{Treatment group: `Control` or `Treated`.}
#'   \item{diet}{Dietary regimen: `Standard` or `HighFiber`.}
#'   \item{age}{Subject age in years.}
#'   \item{depth}{Total sequencing read count per sample.}
#' }
#'
#' @source Synthetically generated with negative binomial dispersion, compositional closure,
#'   and realistic zero inflation mimicking Illumina NovaSeq 16S rRNA amplicon sequencing.
#' @usage data(gut_microbiome)
"gut_microbiome"
