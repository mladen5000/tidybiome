# Synthetic Benchmark Human Gut Microbiome Dataset

A realistic synthetic human gut microbiome benchmark dataset designed
for testing and demonstrating modern microbiome analysis workflows.
Contains 24 samples (12 Control vs. 12 Treated) and 40 amplicon sequence
variants (ASVs) across 4 major bacterial phyla (Bacteroidota,
Firmicutes, Actinobacteriota, and Proteobacteria).

## Usage

``` r
data(gut_microbiome)
```

## Format

A `tidy_microbiome` object with 24 samples and 40 taxa:

- sample_id:

  Unique sample identifier (e.g. `Sample_01`).

- treatment:

  Treatment group: `Control` or `Treated`.

- diet:

  Dietary regimen: `Standard` or `HighFiber`.

- age:

  Subject age in years.

- depth:

  Total sequencing read count per sample.

## Source

Synthetically generated with negative binomial dispersion, compositional
closure, and realistic zero inflation mimicking Illumina NovaSeq 16S
rRNA amplicon sequencing.
