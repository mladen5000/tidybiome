#!/usr/bin/env Rscript

# tidybiome: Modern, Tidy Microbiome Analysis Showcase
# ----------------------------------------------------
message(">>> Loading tidybiome and dependencies...")
library(tidybiome)
library(dplyr)
library(ggplot2)

# 1. Load benchmark dataset
message("\n--- Step 1: Inspecting tidy_microbiome container ---")
data("gut_microbiome", package = "tidybiome")
print(gut_microbiome)

# 2. Tidy data manipulation
message("\n--- Step 2: Tidyverse dplyr verbs in action ---")
gut_clean <- gut_microbiome |>
  filter(depth > 10000) |>
  mutate(log_depth = log10(depth))

cat(sprintf("Filtered to %d high-depth samples.\n", nrow(gut_clean)))

# 3. Compositional Normalization (bioRxiv/CoDA: Robust CLR)
message("\n--- Step 3: Modern Compositional Transformation (rCLR) ---")
gut_clean <- transform_abundance(gut_clean, method = "rclr")
cat("Robust CLR assay added without zero-distortion or pseudocount artifacts.\n")

# 4. Alpha Diversity (Unified Hill Numbers Profile)
message("\n--- Step 4: Alpha Diversity (Hill Profile q=0, 1, 2) ---")
gut_clean <- calc_alpha_diversity(gut_clean, metrics = c("hill", "shannon", "chao1"))
cat("Hill numbers calculated:\n")
print(gut_clean |> select(sample_id, treatment, hill_0, hill_1, hill_2, shannon) |> head(4))

# 5. Beta Diversity (Robust Aitchison Distance)
message("\n--- Step 5: Beta Diversity (Robust Aitchison) ---")
gut_clean <- calc_beta_diversity(gut_clean, metric = "raitchison")
d_mat <- get_distance(gut_clean, "raitchison")
cat(sprintf("Calculated %dx%d Robust Aitchison distance matrix.\n", attr(d_mat, "Size"), attr(d_mat, "Size")))

# 6. Dimensionality Reduction (RPCA Biplot)
message("\n--- Step 6: Dimensionality Reduction (RPCA Biplot) ---")
gut_clean <- calc_ordination(gut_clean, method = "rpca")
ord <- get_ordination(gut_clean, "rpca")
cat(sprintf("RPCA Variance Explained: PC1 = %.1f%%, PC2 = %.1f%%\n",
            ord$variance_explained[1] * 100, ord$variance_explained[2] * 100))

# 7. SOTA Multi-Engine Consensus Differential Abundance
message("\n--- Step 7: Multi-Engine Consensus Differential Abundance ---")
da_results <- calc_differential_abundance(
  gut_clean,
  group = "treatment",
  methods = c("consensus", "linda", "clr_linear", "wilcoxon"),
  fdr_cutoff = 0.05
)

cat("Top Significant Biomarkers (Consensus):\n")
print(da_results |>
        filter(is_significant) |>
        select(taxon_id, Phylum, Genus, log2fc, padj_consensus, agreement_score) |>
        head(6))

# 8. Core Microbiome Discovery
message("\n--- Step 8: Core Microbiome Discovery ---")
core_results <- calc_core_microbiome(gut_clean, method = "threshold", min_prevalence = 0.75)
cat(sprintf("Identified %d core taxa present in >= 75%% of samples.\n", sum(core_results$is_core)))

# 9. Aesthetic Visualizations
message("\n--- Step 9: Generating Publication-Grade Visuals ---")
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

p_comp <- plot_composition(gut_clean, rank = "Phylum", top_n = 4, group_by = "treatment")
p_ord  <- plot_ordination(gut_clean, method = "rpca", color = "treatment", ellipse = TRUE, biplot = TRUE)
p_alph <- plot_alpha(gut_clean, metric = "hill_1", x = "treatment", test = TRUE)
p_volc <- plot_da_volcano(da_results)
p_core <- plot_core(core_results)

message(">>> SUCCESS: Full tidybiome pipeline completed with zero errors!")
