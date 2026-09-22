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
gut_clean <- calc_alpha_diversity(gut_clean, metrics = c("hill", "simplex_variation", "shannon"))
cat("Hill numbers and Simplex Compositional Variation (bioRxiv 2026):\n")
print(gut_clean |> select(sample_id, treatment, hill_0, hill_1, hill_2, simplex_variation) |> head(4))

# 5. Beta Diversity (Robust Aitchison & Optimal Transport Tree-Wasserstein)
message("\n--- Step 5: Beta Diversity (Tree-Wasserstein & Robust Aitchison) ---")
gut_clean <- calc_beta_diversity(gut_clean, metric = "raitchison")
gut_clean <- calc_beta_diversity(gut_clean, metric = "wasserstein")
d_w <- get_distance(gut_clean, "wasserstein")
cat(sprintf("Calculated Optimal Transport Tree-Wasserstein Distance (dim %dx%d).\n",
            attr(d_w, "Size"), attr(d_w, "Size")))

# 5b. Distance-based Test for Homogeneity (DTH, bioRxiv 2025)
message("\n--- Step 5b: Distance-based Test for Homogeneity (DTH, bioRxiv 2025) ---")
dth_res <- dth_test(gut_clean, group = "treatment", metric = "wasserstein", n_perm = 199)
print(dth_res)

# 6. Dimensionality Reduction (RPCA Biplot & Contrastive PCA)
message("\n--- Step 6: Dimensionality Reduction (RPCA Biplot & Contrastive PCA) ---")
gut_clean <- calc_ordination(gut_clean, method = "rpca")
gut_clean <- calc_ordination(gut_clean, method = "cpca", contrast_group = "treatment")
ord <- get_ordination(gut_clean, "rpca")
ord_cpca <- get_ordination(gut_clean, "cpca")
cat(sprintf("RPCA Variance Explained: PC1 = %.1f%%, PC2 = %.1f%%\n",
            ord$variance_explained[1] * 100, ord$variance_explained[2] * 100))
cat(sprintf("Contrastive PCA (cPCA): isolated %d contrastive axes of treatment variance.\n",
            ncol(ord_cpca$samples) - ncol(gut_clean)))

# 7. SOTA Multi-Engine Consensus Differential Abundance with CAFT
message("\n--- Step 7: Consensus Differential Abundance (CAFT, LinDA, CLR, Wilcoxon) ---")
da_results <- calc_differential_abundance(
  gut_clean,
  group = "treatment",
  methods = c("consensus", "caft", "linda", "clr_linear", "wilcoxon"),
  fdr_cutoff = 0.05
)

cat("Top Significant Biomarkers (Consensus with CAFT zero-cell engine):\n")
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

# 10. Ecosystem Connectors & mia Functions
message("\n--- Step 10: Ecosystem Connectors (vegan, phyloseq, mia) & Advanced Screening ---")
# 10a. Vegan export & PERMANOVA
veg <- to_vegan(gut_clean, assay = "counts")
print(veg)

perm_res <- run_permanova(gut_clean, ~ treatment + diet, permutations = 199)
cat("\nTidy PERMANOVA (vegan::adonis2):\n")
print(perm_res)

disp_res <- run_betadisper(gut_clean, group = "treatment", permutations = 199)
print(disp_res)

# 10b. Prevalence & Dominance
prev_tbl <- calc_prevalence(gut_clean)
cat("\nTop 3 Prevalent Taxa:\n")
print(head(prev_tbl[, c("taxon_id", "Phylum", "prevalence", "mean_abundance")], 3))

gut_dom <- calc_dominant(gut_clean, add_to_metadata = TRUE)
cat("\nIdentified Dominant Taxa per Sample (added to metadata):\n")
print(table(gut_dom$dominant_taxa))

# 10c. Sample Divergence from Control
gut_div <- calc_divergence(gut_clean, reference = list(treatment = "Control"), method = "bray")
cat("\nMean Bray-Curtis Divergence to Control Baseline:\n")
cat(sprintf("  * Control Samples:  %.4f\n", mean(gut_div$divergence[gut_div$treatment == "Control"])))
cat(sprintf("  * Treated Samples:  %.4f\n", mean(gut_div$divergence[gut_div$treatment == "Treated"])))

# 10d. Taxon x Metadata Cross-Association
assoc_res <- calc_cross_association(gut_clean, variables = c("age", "depth"))
cat("\nTop 4 Cross-Associations with Host Covariates:\n")
print(head(assoc_res[, c("taxon_id", "variable", "correlation", "padj")], 4))

# 10e. ape Tree Integration & UniFrac Distances
taxa_ids <- attr(gut_clean, "tax_table")$taxon_id
set.seed(123)
tree <- ape::rtree(length(taxa_ids), tip.label = taxa_ids)
gut_with_tree <- set_tree(gut_clean, tree)
cat(sprintf("\nSuccessfully attached ape::phylo tree with %d tips.\n", length(get_tree(gut_with_tree)$tip.label)))

# 11. New Parity & Rigor Features
message("\n--- Step 11: Parity & Rigor Suite (UniFrac, Formula DA, db-RDA, Networks) ---")
# 11a. UniFrac distances
gut_with_tree <- calc_beta_diversity(gut_with_tree, metric = "unifrac")
gut_with_tree <- calc_beta_diversity(gut_with_tree, metric = "wunifrac")
cat("Computed Unweighted UniFrac and Weighted Normalized UniFrac distances.\n")

# 11b. Formula DA with covariate adjustment
da_cov <- calc_differential_abundance(gut_clean, formula = ~ treatment + age, contrast = "treatment")
cat(sprintf("Ran formula-based DA with covariate adjustment (age). Significant taxa: %d\n", sum(da_cov$is_significant)))

# 11c. Distance-based Redundancy Analysis (db-RDA)
dbrda_res <- run_dbrda(gut_clean, ~ treatment + age, distance = "bray")
print(dbrda_res)
p_dbrda <- plot_dbrda(dbrda_res, color = "treatment")

# 11d. Microbial Co-Occurrence Network
net_res <- calc_network(gut_clean, method = "spearman", min_prevalence = 0.3, r_cutoff = 0.3, p_cutoff = 0.05)
print(net_res)
p_net <- plot_network(net_res)

# 11e. Data Integrity and Taxon Filtering
tb_filt_taxa <- filter_taxa(gut_with_tree, Phylum == "Bacteroidota")
cat(sprintf("Filtered taxa by Phylum: remaining taxa = %d, tree tips = %d.\n",
            nrow(attr(tb_filt_taxa, "tax_table")), length(get_tree(tb_filt_taxa)$tip.label)))
cat(sprintf("Container validation passed: %s\n", validate_tidy_microbiome(tb_filt_taxa)))

message("\n>>> SUCCESS: Full tidybiome pipeline completed with zero errors!")


