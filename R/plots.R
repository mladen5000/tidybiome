#' Plot Taxonomic Composition with Intelligent Top-N Grouping
#'
#' @description
#' Generates publication-ready stacked bar charts of taxonomic composition.
#' Automatically aggregates less abundant taxa into an aesthetically styled `"Other"`
#' category at the base of the plot, eliminating cluttered legends.
#'
#' @param tb A `tidy_microbiome` object.
#' @param rank Character string specifying taxonomic rank (e.g. `"Genus"`, `"Phylum"`).
#'   If `NULL`, uses ASV/feature IDs.
#' @param top_n Number of top most abundant taxa to display individually (default: 8).
#' @param group_by Optional character string specifying a sample metadata column to facet by.
#' @param palette Palette name: `"tidybiome"` or `"nature"`.
#'
#' @return A `ggplot2::ggplot` object.
#' @export
#'
#' @examples
#' counts <- matrix(c(100, 20, 5, 2, 5, 10, 80, 70), nrow = 2, byrow = TRUE,
#'                  dimnames = list(c("Tax1", "Tax2"), c("S1", "S2", "S3", "S4")))
#' sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"))
#' tax_table <- data.frame(taxon_id = c("Tax1", "Tax2"), Genus = c("Bacteroides", "Prevotella"))
#' tb <- tidy_microbiome(counts, sample_data, tax_table)
#' p <- plot_composition(tb, rank = "Genus", top_n = 2)
plot_composition <- function(tb,
                             rank = "Genus",
                             top_n = 8,
                             group_by = NULL,
                             palette = "tidybiome") {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  if (!is.null(rank)) {
    tb <- aggregate_taxa(tb, rank = rank)
  }

  df_long <- tidy_abundance(tb, assay = "counts", long = TRUE)
  tax_col <- if (!is.null(rank) && rank %in% colnames(df_long)) rank else "taxon_id"

  # Determine top N taxa by total abundance
  total_by_tax <- tapply(df_long$abundance, df_long[[tax_col]], sum, na.rm = TRUE)
  sorted_taxa <- names(sort(total_by_tax, decreasing = TRUE))
  top_taxa <- head(sorted_taxa, top_n)

  # Pool remaining taxa into "Other"
  df_long$Taxon_Display <- ifelse(df_long[[tax_col]] %in% top_taxa,
                                  as.character(df_long[[tax_col]]),
                                  "Other")

  # Order factors so Other is at the bottom
  all_levels <- c(top_taxa, "Other")
  df_long$Taxon_Display <- factor(df_long$Taxon_Display, levels = rev(all_levels))

  # Build color mapping
  top_colors <- tidybiome_pal(length(top_taxa), palette = palette)
  names(top_colors) <- top_taxa
  color_map <- c(top_colors, Other = "#d3d3d3")

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data$sample_id, y = .data$abundance, fill = .data$Taxon_Display)) +
    ggplot2::geom_col(position = "fill", width = 0.85, color = "white", linewidth = 0.2) +
    ggplot2::scale_y_continuous(labels = function(y) paste0(round(y * 100), "%"), expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = color_map, name = if (!is.null(rank)) rank else "Taxon") +
    ggplot2::labs(
      title = paste("Relative Abundance -", if (!is.null(rank)) rank else "Taxa"),
      x = "Sample",
      y = "Relative Abundance"
    ) +
    theme_tidybiome() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(group_by) && group_by %in% colnames(df_long)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[group_by]]), scales = "free_x")
  }

  p
}

#' Plot Microbiome Ordination and Compositional Biplots
#'
#' @description
#' Generates elegant 2D ordination plots (RPCA, PCoA, PCA) with automatic variance
#' explained axis annotations, confidence ellipses, and optional taxon loading biplot vectors.
#'
#' @param tb A `tidy_microbiome` object. Run [calc_ordination()] beforehand or specify `method`.
#' @param method Character string specifying ordination method (`"rpca"`, `"pcoa"`, or `"pca"`).
#' @param color Optional character string naming the metadata column to color points by.
#' @param shape Optional character string naming the metadata column for point shapes.
#' @param ellipse Logical. If `TRUE` and `color` is provided, draws 95% confidence ellipses around groups.
#' @param biplot Logical. If `TRUE` (and `method = "rpca"`), overlays top taxon loading vectors as arrows.
#' @param top_taxa Integer number of top driving taxa to annotate if `biplot = TRUE` (default: 5).
#' @param palette Palette name: `"tidybiome"` or `"nature"`.
#'
#' @return A `ggplot2::ggplot` object.
#' @export
#'
#' @examples
#' counts <- matrix(c(100, 20, 5, 2, 5, 10, 80, 70, 1, 1, 50, 40), nrow = 3, byrow = TRUE,
#'                  dimnames = list(c("Tax1", "Tax2", "Tax3"), c("S1", "S2", "S3", "S4")))
#' sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"))
#' tb <- tidy_microbiome(counts, sample_data)
#' tb <- calc_ordination(tb, method = "rpca")
#' p <- plot_ordination(tb, method = "rpca", color = "group")
plot_ordination <- function(tb,
                            method = c("rpca", "pcoa", "pca"),
                            color = NULL,
                            shape = NULL,
                            ellipse = TRUE,
                            biplot = FALSE,
                            top_taxa = 5,
                            palette = "tidybiome") {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  method <- match.arg(method)
  ord_list <- attr(tb, "metadata")$ordinations
  if (!method %in% names(ord_list)) {
    tb <- calc_ordination(tb, method = method)
    ord_list <- attr(tb, "metadata")$ordinations
  }

  ord <- ord_list[[method]]
  sample_df <- ord$samples
  var_exp <- ord$variance_explained

  coord_cols <- grep("^(PC|PCoA|MDS)[0-9]+$", colnames(sample_df), value = TRUE)
  if (length(coord_cols) < 2) {
    coord_cols <- c("PC1", "PC2")
  }
  x_col <- coord_cols[1]
  y_col <- coord_cols[2]

  x_lab <- paste0(x_col, " (", round(var_exp[1] * 100, 1), "%)")
  y_lab <- paste0(y_col, " (", round(var_exp[2] * 100, 1), "%)")

  aes_args <- list(x = rlang::sym(x_col), y = rlang::sym(y_col))
  if (!is.null(color) && color %in% colnames(sample_df)) {
    aes_args$color <- rlang::sym(color)
  }
  if (!is.null(shape) && shape %in% colnames(sample_df)) {
    aes_args$shape <- rlang::sym(shape)
  }

  p <- ggplot2::ggplot(sample_df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(size = 3.5, alpha = 0.85) +
    scale_color_tidybiome(palette = palette) +
    ggplot2::labs(
      title = paste("Microbiome Ordination -", toupper(method)),
      x = x_lab,
      y = y_lab
    ) +
    theme_tidybiome()

  # Ellipses
  if (ellipse && !is.null(color) && color %in% colnames(sample_df)) {
    grp_counts <- table(sample_df[[color]])
    if (all(grp_counts >= 3)) {
      p <- p + ggplot2::stat_ellipse(ggplot2::aes(fill = .data[[color]]),
                                     geom = "polygon", alpha = 0.15, linetype = "dashed",
                                     show.legend = FALSE) +
        scale_fill_tidybiome(palette = palette)
    }
  }

  # Biplot taxon vectors for RPCA
  if (biplot && !is.null(ord$taxa)) {
    taxa_df <- ord$taxa
    taxa_df$mag <- sqrt(taxa_df[[x_col]]^2 + taxa_df[[y_col]]^2)
    top_t <- head(taxa_df[order(-taxa_df$mag), ], top_taxa)

    # Scale arrows to sample range
    x_range <- range(sample_df[[x_col]], na.rm = TRUE)
    y_range <- range(sample_df[[y_col]], na.rm = TRUE)
    sf <- 0.7 * min(diff(x_range) / max(abs(taxa_df[[x_col]])),
                    diff(y_range) / max(abs(taxa_df[[y_col]])))

    top_t$x_end <- top_t[[x_col]] * sf
    top_t$y_end <- top_t[[y_col]] * sf

    label_col <- if ("Genus" %in% colnames(top_t)) "Genus" else "taxon_id"

    p <- p +
      ggplot2::geom_segment(data = top_t,
                            ggplot2::aes(x = 0, y = 0, xend = .data$x_end, yend = .data$y_end),
                            arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
                            color = "#2b2d42", linewidth = 0.7, inherit.aes = FALSE) +
      ggplot2::geom_text(data = top_t,
                         ggplot2::aes(x = .data$x_end * 1.1, y = .data$y_end * 1.1, label = .data[[label_col]]),
                         color = "#2b2d42", fontface = "bold", size = 3.2, inherit.aes = FALSE)
  }

  p
}

#' Plot Alpha Diversity with Statistical Comparisons
#'
#' @description
#' Generates publication-grade hybrid raincloud/boxplot visuals for alpha diversity metrics
#' (including Hill numbers \eqn{q = 0, 1, 2}), with overlaid sample points and automated
#' two-group statistical significance annotations.
#'
#' @param tb A `tidy_microbiome` object.
#' @param metric Character string naming the alpha diversity metric (e.g. `"hill_1"`, `"hill_0"`, `"shannon"`).
#' @param x Character string naming the metadata column for the horizontal grouping axis.
#' @param color Optional character string naming the metadata column to color by.
#' @param test Logical. If `TRUE` and `x` has 2 groups, runs a Wilcoxon rank-sum test and displays the p-value.
#' @param palette Palette name: `"tidybiome"` or `"nature"`.
#'
#' @return A `ggplot2::ggplot` object.
#' @export
#'
#' @examples
#' counts <- matrix(c(25, 25, 25, 25, 10, 5, 2, 80), nrow = 4, ncol = 2,
#'                  dimnames = list(paste0("T", 1:4), c("S1", "S2")))
#' sample_data <- data.frame(sample_id = c("S1", "S2"), treatment = c("Control", "Treated"))
#' tb <- tidy_microbiome(counts, sample_data)
#' tb <- calc_alpha_diversity(tb, metrics = "hill")
#' p <- plot_alpha(tb, metric = "hill_1", x = "treatment")
plot_alpha <- function(tb,
                       metric = "hill_1",
                       x = NULL,
                       color = NULL,
                       test = TRUE,
                       palette = "tidybiome") {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  sample_df <- tibble::as_tibble(tb)
  if (!metric %in% colnames(sample_df)) {
    # Auto-calculate
    calc_group <- if (grepl("^hill", metric)) "hill" else metric
    tb <- calc_alpha_diversity(tb, metrics = calc_group)
    sample_df <- tibble::as_tibble(tb)
  }

  if (is.null(x)) {
    x <- "All"
    sample_df$All <- "Samples"
  }

  if (is.null(color)) {
    color <- x
  }

  p <- ggplot2::ggplot(sample_df, ggplot2::aes(x = .data[[x]], y = .data[[metric]], color = .data[[color]])) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[color]]),
                          width = 0.35, alpha = 0.25, outlier.shape = NA, linewidth = 0.6) +
    ggplot2::geom_jitter(width = 0.15, size = 3, alpha = 0.8) +
    scale_color_tidybiome(palette = palette) +
    scale_fill_tidybiome(palette = palette) +
    ggplot2::labs(
      title = paste("Alpha Diversity -", metric_display_name(metric)),
      x = x,
      y = metric_display_name(metric)
    ) +
    theme_tidybiome() +
    ggplot2::theme(legend.position = if (x == color) "none" else "right")

  # Add significance test subtitle if 2 groups
  if (test && length(unique(sample_df[[x]])) == 2) {
    grps <- split(sample_df[[metric]], sample_df[[x]])
    wt <- tryCatch(stats::wilcox.test(grps[[1]], grps[[2]]), error = function(e) NULL)
    if (!is.null(wt)) {
      p_val <- wt$p.value
      p_text <- if (p_val < 0.001) "p < 0.001" else paste0("p = ", round(p_val, 3))
      p <- p + ggplot2::labs(subtitle = paste0("Wilcoxon rank-sum test: ", p_text))
    }
  }

  p
}

#' Plot Consensus Differential Abundance Volcano Plot
#'
#' @description
#' Plots effect size (\eqn{\log_2\text{FC}}) versus statistical significance
#' (\eqn{-\log_{10}(p_{\text{adj}})}), sizing points by their multi-engine consensus
#' agreement score and labeling top biomarker taxa.
#'
#' @param da_res A data frame or tibble produced by [calc_differential_abundance()].
#' @param fdr_cutoff False discovery rate significance threshold (default: 0.05).
#' @param fc_cutoff Absolute log2 fold-change cutoff (default: 1.0).
#' @param label_top Number of top significant taxa to label with names (default: 8).
#'
#' @return A `ggplot2::ggplot` object.
#' @export
#'
#' @examples
#' counts <- matrix(c(100, 120, 110, 5, 8, 6, 50, 45, 55, 48, 52, 50), nrow = 2, byrow = TRUE,
#'                  dimnames = list(c("DiffTaxon", "NullTaxon"), paste0("S", 1:6)))
#' sample_data <- data.frame(sample_id = paste0("S", 1:6), group = c("A", "A", "A", "B", "B", "B"))
#' tb <- tidy_microbiome(counts, sample_data)
#' da_df <- calc_differential_abundance(tb, group = "group")
#' p <- plot_da_volcano(da_df)
plot_da_volcano <- function(da_res,
                            fdr_cutoff = 0.05,
                            fc_cutoff = 1.0,
                            label_top = 8) {
  df <- as.data.frame(da_res)
  padj_col <- if ("padj_consensus" %in% colnames(df)) "padj_consensus" else "padj_linda"
  log2fc_col <- "log2fc"

  df$neg_log_p <- -log10(pmax(1e-15, df[[padj_col]]))

  # Assign status
  df$Status <- "Not Significant"
  sig_up   <- df[[padj_col]] <= fdr_cutoff & df[[log2fc_col]] >= fc_cutoff
  sig_down <- df[[padj_col]] <= fdr_cutoff & df[[log2fc_col]] <= -fc_cutoff
  df$Status[sig_up]   <- "Enriched"
  df$Status[sig_down] <- "Depleted"

  df$Status <- factor(df$Status, levels = c("Enriched", "Depleted", "Not Significant"))

  # Label column
  label_var <- if ("Genus" %in% colnames(df)) "Genus" else "taxon_id"
  df$Label <- ifelse(is.na(df[[label_var]]), df$taxon_id, as.character(df[[label_var]]))

  # Pick top biomarkers
  top_sig <- head(df[df$Status != "Not Significant", ][order(df[df$Status != "Not Significant", ][[padj_col]]), ], label_top)

  size_aes <- if ("agreement_score" %in% colnames(df)) "agreement_score" else NULL

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[log2fc_col]], y = .data$neg_log_p)) +
    ggplot2::geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "#888888", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "#888888", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(color = .data$Status, size = if (!is.null(size_aes)) .data[[size_aes]] else NULL),
                        alpha = 0.75) +
    ggplot2::scale_color_manual(values = c(Enriched = "#D55E00", Depleted = "#0072B2", `Not Significant` = "#999999")) +
    ggplot2::labs(
      title = "Consensus Differential Abundance Volcano",
      subtitle = sprintf("Cutoffs: FDR \u2264 %.2f, |log2FC| \u2265 %.1f", fdr_cutoff, fc_cutoff),
      x = expression(bold(log[2]~Fold~Change)),
      y = expression(bold(-log[10]~Adjusted~italic(p)-value)),
      size = "Agreement"
    ) +
    theme_tidybiome()

  # Add text labels for top taxa
  if (nrow(top_sig) > 0) {
    p <- p + ggplot2::geom_text(data = top_sig,
                                ggplot2::aes(label = .data$Label),
                                vjust = -0.6, size = 3.2, fontface = "bold", color = "#222222")
  }

  p
}

#' Plot Core Microbiome Landscape
#'
#' @description
#' Visualizes the core microbiome distribution across prevalence and abundance.
#'
#' @param core_res A data frame or tibble produced by [calc_core_microbiome()].
#' @return A `ggplot2::ggplot` object.
#' @export
plot_core <- function(core_res) {
  df <- as.data.frame(core_res)
  label_var <- if ("Genus" %in% colnames(df)) "Genus" else "taxon_id"
  df$Label <- ifelse(is.na(df[[label_var]]), df$taxon_id, as.character(df[[label_var]]))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$mean_abundance, y = .data$prevalence, color = .data$is_core)) +
    ggplot2::geom_point(size = 3.5, alpha = 0.8) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#009E73", `FALSE` = "#999999"), name = "Core Member") +
    ggplot2::scale_x_log10(labels = function(x) paste0(round(x * 100, 2), "%")) +
    ggplot2::scale_y_continuous(labels = function(y) paste0(round(y * 100), "%")) +
    ggplot2::labs(
      title = "Core Microbiome Landscape",
      x = "Mean Relative Abundance (Log Scale)",
      y = "Prevalence (% of samples)"
    ) +
    theme_tidybiome()
}

#' Plot Distance-Based Redundancy Analysis (db-RDA) Biplot
#'
#' @param dbrda_res An object produced by [run_dbrda()].
#' @param color Optional metadata variable to color sample points.
#' @param shape Optional metadata variable for point shapes.
#' @param palette Palette name: `"tidybiome"` or `"nature"`.
#'
#' @return A `ggplot2::ggplot` object.
#' @export
plot_dbrda <- function(dbrda_res, color = NULL, shape = NULL, palette = "tidybiome") {
  if (!inherits(dbrda_res, "tidybiome_dbrda")) {
    stop("`dbrda_res` must be an object produced by `run_dbrda()`.", call. = FALSE)
  }
  samples_df <- dbrda_res$samples
  var_exp <- dbrda_res$variance_explained

  x_lab <- if (length(var_exp) >= 1) paste0("dbRDA1 (", round(var_exp[1] * 100, 1), "%)") else "dbRDA1"
  y_lab <- if (length(var_exp) >= 2) paste0("dbRDA2 (", round(var_exp[2] * 100, 1), "%)") else "dbRDA2"

  aes_args <- list(x = rlang::sym("dbRDA1"), y = rlang::sym("dbRDA2"))
  if (!is.null(color) && color %in% colnames(samples_df)) {
    aes_args$color <- rlang::sym(color)
  }
  if (!is.null(shape) && shape %in% colnames(samples_df)) {
    aes_args$shape <- rlang::sym(shape)
  }

  p <- ggplot2::ggplot(samples_df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(size = 3.5, alpha = 0.85) +
    scale_color_tidybiome(palette = palette) +
    ggplot2::labs(
      title = "Constrained Ordination - db-RDA",
      x = x_lab,
      y = y_lab
    ) +
    theme_tidybiome()

  # Biplot vectors
  biplot_df <- dbrda_res$biplot
  if (nrow(biplot_df) > 0 && all(c("dbRDA1", "dbRDA2") %in% colnames(biplot_df))) {
    x_range <- range(samples_df$dbRDA1, na.rm = TRUE)
    y_range <- range(samples_df$dbRDA2, na.rm = TRUE)
    sf <- 0.7 * min(diff(x_range) / max(1e-6, max(abs(biplot_df$dbRDA1))),
                    diff(y_range) / max(1e-6, max(abs(biplot_df$dbRDA2))))

    biplot_df$x_end <- biplot_df$dbRDA1 * sf
    biplot_df$y_end <- biplot_df$dbRDA2 * sf

    p <- p +
      ggplot2::geom_segment(data = biplot_df,
                            ggplot2::aes(x = 0, y = 0, xend = .data$x_end, yend = .data$y_end),
                            arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
                            color = "#2b2d42", linewidth = 0.8, inherit.aes = FALSE) +
      ggplot2::geom_text(data = biplot_df,
                         ggplot2::aes(x = .data$x_end * 1.15, y = .data$y_end * 1.15, label = .data$term),
                         color = "#2b2d42", fontface = "bold", size = 3.5, inherit.aes = FALSE)
  }
  p
}

#' Plot Microbial Co-Occurrence Network
#'
#' @param net An object produced by [calc_network()].
#' @param color_by Taxonomy column to color nodes by (e.g. `"Phylum"`). Defaults to `"Phylum"`.
#' @param min_degree Minimum degree for a node to be displayed (default: 1).
#' @param palette Palette name: `"tidybiome"` or `"nature"`.
#'
#' @return A `ggplot2::ggplot` object.
#' @export
plot_network <- function(net, color_by = "Phylum", min_degree = 1, palette = "tidybiome") {
  if (!inherits(net, "tidybiome_network")) {
    stop("`net` must be an object produced by `calc_network()`.", call. = FALSE)
  }
  nodes <- net$nodes
  edges <- net$edges

  nodes_sub <- nodes[nodes$degree >= min_degree, ]
  if (nrow(nodes_sub) == 0) {
    nodes_sub <- nodes
  }

  n_nodes <- nrow(nodes_sub)
  theta <- seq(0, 2 * pi, length.out = n_nodes + 1)[seq_len(n_nodes)]
  nodes_sub$x <- cos(theta)
  nodes_sub$y <- sin(theta)

  coord_map_x <- stats::setNames(nodes_sub$x, nodes_sub$taxon_id)
  coord_map_y <- stats::setNames(nodes_sub$y, nodes_sub$taxon_id)

  edges_sub <- edges[edges$from %in% nodes_sub$taxon_id & edges$to %in% nodes_sub$taxon_id, ]
  edges_sub$x_start <- coord_map_x[edges_sub$from]
  edges_sub$y_start <- coord_map_y[edges_sub$from]
  edges_sub$x_end   <- coord_map_x[edges_sub$to]
  edges_sub$y_end   <- coord_map_y[edges_sub$to]

  color_col <- if (!is.null(color_by) && color_by %in% colnames(nodes_sub)) color_by else "taxon_id"

  p <- ggplot2::ggplot()

  if (nrow(edges_sub) > 0) {
    p <- p + ggplot2::geom_segment(
      data = edges_sub,
      ggplot2::aes(x = .data$x_start, y = .data$y_start, xend = .data$x_end, yend = .data$y_end,
                   color = .data$direction, linewidth = .data$weight),
      alpha = 0.5
    ) +
    ggplot2::scale_color_manual(values = c(positive = "#0072B2", negative = "#D55E00"), name = "Association") +
    ggplot2::scale_linewidth_continuous(range = c(0.4, 1.5), guide = "none")
  }

  label_var <- if ("Genus" %in% colnames(nodes_sub)) "Genus" else "taxon_id"
  nodes_sub$Label <- ifelse(is.na(nodes_sub[[label_var]]), nodes_sub$taxon_id, as.character(nodes_sub[[label_var]]))

  p <- p +
    ggplot2::geom_point(
      data = nodes_sub,
      ggplot2::aes(x = .data$x, y = .data$y, size = .data$mean_abundance, fill = .data[[color_col]]),
      shape = 21, color = "white", stroke = 1
    ) +
    ggplot2::geom_text(
      data = nodes_sub,
      ggplot2::aes(x = .data$x * 1.15, y = .data$y * 1.15, label = .data$Label),
      size = 3.0, fontface = "bold"
    ) +
    ggplot2::scale_size_continuous(range = c(3, 8), name = "Abundance") +
    ggplot2::labs(title = "Microbial Co-Occurrence Network") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5))

  p
}

metric_display_name <- function(metric) {
  switch(
    metric,
    hill_0 = "Species Richness (Hill q=0)",
    hill_1 = "Exponential Shannon (Hill q=1)",
    hill_2 = "Inverse Simpson (Hill q=2)",
    shannon = "Shannon Diversity Index (H)",
    simpson = "Gini-Simpson Index (1-D)",
    inv_simpson = "Inverse Simpson Index (1/D)",
    observed = "Observed Species Richness",
    chao1 = "Chao1 Richness Estimator",
    faith_pd = "Faith's Phylogenetic Diversity",
    metric
  )
}

#' Plot Compositional Microbiome Heatmap
#'
#' @description
#' Generates an aesthetic, publication-ready abundance heatmap for top taxa across samples,
#' featuring hierarchical clustering of samples and taxa, compositional scaling,
#' and optional metadata grouping.
#'
#' @param tb A `tidy_microbiome` object.
#' @param rank Optional taxonomic rank to aggregate by before plotting (e.g. `"Genus"`).
#' @param top_n Number of top abundant taxa to display (default: 25).
#' @param assay Assay to extract and plot. Defaults to `"counts"`.
#' @param scale Scaling applied to values: `"log10"`, `"relabundance"`, `"rclr"`, or `"none"`.
#' @param cluster_samples Logical. If `TRUE` (default), clusters samples via hierarchical clustering.
#' @param cluster_taxa Logical. If `TRUE` (default), clusters taxa via hierarchical clustering.
#' @param annotation_col Optional sample metadata column to group or facet samples by.
#' @param palette Palette option: `"viridis"`, `"magma"`, or `"plasma"`.
#'
#' @return A `ggplot2::ggplot` object.
#' @export
#'
#' @examples
#' counts <- matrix(c(100, 20, 5, 2, 5, 10, 80, 70), nrow = 2, byrow = TRUE,
#'                  dimnames = list(c("Tax1", "Tax2"), c("S1", "S2", "S3", "S4")))
#' sample_data <- data.frame(sample_id = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"))
#' tax_table <- data.frame(taxon_id = c("Tax1", "Tax2"), Genus = c("Bacteroides", "Prevotella"))
#' tb <- tidy_microbiome(counts, sample_data, tax_table)
#' p <- plot_heatmap(tb, top_n = 2)
plot_heatmap <- function(tb,
                         rank = NULL,
                         top_n = 25,
                         assay = "counts",
                         scale = c("log10", "relabundance", "rclr", "none"),
                         cluster_samples = TRUE,
                         cluster_taxa = TRUE,
                         annotation_col = NULL,
                         palette = "viridis") {
  if (!inherits(tb, "tidy_microbiome")) {
    stop("`tb` must be a `tidy_microbiome` object.", call. = FALSE)
  }

  scale <- match.arg(scale)

  if (!is.null(rank)) {
    tb <- aggregate_taxa(tb, rank = rank)
  }

  mat <- assay(tb, assay)
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop("Assay matrix is empty.", call. = FALSE)
  }

  # Filter to top N taxa by total abundance
  top_n <- min(top_n, nrow(mat))
  row_sums <- rowSums(mat, na.rm = TRUE)
  top_idx <- order(row_sums, decreasing = TRUE)[seq_len(top_n)]
  mat_sub <- mat[top_idx, , drop = FALSE]

  # Apply scaling
  mat_scaled <- switch(
    scale,
    log10 = log10(mat_sub + 1),
    relabundance = {
      cs <- colSums(mat_sub)
      cs[cs == 0] <- 1
      sweep(mat_sub, 2, cs, "/")
    },
    rclr = calc_rclr_matrix(mat_sub),
    none = mat_sub
  )

  # Cluster samples
  if (cluster_samples && ncol(mat_scaled) > 2) {
    dist_s <- stats::dist(t(mat_scaled))
    hc_s <- stats::hclust(dist_s)
    sample_order <- colnames(mat_scaled)[hc_s$order]
  } else {
    sample_order <- colnames(mat_scaled)
  }

  # Cluster taxa
  if (cluster_taxa && nrow(mat_scaled) > 2) {
    dist_t <- stats::dist(mat_scaled)
    hc_t <- stats::hclust(dist_t)
    taxa_order <- rownames(mat_scaled)[hc_t$order]
  } else {
    taxa_order <- rownames(mat_scaled)
  }

  n_samp <- ncol(mat_scaled)
  n_taxa <- nrow(mat_scaled)

  df <- tibble::tibble(
    taxon_id = factor(rep(rownames(mat_scaled), times = n_samp), levels = taxa_order),
    sample_id = factor(rep(colnames(mat_scaled), each = n_taxa), levels = sample_order),
    value = as.vector(mat_scaled)
  )

  # Join metadata if annotation_col requested
  meta <- tibble::as_tibble(tb)
  if (!is.null(annotation_col) && annotation_col %in% colnames(meta)) {
    df <- dplyr::left_join(df, meta[, c("sample_id", annotation_col)], by = "sample_id")
  }

  legend_label <- switch(
    scale,
    log10 = "log10(Counts + 1)",
    relabundance = "Relative Abundance",
    rclr = "Robust CLR",
    none = "Abundance"
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$sample_id, y = .data$taxon_id, fill = .data$value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::scale_fill_viridis_c(option = palette, name = legend_label) +
    ggplot2::labs(
      title = paste("Abundance Heatmap", if (!is.null(rank)) paste0("(", rank, ")") else ""),
      x = "Sample",
      y = if (!is.null(rank)) rank else "Taxon"
    ) +
    theme_tidybiome() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(annotation_col) && annotation_col %in% colnames(df)) {
    p <- p + ggplot2::facet_grid(
      cols = ggplot2::vars(.data[[annotation_col]]),
      scales = "free_x",
      space = "free_x"
    )
  }

  p
}
