#' Modern Minimalist Publication Theme for tidybiome
#'
#' @param base_size Base font size (default: 11).
#' @param base_family Base font family.
#'
#' @return A `ggplot2` theme object.
#' @export
theme_tidybiome <- function(base_size = 11, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Typography
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.15), margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(color = "#555555", size = ggplot2::rel(0.95), margin = ggplot2::margin(b = 10)),
      plot.caption = ggplot2::element_text(color = "#777777", size = ggplot2::rel(0.8), margin = ggplot2::margin(t = 8)),
      axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.95)),
      axis.text = ggplot2::element_text(color = "#333333", size = ggplot2::rel(0.85)),

      # Grid and background
      panel.grid.major = ggplot2::element_line(color = "#ebebeb", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),

      # Facet strips
      strip.background = ggplot2::element_rect(fill = "#f4f4f4", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9), color = "#222222"),

      # Legend
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9)),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.85)),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank()
    )
}

#' Curated Color Palettes for tidybiome
#'
#' @param n Number of colors desired.
#' @param palette Palette name: `"tidybiome"` (Okabe-Ito colorblind-safe) or `"nature"` (Editorial journal palette).
#'
#' @return A character vector of hexadecimal color codes.
#' @export
tidybiome_pal <- function(n = NULL, palette = c("tidybiome", "nature")) {
  palette <- match.arg(palette)

  okabe_ito <- c(
    "#0072B2", "#D55E00", "#009E73", "#E69F00",
    "#CC79A7", "#56B4E9", "#F0E442", "#495057",
    "#2b2d42", "#8d99ae"
  )

  nature <- c(
    "#3B4992", "#EE0000", "#008B45", "#631879",
    "#008280", "#BB0021", "#5F559B", "#A20056",
    "#808180", "#1B1919"
  )

  cols <- if (palette == "tidybiome") okabe_ito else nature

  if (is.null(n)) {
    return(cols)
  }
  if (n <= length(cols)) {
    return(cols[seq_len(n)])
  }
  # Interpolate if more colors requested
  grDevices::colorRampPalette(cols)(n)
}

#' @rdname tidybiome_pal
#' @param ... Arguments passed to [ggplot2::discrete_scale()].
#' @export
scale_color_tidybiome <- function(palette = "tidybiome", ...) {
  ggplot2::discrete_scale(
    aesthetics = "colour",
    palette = function(n) tidybiome_pal(n, palette = palette),
    ...
  )
}

#' @rdname tidybiome_pal
#' @export
scale_fill_tidybiome <- function(palette = "tidybiome", ...) {
  ggplot2::discrete_scale(
    aesthetics = "fill",
    palette = function(n) tidybiome_pal(n, palette = palette),
    ...
  )
}
