library(ggplot2)

plot_histogram <- function(
  df,
  x_var,
  x_label,
  y_label = NULL,
  title = NULL,
  fill_label = NULL,
  x_limits = waiver(),
  x_breaks = waiver(),
  x_minor_breaks = waiver(),
  binwidth = 0.025,
  histogram_position = "dodge",
  use_polar = FALSE,
  outfile = NULL
) {
  histogram <- ggplot(df, aes(x = .data[[x_var]], color = condition, fill = condition)) +
    geom_histogram(position = histogram_position, alpha = 0.2, binwidth = binwidth) +
    facet_wrap(~fix_at) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      strip.text.x = element_text(size = 12),
      strip.text.y = element_text(size = 12),
      legend.position = "bottom",
      legend.direction = "horizontal",
      panel.spacing = unit(2, "lines")
    ) +
    scale_color_manual(values = c('red', 'blue')) +
    scale_fill_manual(values = c('red', 'blue'))

    # Polar/radial plot
    if (use_polar) {
      histogram <- histogram + coord_polar(start = 300, direction = -1)
    }

    # Add optional plot labels
    label_args <- list()
    if (!is.null(x_label)) label_args$x <- x_label
    if (!is.null(y_label)) label_args$y <- y_label
    if (!is.null(title)) label_args$title <- title
    if (!is.null(fill_label)) label_args$fill <- fill_label
    if (length(label_args) > 0) {
      histogram <- histogram + do.call(labs, label_args)
    }

    # Add x scale with/without limits and breaks
    scale_args <- list(name = x_label)

    # Add limits if specified
    if (!inherits(x_limits, "waiver")) {
      scale_args$limits <- x_limits
    }

    # Add breaks if specified
    if (!inherits(x_breaks, "waiver")) {
      scale_args$breaks <- x_breaks
    }

    # Add minor_breaks if specified
    if (!inherits(x_minor_breaks, "waiver")) {
      scale_args$minor_breaks <- x_minor_breaks
    }

    # Add the scale with do.call
    histogram <- histogram + do.call(scale_x_continuous, scale_args)

  if (!is.null(outfile)) {
    ggsave(
      filename = outfile,
      plot = histogram,
      width = 8,
      height = 6,
      dpi = 300
    )
  }

  invisible(histogram)
}
