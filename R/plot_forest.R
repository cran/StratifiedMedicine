### Forest Plot ###
plot_forest <- function(x, estimand_filter = NULL, ref_line = NULL) {

  # Combine parameter-estimates with rules ##
  parm <- x$param.dat

  # Join rules to parameter data
  if (!is.null(x$Rules)) {
    rules <- rbind(data.frame(Subgrps = "ovrl", Rules = "Overall",
                              stringsAsFactors = FALSE),
                   x$Rules)
    plot.dat <- dplyr::left_join(parm, rules, by = "Subgrps")
  } else {
    plot.dat <- parm
    plot.dat$Rules <- ifelse(plot.dat$Subgrps == "ovrl", "Overall",
                             paste("Node", plot.dat$Subgrps))
  }

  # Auto-detect the treatment-difference estimand (contains "diff", "-", or "vs")
  if (is.null(estimand_filter)) {
    diff_idx <- grepl("diff|\\-|vs|HR|RMST", plot.dat$estimand, ignore.case = TRUE)
    if (any(diff_idx)) {
      estimand_filter <- unique(plot.dat$estimand[diff_idx])[1]
    } else {
      # Fallback: use the last estimand (often the contrast)
      estimand_filter <- utils::tail(unique(plot.dat$estimand), 1)
    }
  }
  plot.dat <- plot.dat[plot.dat$estimand == estimand_filter, ]

  if (nrow(plot.dat) == 0) {
    stop("No rows match estimand_filter = '", estimand_filter, "'")
  }

  # Build display columns
  plot.dat$est_plot <- plot.dat$est0
  plot.dat$lcl_plot <- plot.dat$LCL0
  plot.dat$ucl_plot <- plot.dat$UCL0
  plot.dat$is_overall <- plot.dat$Subgrps == "ovrl"
  plot.dat$Rules_N <- paste0(plot.dat$Rules, "  (N=", plot.dat$N, ")")
  plot.dat$ci_label <- sprintf("%.2f [%.2f, %.2f]",
                               plot.dat$est_plot, plot.dat$lcl_plot,
                               plot.dat$ucl_plot)

  # Order: overall at bottom (first in rev), then subgroups by node id
  row_order <- order(!plot.dat$is_overall, as.character(plot.dat$Subgrps))
  ordered_labels <- plot.dat$Rules_N[row_order]
  plot.dat$Rules_N <- factor(plot.dat$Rules_N, levels = rev(ordered_labels))

  # Default reference line
  if (is.null(ref_line)) {
    ref_line <- ifelse(grepl("HR", estimand_filter, ignore.case = TRUE), 1, 0)
  }

  # Plot
  res <- ggplot2::ggplot(plot.dat,
                         ggplot2::aes(x = Rules_N, y = est_plot,
                                      ymin = lcl_plot, ymax = ucl_plot)) +
    ggplot2::geom_hline(yintercept = ref_line, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_pointrange(
      ggplot2::aes(colour = is_overall, size = is_overall),
      show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c("TRUE" = "black",
                                            "FALSE" = "steelblue")) +
    ggplot2::scale_size_manual(values = c("TRUE" = 0.8, "FALSE" = 0.5)) +
    ggplot2::geom_text(ggplot2::aes(label = ci_label),
                       hjust = -0.15, size = 3) +
    ggplot2::coord_flip() +
    ggplot2::labs(y = paste0(estimand_filter, "  (95% CI)"),
                  x = NULL,
                  title = "Forest Plot") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.y  = ggplot2::element_text(size = 10),
      axis.text.x  = ggplot2::element_text(face = "bold"),
      axis.title   = ggplot2::element_text(size = 12, face = "bold"),
      panel.grid.major.y = ggplot2::element_blank())

  return(res)
}
