graphing_entropy <- function(
    df,
    cluster = NULL,
    legend_title = "Cluster",
    palette = NULL,
    label_size_pt = 10,
    base_size = 18,
    point_alpha = 0.6,
    point_size = 1.5
) {
  
  need_cols <- c("Genus", "I_E", "I_PG", "I_LH")
  
  # remomve NA
  keep_idx <- complete.cases(df[, need_cols])
  df <- df[keep_idx, need_cols, drop = FALSE]
  
  # cluster
  has_cluster <- !is.null(cluster)
  if (has_cluster) {
    df$Cluster <- factor(cluster[keep_idx])
  }
  
  # size was weired. As gpt-5 suggested, this is due to the 1 pt in font size is not equal to 1 mm in point size
  pt2gg <- function(pt) pt / 2.845
  
  # Unified drawing function
  make_one <- function(xvar, yvar, xlab_expr, ylab_expr) {
    aes_base <- ggplot2::aes(.data[[xvar]], .data[[yvar]])
    if (has_cluster) aes_base$colour <- quote(Cluster)
    p <- ggplot2::ggplot(df, aes_base) +
      ggplot2::geom_point(alpha = point_alpha, size = point_size) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = Genus),
        size = pt2gg(label_size_pt),
        max.overlaps = Inf,
        box.padding = 0.2, point.padding = 0.1, segment.size = 0.2,
        show.legend = FALSE
      ) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::labs(x = xlab_expr, y = ylab_expr) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(color = "black"),
        axis.line = ggplot2::element_line(size = 0.8)
      )
    
    # color setting
    if (has_cluster) {
      if (is.null(palette)) {
        p <- p + ggplot2::scale_color_brewer(palette = "Set2", name = legend_title)
      } else {
        p <- p + ggplot2::scale_color_manual(values = palette, name = legend_title)
      }
    }
    
    p
  }
  
  # label
  label_PG <- bquote(I[PG])
  label_E  <- bquote(I[E])
  label_LH <- bquote(I[LH])
  
  list(
    IE_IPG  = make_one("I_PG", "I_E",  label_PG, label_E),
    IE_ILH  = make_one("I_LH", "I_E",  label_LH, label_E),
    IPG_ILH = make_one("I_LH", "I_PG", label_LH, label_PG)
  )
}
