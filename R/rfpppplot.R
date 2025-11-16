rfpppplot <- function(rfresult, what = c("density", "importance"),
                      save_path = NULL, width = 6, height = 6, dpi = 600) {
  what <- match.arg(what)
  
  if (what == "density") { #the first part -- point density plot
    y_te <- rfresult$y_test
    y_pred <- rfresult$y_pred
    R2 <- as.numeric(rfresult$test_metrics["R2"])
    RMSE <- as.numeric(rfresult$test_metrics["RMSE"])
    dfp <- data.frame(Observed = y_te, Predicted = y_pred)
    
    x_min <- min(dfp$Observed, dfp$Predicted, na.rm = TRUE)
    x_max <- max(dfp$Observed, dfp$Predicted, na.rm = TRUE)
    pad   <- 0.02 * (x_max - x_min + 1e-12)
    label_txt <- sprintf("R² = %.3f\nRMSE = %.3f", R2, RMSE)
    
    # point+density line
    p <- ggplot2::ggplot(dfp, ggplot2::aes(Observed, Predicted)) +  # modifying graphs
      ggplot2::geom_point(color = "grey40", alpha = 0.4, size = 0.8) +
      ggplot2::stat_density_2d(
        color = "black",
        size = 0.5,
        bins = 8
      ) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", color = "gray50") +
      ggplot2::annotate("text",
                        x = x_max - pad, y = x_min + pad,   
                        label = label_txt,
                        hjust = 1, vjust = 0,
                        size = 5, color = "black") +
      ggplot2::coord_equal(xlim = c(x_min, x_max), ylim = c(x_min, x_max)) +
      ggplot2::labs(x = "Observed difference", y = "Predicted difference") +
      ggplot2::theme_bw(base_size = 18) +
      ggplot2::theme(
        plot.title = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        axis.text  = ggplot2::element_text(color = "black"),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.line  = ggplot2::element_line(size = 0.8, color = "black"),
        axis.ticks = ggplot2::element_line(color = "black"),
        axis.ticks.length = ggplot2::unit(0.2, "cm"),
        panel.border = ggplot2::element_blank(),
        legend.position = "none",
        axis.line.x.top = ggplot2::element_blank(),
        axis.line.y.right = ggplot2::element_blank()
      )
    
    if (!is.null(save_path))
      ggplot2::ggsave(save_path, p, width = width, height = height, dpi = dpi)
    
    return(p)
  }
  

  if (what == "importance") { #the second part -- importance plot
    imp_mat <- randomForest::importance(rfresult$rf_imp, scale = TRUE)
    imp_df <- as.data.frame(imp_mat, stringsAsFactors = FALSE)
    imp_df$Variable <- rownames(imp_df)
    rownames(imp_df) <- NULL
    
    score_col <- "%IncMSE"
    
    imp_df[[score_col]] <- as.numeric(imp_df[[score_col]])
    imp_df <- imp_df[order(imp_df[[score_col]], decreasing = TRUE), ]
    
    p <- ggplot2::ggplot(
      imp_df,
      ggplot2::aes(x = stats::reorder(Variable, !!as.name(score_col)),
                   y = .data[[score_col]])
    ) +
      ggplot2::geom_point(size = 3, color = "black") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = score_col) +
      ggplot2::theme_bw(base_size = 18) +
      ggplot2::theme(
        plot.title = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        axis.text  = ggplot2::element_text(color = "black"),
        axis.title = ggplot2::element_text(face = "bold"),
        axis.line  = ggplot2::element_line(size = 0.8, color = "black"),
        axis.ticks = ggplot2::element_line(color = "black"),
        panel.border = ggplot2::element_blank()
      )
    
    if (!is.null(save_path))
      ggplot2::ggsave(save_path, p, width = width, height = height, dpi = dpi)
    
    return(p)
  }
}