clustersmapping <- function(lon,lat, cluster,
                            xlim = c(-130, -30),
                            ylim = c(-55, 50),
                            point_size = 2,
                            font_size =  18,
                            palette = "Set2") {

  df_lonlatclus <- data.frame(
    lon = lon,
    lat = lat,
    cluster = as.factor(cluster)
  )
  
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  bbox  <- sf::st_bbox(
    c(xmin = min(xlim), xmax = max(xlim),
      ymin = min(ylim), ymax = max(ylim)),
    crs = sf::st_crs(4326)
  )
  world_crop <- tryCatch(sf::st_crop(world, bbox), error = function(e) world)
  
  df_lonlatclus_in <- subset(df_lonlatclus,
                  lon >= min(xlim) & lon <= max(xlim) &
                    lat >= min(ylim) & lat <= max(ylim))
  
  p <- ggplot() +
    geom_sf(data = world_crop, fill = "grey95", color = "grey70", linewidth = 0.2) +
    geom_point(data = df_lonlatclus_in, aes(x = lon, y = lat, color = cluster),
               size = point_size, alpha = 0.9) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    scale_color_brewer(palette = palette, na.value = "black", name = "Cluster") +  
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme_void(base_size = font_size) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = font_size),
      legend.text = element_text(size = font_size),
      plot.title = element_blank()
    )
  
  return(p)
}
