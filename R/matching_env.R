matching_env <- function(lon, lat, csv_path, out_cols = c(6, 11)) {

  env_featuretable <- read.csv(csv_path, check.names = FALSE)
  
  # find the lon/lat colunmes
  col_lat <- names(env_featuretable)[names(env_featuretable) == "latitude"]
  col_lon <- names(env_featuretable)[names(env_featuretable) == "longitude"]
  lat_envt <- env_featuretable[[col_lat]]
  ref_envt <- env_featuretable[[col_lon]]
  
  #colume needed to be extracted
  keep_cols <- names(env_featuretable)[out_cols[1]:out_cols[2]]
  
  # initial setting
  n <- length(lon)
  matched_vals <- matrix(NA, nrow = n, ncol = length(keep_cols))

  #find the nearest point. match the lon/lat between env table and lizards table
  for (i in seq_len(n)) {
    d <- (lat_envt - lat[i])^2 + (ref_envt - lon[i])^2
    j <- which.min(d)
    matched_vals[i, ] <- as.numeric(env_featuretable[j, keep_cols])
  }
  
  result <- as.data.frame(matched_vals)
  colnames(result) <- keep_cols
  return(result)
}
