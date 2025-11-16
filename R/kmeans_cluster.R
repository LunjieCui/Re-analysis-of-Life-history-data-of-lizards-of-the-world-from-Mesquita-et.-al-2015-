kmeans_cluster <- function(
    X, K, nstart = 20, seed = 200202) {

  X <- as.data.frame(X, stringsAsFactors = FALSE)
  
  # z-score
  Z <-scale(X)
  
  set.seed(seed) # repeatable
  km <- stats::kmeans(Z, centers = K, nstart = nstart, algorithm = "Hartigan-Wong")
  
 message(sprintf("KMeans done!!!：K=%d；samples per cluster = %s",
                               K, paste(tabulate(km$cluster, nbins = K), collapse = ", ")))
  
  list(
    cluster = km$cluster,
    centers = km$centers,
    Z = Z,
    model = km
  )
}
