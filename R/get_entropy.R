get_entropy <- function(
    species_genus_family, env, X,
    phylo_path = "data/phylo_distances.csv"
) {
  
  # as shown in the course "Specie Diversity & distribution", we calculate similarity, then calculate the weight
  dist_to_similarity <- function(D) {
    D <- as.matrix(D)
    if (nrow(D) == 1L) return(matrix(1, 1, 1, dimnames = dimnames(D)))
    mx <- max(D[upper.tri(D)], na.rm = TRUE)
    if (!is.finite(mx) || mx <= 0) {
      Z <- matrix(1, nrow(D), ncol(D)); dimnames(Z) <- dimnames(D); return(Z)
    }
    Z <- 1 - (D / mx)
    diag(Z) <- 1
    Z
  }
  
  simpson_weighted <- function(Z, p = NULL) {
    Z <- as.matrix(Z)
    n <- nrow(Z)
    if (is.null(p)) p <- rep(1 / n, n)
    p <- as.numeric(p)
    p <- p / sum(p)
    val <- as.numeric(1 - t(p) %*% Z %*% p)
    pmax(0, pmin(1, val))
  }
  
  env_scaled <- as.data.frame(scale(env))
  phylo_mat <- as.matrix(read.csv(phylo_path, row.names = 1, check.names = FALSE))
  if (!all(rownames(phylo_mat) == colnames(phylo_mat))) {
    sp_inter <- intersect(rownames(phylo_mat), colnames(phylo_mat)) # find the species that exsit both in rowname and columnname
    phylo_mat <- phylo_mat[sp_inter, sp_inter, drop = FALSE]
  }
  uptri <- phylo_mat[upper.tri(phylo_mat)] #find the "distance" in the upper triangle
  rng <- range(uptri, na.rm = TRUE)
  phylo_scaled <- (phylo_mat - rng[1]) / (rng[2] - rng[1])

  
  # genus name and species name
  sp_raw <- as.character(species_genus_family[[1]]) 
  genus  <- as.character(species_genus_family[[2]])
  sp_under <- gsub("\\s+", "_", sp_raw)
  
  genuslel <- sort(unique(genus)) #find all the genus name
  out <- data.frame(
    Genus = genuslel,
    I_E   = NA_real_, #will be filled with the environmental functional entropy
    I_PG  = NA_real_, #will be filled with the phylogenetic entropy
    I_LH  = NA_real_, #will be filled with the life history functional entropy
    stringsAsFactors = FALSE
  )
  
  for (g in seq_along(genuslel)) { #sort by each genus
    gname <- genuslel[g]
    idx   <- which(genus == gname)
    
    # I_E
    if (length(idx) >= 2) {
      D_env <- dist(as.matrix(env_scaled[idx, , drop = FALSE]), method = "euclidean")
      Z_env <- dist_to_similarity(as.matrix(D_env))
      p_env <- rep(1 / length(idx), length(idx))
      out$I_E[g] <- simpson_weighted(Z_env, p_env)
    }
    
    # I_PG
    spp_g <- unique(sp_under[idx])
    spp_g <- spp_g[spp_g %in% rownames(phylo_scaled)]
    if (length(spp_g) >= 2) {
      Z_pg <- 1 - phylo_scaled[spp_g, spp_g, drop = FALSE]
      diag(Z_pg) <- 1
      tab <- table(sp_under[idx])
      p_pg <- as.numeric(tab[spp_g])
      p_pg <- p_pg / sum(p_pg)
      out$I_PG[g] <- simpson_weighted(Z_pg, p_pg)
    }
    
    # I_LH
    if (length(idx) >= 2) {
      D_lh <- dist(as.matrix(X[idx, , drop = FALSE]), method = "euclidean")
      Z_lh <- dist_to_similarity(as.matrix(D_lh))
      p_lh <- rep(1 / length(idx), length(idx))
      out$I_LH[g] <- simpson_weighted(Z_lh, p_lh)
    }
  }
  
  fl <- na.omit(out)
  
  return(fl)
}