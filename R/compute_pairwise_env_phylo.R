compute_pairwise_env_phylo <- function(env, species_genus_family, phylo_path, trait_feature) {
  phylo_mat <- as.data.frame(read.csv(phylo_path, row.names = 1, check.names = FALSE))
  phylo_species <- rownames(phylo_mat)
  
  # suitable for dataframe and other
  trait_species <- if (is.data.frame(species_genus_family)) species_genus_family[[1]] else species_genus_family[,1]
  
  cbs <- t(combn(nrow(env), 2))
  n_pairs <- nrow(cbs) # number of pair
  n_feats <- ncol(env) #number of features in env
  

  diff_env <- matrix(NA_real_, nrow = n_pairs, ncol = n_feats) #empty matrix to contain the result of distance
  phylo_dist <- rep(NA_real_, n_pairs)
  trait_feature_dist <- rep(NA_real_, n_pairs)
  colnames(diff_env) <- paste0("Δ", names(env)) #name the col
  
  
  # zscore of the trait features
  trait_feature <- as.data.frame(trait_feature)
  mu <- vapply(trait_feature, mean, numeric(1), na.rm = TRUE) #as suggested by GPT-5, using vapply()here is more stable, making sure the output is numeric
  sdv <- vapply(trait_feature, sd,numeric(1), na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1 #cannot "delete 0" or “delete NA”
  trait_feature_scaled <- (X - rep(mu, each=nrow(X))) / rep(sdv, each=nrow(X))
  
  # calculate the distance matrix
  for (k in seq_len(n_pairs)) {
    i <- cbs[k, 1]; j <- cbs[k, 2]
    
    # absolute differences of environmental factors between records
    diff_env[k, ] <- abs(as.numeric(env[i, ]) - as.numeric(env[j, ]))
    
    # Euclidean distance of life history trait between records 
    v <- trait_feature_scaled[i, ] - trait_feature_scaled[j, ]
    trait_feature_dist[k] <- sqrt(sum(v * v))
    
    # phylogenetic distance between species
    sp1 <- gsub(" ", "_", trait_species[i], perl = TRUE) #in phylogenetic table, the name of species is "Genus_"
    sp2 <- gsub(" ", "_", trait_species[j], perl = TRUE)
    if (sp1 %in% phylo_species && sp2 %in% phylo_species) {
      phylo_dist[k] <- phylo_mat[sp1, sp2]
    } else {
      phylo_dist[k] <- NA_real_
    }
  }
  
  # return the result with species names, env table distance (differences), phylogenetic distance, life history trait distances
  result <- data.frame(
    Species1  = trait_species[cbs[, 1]],
    Species2  = trait_species[cbs[, 2]],
    diff_env,
    phylogeny = phylo_dist,
    trait_distance = trait_feature_dist,
    check.names = FALSE
  )
  
  # remove the row contains na
  result <- na.omit(result)
  
  return(result)
}
