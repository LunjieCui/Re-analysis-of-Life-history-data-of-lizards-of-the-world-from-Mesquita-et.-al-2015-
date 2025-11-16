#required package: ggplot2,ggrepel, scales, caret, randomForest, cowplot, sf, rnaturalearth, rnaturalearthdata, ranger, doParallel
rm(list = ls())
graphics.off()
cat("\014")
gc()

library(ggplot2)
library(ggrepel)
library(scales)
library(caret)
library(randomForest)
library(cowplot)
library(rnaturalearth)
library(rnaturalearthdata)
library(ranger)
library(doParallel)

source("R/prep_data.R")
source("R/kmeans_cluster.R")
source("R/matching_env.R")
source("R/rf111.R")
source("R/clustersmapping.R")
source("R/compute_pairwise_env_phylo.R")
source("R/rfpppplot.R")
source("R/get_entropy.R")
source("R/graphing_entropy.R")

##----data cleaning and preparation-----
data_1p <- prep_data(
  path = "data/lizard.csv",
  coverage_threshold = 0.5
)

labels <- data_1p$labels
X <- data_1p$X
y <- data_1p$y

lon <- as.numeric(trimws(labels[,5]))
lat <- as.numeric(trimws(labels[,6]))
species_genus_family <- labels[,1:3]


##----cluster mapping-----
# 1) kmeans cluster
km <- kmeans_cluster(X, K = 4, nstart = 50, seed = 200202)


p <- clustersmapping(lon, lat, km$cluster)
print(p)
ggsave("cluster_map.png", p, width = 8, height = 6, dpi = 600)



##----Formation of distance matrix-----
#step1: environmental feature matching
env <- matching_env(
  lon = lon,
  lat = lat,
  csv_path = "data/TerraClimate_Lizard_Table_1970_2020_withElevation.csv",
  out_cols = c(6, 11)
)

#step2: generation of the distance metrix x (+genetic distance) and y
species_genus_family <- as.data.frame(species_genus_family)
Y_ori <- as.data.frame(X)
matrixXY <- compute_pairwise_env_phylo(
  env = env,
  species_genus_family = species_genus_family,
  phylo_path = "data/phylo_distances.csv",
  trait_feature = Y_ori
)
matrixX <- matrixXY[3:9]
matrixY <- matrixXY[10]

#step3: training & estimation & feature importance
#warning!~~~!  it takes a large amount of tiiime
rf_reg <- rf111(
  matrixX, matrixY,
  cv = 3,
  ntree = 300,
)

# graphing 
p_density <- rfpppplot(rf_reg, what = "density")
print(p_density)
ggplot2::ggsave("rf_density_contour.png", p_density, width = 6, height = 6, dpi = 300)

p_imp <- rfpppplot(rf_reg, what = "importance")
print(p_imp)
ggplot2::ggsave("rf_importance.png", p_imp, width = 6, height = 6, dpi = 300)



##----cluster at the genus scale-----
X_df <-as.data.frame(X)
entropy_matrix <- get_entropy(species_genus_family, env, X_df,
                               phylo_path = "data/phylo_distances.csv")

km_entropy <- kmeans_cluster(entropy_matrix[2:4], K=4, nstart = 50, seed = 43)

plots <- graphing_entropy(
  df = entropy_matrix,
  cluster = km_entropy$cluster,
  legend_title = "Cluster",
  # palette = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"), # 可选：自定义颜色
  label_size_pt = 10,
  base_size = 18
)

print(plots$IE_vs_IPG)
ggplot2::ggsave("ie_ipg.png", plots$IE_IPG, width = 12, height = 6, dpi = 600)
print(plots$IE_vs_ILH)
ggplot2::ggsave("ie_ilh.png", plots$IE_ILH, width = 12, height = 6, dpi = 600)
print(plots$IPG_vs_ILH)
ggplot2::ggsave("IPG_ILH.png", plots$IPG_ILH, width=12, height=6, dpi=600)
