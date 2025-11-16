prep_data <- function(
    path,
    drop_colsnum = c(8, 9, 11, 12, 16, 21, 22, 24),
    class_type = c("Mode.of.reproduction", "Foraging.Mode"),
    coverage_threshold = 0.5
) {

  tb <- read.csv(path, stringsAsFactors = FALSE)
  tb <- as.data.frame(tb, stringsAsFactors = FALSE)
  
  #drop those useless cols
  keep_cols <- setdiff(seq_len(ncol(tb)), drop_colsnum)
  tb <- tb[, keep_cols, drop = FALSE]
  
  # turn "" in char into NA
  char_cols_all <- names(tb)[vapply(tb, is.character, logical(1))]
  for (cn in char_cols_all) {
    tb[[cn]][tb[[cn]] == ""] <- NA
  }
  
  #labels: 1-6, features：7-
  labelcols_num  <- 1:6
  predictor_cols <- names(tb)[7:(ncol(tb)-1)]
  
  #delete some columes with too many NA to increase the data coverage
  coverage_ratio <- function(cols) {
    if (length(cols) == 0) return(0)
    mean(complete.cases(tb[, cols, drop = FALSE]))   #if the row contain no NA, return true
  }
  
  S <- predictor_cols   #initial setting 
  cov_S <- coverage_ratio(S)
  
  if (length(S) > 0 && cov_S < coverage_threshold) {     
    improved <- TRUE
    while (improved && length(S) > 0 && cov_S < coverage_threshold) {   #if there is improvement (improve == True), and there are still columes to delete and it has not met the thretholds
      best_cov <- -Inf
      best_idx <- NA
      
      for (i in seq_along(S)) {
        candidate <- setdiff(S, S[i]) # try to delete S[i]
        cov_i <- coverage_ratio(candidate) #calculate the coverage after deleting
        if (cov_i > best_cov) { #if the coverage is the largest, undate it as the best solution (so called "Greedy Algorithm")
          best_cov <- cov_i
          best_idx <- i
        }
      }
      if (is.na(best_idx)) {
        improved <- FALSE
      } else {dropcol <- S[best_idx]
        S <- setdiff(S, dropcol)
        cov_S <- best_cov
      }
    }
  }
  
  keep_rows <- rep(TRUE, nrow(tb)) # set all as true firstly
  keep_rows <- keep_rows & complete.cases(tb[, S, drop = FALSE]) # delete the NA in feature matrix
  keep_rows <- keep_rows & complete.cases(tb[, labelcols_num, drop = FALSE]) #delete the NA in labels
  
  tb_kept <- tb[keep_rows, , drop = FALSE]
  

  labels_mat <- as.matrix(tb_kept[, labelcols_num, drop = FALSE]) #select the labels part
  X_raw <- tb_kept[, S, drop = FALSE] #select the feature part
  
  
  #change categorized features (only two categories) as {0,1} (all the selected features contain only two categories)
  binary_map <- list()
  for (var in class_type) {
    if (var %in% names(X_raw)) {
      unique_vals <- unique(na.omit(X_raw[[var]]))
      if (length(unique_vals) == 2) {
        levels_sorted    <- sort(unique_vals)
        X_raw[[var]]     <- ifelse(X_raw[[var]] == levels_sorted[1], 0, 1)
        binary_map[[var]] <- setNames(c(0, 1), levels_sorted)
      } else {
        X_raw[[var]] <- as.factor(X_raw[[var]])
      }
    }
  }
  
  X_mm <- as.matrix(X_raw)
  storage.mode(X_mm) <- "double"
  
  # return a df with labels and X
  list(
    labels = labels_mat,
    X = X_mm
  )
}
