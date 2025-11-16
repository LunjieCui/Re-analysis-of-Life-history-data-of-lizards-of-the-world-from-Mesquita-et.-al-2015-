rf111 <- function(
    X, y, test_size = 0.2,
    seed = 114514,
    cv = 3,
    ntree = 300,
    n_cores = max(1, parallel::detectCores()-1)
) {
  
  # as suggested by gpt-5, parallel computing shorten the time and avoid stuck
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  X <- as.data.frame(X)
  y <- y[,1]
  y <- as.numeric(y)
  set.seed(seed)
  n <- length(y)
  idx_tr <- caret::createDataPartition(y, p = 1 - test_size, list = FALSE) #Stratified sampling
  
  X_tr <- X[idx_tr, , drop = FALSE]; y_tr <- y[idx_tr]
  X_te <- X[-idx_tr, , drop = FALSE]; y_te <- y[-idx_tr]
  
  # CV
  ctrl <- caret::trainControl(method = "cv", number = cv) #setting of the k-fold
  
  fitrf <- caret::train(
    x = X_tr, y = y_tr,
    method = "ranger",
    trControl = ctrl,
    metric = "RMSE",
    tuneLength = 3,
    num.trees = ntree,
    importance = "impurity"
  )
  
  # print the best hyperparameters
  message("CV best hyperparameters："); print(fitrf$bestTune)
  message("CV accuracy：")
  print(fitrf$results[order(fitrf$results$mtry), c("mtry","RMSE","Rsquared")])
  
  # prediction and estimation using testing set
  y_pred <- predict(fitrf, newdata = X_te)
  metr <- caret::postResample(pred = y_pred, obs = y_te)
  RMSE <- unname(metr["RMSE"]); R2 <- unname(metr["Rsquared"])
  
  # importance
  best_mtry <- fitrf$bestTune$mtry
  rf_imp <- randomForest::randomForest(
    x = X_tr, y = y_tr,
    ntree = ntree, mtry = best_mtry,
    importance = TRUE
  )
  
  imp_mat <- randomForest::importance(rf_imp, scale = TRUE)
  imp_df <- data.frame(
    Variable = rownames(imp_mat),
    `%IncMSE` = imp_mat[, "%IncMSE"],
    IncNodePurity = imp_mat[, "IncNodePurity"],
    row.names = NULL, check.names = FALSE
  )
  
  list(
    model = fitrf,
    train_index = idx_tr,
    y_test = y_te,
    y_pred = y_pred,
    test_metrics = c(RMSE = RMSE, R2 = R2),
    rf_imp = rf_imp,
    imp_table = imp_df
  )
}
