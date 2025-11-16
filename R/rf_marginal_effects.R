# 文件：R/rf_ame_beeswarm.R
rf_marginal_effects <- function(
    model, X,
    target_class = NULL,
    group_onehot = TRUE,
    top_n = NULL,
    delta_frac = 0.05,
    baseline_levels = NULL,
    seed = 123,
    verbose = TRUE
) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("请先安装 ggplot2：install.packages('ggplot2')")
  
  # ---- 类别水平 ----
  cls_levels <- tryCatch(levels(model), error = function(e) NULL)
  if (is.null(cls_levels)) cls_levels <- tryCatch(model$finalModel$classes, error = function(e) NULL)
  if (is.null(cls_levels)) stop("无法从模型获取类别水平；请确认是分类模型且支持 type='prob'。")
  
  if (is.null(target_class)) target_classes <- cls_levels else {
    if (!all(target_class %in% cls_levels)) stop("target_class 不在模型类别中")
    target_classes <- target_class
  }
  
  # ---- 统一 X 类型 ----
  if (!is.data.frame(X)) X <- as.data.frame(X, stringsAsFactors = FALSE)
  if (anyNA(X)) stop("X 中存在 NA，请先清洗/插补。")
  is_num <- vapply(X, is.numeric,   logical(1))
  is_fac <- vapply(X, is.factor,    logical(1))
  is_chr <- vapply(X, is.character, logical(1))
  if (any(is_chr)) { for (cn in names(X)[is_chr]) X[[cn]] <- as.factor(X[[cn]]) ; is_fac <- vapply(X, is.factor, logical(1)) }
  is_binary_num <- is_num & vapply(X, function(v) all(na.omit(unique(v)) %in% c(0,1)), logical(1))
  is_binary_fac <- is_fac & vapply(X, function(v) length(levels(v)) == 2, logical(1))
  is_cont       <- is_num & !is_binary_num
  is_multifac   <- is_fac & vapply(X, function(v) length(levels(v)) > 2, logical(1))
  
  # 连续变量步长
  deltas <- rep(NA_real_, ncol(X))
  deltas[is_cont] <- vapply(X[is_cont], function(v){
    sdv <- stats::sd(v); ifelse(is.finite(sdv) && sdv > 0, sdv * delta_frac, NA_real_)
  }, numeric(1))
  
  # 概率预测器
  predict_prob_class <- function(newdata, class_name) {
    pr <- stats::predict(model, newdata = newdata, type = "prob")
    if (is.null(dim(pr))) stop("模型不支持 type='prob' 概率输出。")
    as.numeric(pr[, class_name])
  }
  
  set.seed(seed)
  pts <- list()
  
  # ---- 逐类别构造点云（Effect + 原始取值 Value）----
  for (cls in target_classes) {
    if (verbose) message(sprintf("计算 '%s' 的样本级 AME...", cls))
    
    # 连续
    for (v in names(X)[is_cont]) {
      d <- deltas[match(v, names(X))]; if (!is.finite(d) || d <= 0) next
      Xp <- X; Xn <- X
      vmin <- min(X[[v]], na.rm = TRUE); vmax <- max(X[[v]], na.rm = TRUE)
      Xp[[v]] <- pmin(vmax, X[[v]] + d); Xn[[v]] <- pmax(vmin, X[[v]] - d)
      p_plus  <- predict_prob_class(Xp, cls)
      p_minus <- predict_prob_class(Xn, cls)
      eff <- (p_plus - p_minus) / (2 * d)
      pts[[length(pts)+1]] <- data.frame(Class = cls, Variable = v, Effect = eff, Value = X[[v]], stringsAsFactors = FALSE)
    }
    
    # 二值数值
    for (v in names(X)[is_binary_num]) {
      X0 <- X; X1 <- X
      X0[[v]] <- 0; X1[[v]] <- 1
      diff <- predict_prob_class(X1, cls) - predict_prob_class(X0, cls)
      pts[[length(pts)+1]] <- data.frame(Class = cls, Variable = v, Effect = diff, Value = X[[v]], stringsAsFactors = FALSE)
    }
    
    # 二水平因子
    for (v in names(X)[is_binary_fac]) {
      lv <- levels(X[[v]])
      X0 <- X; X1 <- X
      X0[[v]] <- factor(lv[1], levels = lv); X1[[v]] <- factor(lv[2], levels = lv)
      diff <- predict_prob_class(X1, cls) - predict_prob_class(X0, cls)
      val  <- as.integer(X[[v]] == lv[2])
      pts[[length(pts)+1]] <- data.frame(Class = cls, Variable = v, Effect = diff, Value = val, stringsAsFactors = FALSE)
    }
    
    # 多水平因子
    for (v in names(X)[is_multifac]) {
      lv <- levels(X[[v]])
      base <- if (!is.null(baseline_levels) && !is.null(baseline_levels[[v]])) baseline_levels[[v]] else lv[1]
      if (!(base %in% lv)) base <- lv[1]
      Xb <- X; Xb[[v]] <- factor(base, levels = lv)
      pb <- predict_prob_class(Xb, cls)
      for (alt in lv) {
        Xa <- X; Xa[[v]] <- factor(alt, levels = lv)
        pa <- predict_prob_class(Xa, cls)
        eff <- pa - pb
        val <- as.integer(X[[v]] == alt)
        pts[[length(pts)+1]] <- data.frame(
          Class = cls,
          Variable = paste0(v, "=", alt),
          Effect = eff,
          Value  = val,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  pts <- do.call(rbind, pts)
  
  # ---- 分组名（是否聚合 one-hot 到因子层面）----
  factor_name_from <- function(vec_names, X_colnames) {
    split_positions <- function(s) which(strsplit(s, "")[[1]] %in% c(".", "_", ":"))
    prefix_candidates <- lapply(vec_names, function(s){
      pos <- split_positions(s); c(if (length(pos)) substring(s,1,pos) else character(0), s)
    })
    all_candidates <- unique(unlist(prefix_candidates))
    count_by_prefix <- setNames(vapply(all_candidates, function(pf) sum(startsWith(X_colnames, pf)), integer(1)), all_candidates)
    vapply(seq_along(vec_names), function(i){
      cands <- prefix_candidates[[i]]; cands <- cands[count_by_prefix[cands] >= 2]
      if (length(cands) == 0) vec_names[i] else cands[which.max(nchar(cands))]
    }, character(1))
  }
  
  key_col <- "Variable"
  if (group_onehot) {
    pts$Factor <- factor_name_from(pts$Variable, colnames(X))
    key_col <- "Factor"
  }
  pts$Key <- if (group_onehot) pts$Factor else pts$Variable
  
  # ---- 选择每类 Top N（按 |Effect| 平均强度）----
  if (!is.null(top_n)) {
    ord <- stats::aggregate(abs(Effect) ~ Key + Class, data = pts, FUN = mean)
    keep <- do.call(rbind, lapply(split(ord, ord$Class), function(d){
      d[order(d$`abs(Effect)`, decreasing = TRUE), ][seq_len(min(nrow(d), top_n)), , drop = FALSE]
    }))
    pts <- merge(pts, keep[, c("Key", "Class")], by = c("Key", "Class"))
  }
  
  # ---- 颜色标准化（每个 Key×Class 内 min-max 到 0~1）----
  pts <- do.call(rbind, lapply(split(pts, list(pts$Key, pts$Class), drop = TRUE), function(d){
    rng <- range(d$Value, na.rm = TRUE)
    d$ValueScaled <- if (is.finite(rng[1]) && diff(rng) > 0) (d$Value - rng[1]) / diff(rng) else 0.5
    d
  }))
  
  # ---- 计算 rank：|Effect| 越大 rank 越小；并做每类内排序 ----
  ord2 <- stats::aggregate(abs(Effect) ~ Key + Class, data = pts, FUN = mean)
  ord2 <- do.call(rbind, lapply(split(ord2, ord2$Class), function(d){
    d$rank <- rank(-d$`abs(Effect)`, ties.method = "first"); d
  }))
  pts <- merge(pts, ord2[, c("Key", "Class", "rank")], by = c("Key", "Class"), all.x = TRUE)
  
  # ---- 固定子图顺序 + 子图标题加序号 (a)/(b)/(c) ----
  pts$Class <- factor(pts$Class, levels = target_classes)  # 控制facet顺序
  strip_labs <- setNames(
    paste0("(", letters[seq_along(target_classes)], ") ", target_classes),
    target_classes
  )
  
  # ---- 作图（按 |AME| 从高到低：上→下）----
  library(ggplot2)
  
  # 保持 facet 顺序与编号
  pts$Class <- factor(pts$Class, levels = target_classes)
  strip_labs <- setNames(
    paste0("(", letters[seq_along(target_classes)], ") ", target_classes),
    target_classes
  )
  
  p <- ggplot(pts, aes(x = Effect, y = reorder(Key, -rank), color = ValueScaled)) +
    geom_point(alpha = 0.55, size = 0.9, position = position_jitter(height = 0.2)) +
    scale_color_gradient(low = "#2166AC", high = "#B2182B", name = "Feature value") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    facet_wrap(~Class, scales = "free_y", labeller = ggplot2::labeller(Class = strip_labs)) +
    labs(
      title = paste0("SHAP-like AME Beeswarm (", if (group_onehot) "by Factor" else "by Feature", ")"),
      x = "Signed AME (Δ probability)", y = NULL
    ) +
    theme_minimal(base_size = 12)
  
  
  list(
    points = pts,
    plot   = p
  )
}
