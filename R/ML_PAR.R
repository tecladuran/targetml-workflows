# Stratified Nested Cross-Validation Utilities

# Parallelized nested cross-validation workflows for PLS-DA and Random Forest models,
# including optional recursive feature elimination (RFE) and stratified fold creation.
# Designed for robust performance evaluation in high-dimensional GC-IMS analyses.


library(BiocParallel)

# MAKE STRATIFIED FOLDS -------------------------------------------------------

make_folds <- function(..., k) {
  set.seed(1)
  folds <- vector("list", k)
  
  # Capture all input variables as a data frame
  vars <- list(...)
  if (length(vars) == 0) stop("You must provide at least one variable for stratification.")
  
  # Combine all variables into one data frame
  df <- as.data.frame(vars)
  
  # Create stratification factor (interaction between all variables)
  strata <- interaction(df, drop = TRUE)
  
  # Create folds maintaining proportional representation of each stratum
  for (lvl in levels(strata)) {
    idx <- which(strata == lvl)
    idx <- sample(idx)
    split_idx <- cut(seq_along(idx), breaks = k, labels = FALSE)
    for (i in 1:k) {
      folds[[i]] <- c(folds[[i]], idx[split_idx == i])
    }
  }
  
  return(folds)
}

# PARALLELIZED NESTED CV FUNCTIONS --------------------------------------------


plsda_nested_cv_rfe <- function(X, y, positive_class, negative_class, 
                                ncomp_max = 10, outer_folds = 5, inner_folds = 4) {
  
  y <- factor(y, levels = c(negative_class, positive_class))
  y_num <- ifelse(y == positive_class, 1, 0)
  
  outer_folds_idx <- make_folds(y, k=outer_folds)
  
  # --- PARALLEL OUTER LOOP ---
  results <- bplapply(seq_along(outer_folds_idx), function(o) {
    
    test_id <- outer_folds_idx[[o]]
    train_id <- setdiff(seq_along(y), test_id)
    
    X_train <- X[train_id, , drop = FALSE]
    y_train <- y[train_id]
    y_train_num <- y_num[train_id]
    X_test  <- X[test_id, , drop = FALSE]
    y_test  <- y[test_id]
    y_test_num <- y_num[test_id]
    
    features <- colnames(X_train)
    best_acc <- 0
    best_features <- features
    
    internal_acc_history <- data.frame()
    external_acc_history <- data.frame()
    
    repeat {
      inner_idx <- make_folds(y_train, k=inner_folds)
      acc_mat <- matrix(NA, nrow = inner_folds, ncol = ncomp_max)
      
      # --- INNER CV --------
      for (f in seq_along(inner_idx)) {
        val_id <- inner_idx[[f]]
        tr_id  <- setdiff(seq_along(y_train), val_id)
        
        Xi_tr <- X_train[tr_id, features, drop = FALSE]
        yi_tr <- y_train_num[tr_id]
        Xi_val <- X_train[val_id, features, drop = FALSE]
        yi_val <- y_train_num[val_id]
        
        model <- plsr(y ~ ., data = data.frame(y = yi_tr, Xi_tr),
                      ncomp = ncomp_max, method = "oscorespls", scale = TRUE)
        
        preds <- predict(model, newdata = Xi_val, ncomp = 1:ncomp_max)
        preds_mat <- matrix(preds, nrow = length(yi_val), ncol = ncomp_max)
        
        acc_mat[f, ] <- colMeans((preds_mat > 0.5) == yi_val)
      }
      
      mean_acc <- colMeans(acc_mat, na.rm = TRUE)
      opt_ncomp <- which.max(mean_acc)
      
      internal_acc_history <- rbind(
        internal_acc_history,
        data.frame(Nfeatures = length(features), Accuracy = max(mean_acc))
      )
      
      # Train full model
      model_full <- plsr(y ~ ., data = data.frame(y = y_train_num, 
                                                  X_train[, features, drop = FALSE]),
                         ncomp = opt_ncomp, method = "oscorespls", scale = TRUE)
      
      test_pred <- predict(model_full, X_test[, features, drop = FALSE], ncomp = opt_ncomp)
      test_class <- as.integer(test_pred > 0.5)
      acc_ext <- mean(test_class == y_test_num)
      
      external_acc_history <- rbind(
        external_acc_history,
        data.frame(Nfeatures = length(features), Accuracy = acc_ext)
      )
      
      if (max(mean_acc) > best_acc) {
        best_acc <- max(mean_acc)
        best_features <- features
        best_ncomp <- opt_ncomp
      }
      
      if (length(features) <= 10) break
      
      vip <- plsVarSel::VIP(model_full, opt.comp = opt_ncomp)
      drop_n <- max(1, ceiling(length(features) * 0.05))
      drop_feat <- names(sort(vip))[1:drop_n]
      features <- setdiff(features, drop_feat)
    }
    
    # Final model
    final_model <- plsr(
      y ~ ., 
      data = data.frame(y = y_train_num, X_train[, best_features, drop = FALSE]),
      ncomp = best_ncomp, method = "oscorespls", scale = TRUE
    )
    
    test_pred <- predict(final_model, X_test[, best_features, drop = FALSE], ncomp = best_ncomp)
    
    predictions <- data.frame(
      SampleIndex = test_id,
      TrueLabel = y_test_num,
      PredProb = as.numeric(test_pred),
      PredClass = as.integer(test_pred > 0.5),
      ncomp = best_ncomp
    )
    
    list(
      ncomp = best_ncomp,
      best_features = best_features,
      internal_acc_history = internal_acc_history,
      external_acc_history = external_acc_history,
      predictions = predictions
    )
  })
  
  # --- Collate predictions globally ---
  preds_all <- do.call(rbind, lapply(results, function(x) x$predictions))
  overall_acc <- mean(preds_all$TrueLabel == preds_all$PredClass)
  
  return(list(
    results = results,
    predictions = preds_all,
    overall_accuracy = overall_acc
  ))
}



rf_nested_cv_rfe <- function(X, y, positive_class, negative_class,
                             outer_folds = 5, inner_folds = 4,
                             ntree = 500, mtry_grid = c(2, 4, 6),
                             drop_rate = 0.05, strat_vars = NULL) {
  
  library(randomForest)
  
  y <- factor(y, levels = c(negative_class, positive_class))
  y_num <- ifelse(y == positive_class, 1, 0)
  
  if (is.null(strat_vars)) {
    outer_folds_idx <- make_folds(y, k = outer_folds)
  } else {
    outer_folds_idx <- do.call(make_folds, c(list(y), strat_vars, k = outer_folds))
  }
  
  results <- bplapply(seq_along(outer_folds_idx), function(o) {
    
    test_id <- outer_folds_idx[[o]]
    train_id <- setdiff(seq_along(y), test_id)
    
    X_train <- X[train_id, , drop = FALSE]
    y_train <- y[train_id]
    X_test  <- X[test_id, , drop = FALSE]
    y_test  <- y[test_id]
    y_test_num <- y_num[test_id]
    
    features <- colnames(X_train)
    best_acc <- 0
    best_features <- features
    best_mtry <- mtry_grid[1]
    
    internal_acc_history <- data.frame()
    external_acc_history <- data.frame()
    
    repeat {
      inner_idx <- make_folds(y_train, k = inner_folds)
      acc_mat <- matrix(NA, nrow = inner_folds, ncol = length(mtry_grid))
      
      # --- Inner CV --------
      for (f in seq_along(inner_idx)) {
        val_id <- inner_idx[[f]]
        tr_id <- setdiff(seq_along(y_train), val_id)
        
        Xi_tr <- X_train[tr_id, features, drop = FALSE]
        yi_tr <- y_train[tr_id]
        Xi_val <- X_train[val_id, features, drop = FALSE]
        yi_val <- y_train[val_id]
        
        for (m in seq_along(mtry_grid)) {
          model <- randomForest(
            x = Xi_tr, y = yi_tr,
            ntree = ntree, mtry = mtry_grid[m]
          )
          acc_mat[f, m] <- mean(predict(model, Xi_val) == yi_val)
        }
      }
      
      mean_acc <- colMeans(acc_mat, na.rm = TRUE)
      opt_mtry <- mtry_grid[which.max(mean_acc)]
      
      internal_acc_history <- rbind(
        internal_acc_history,
        data.frame(Nfeatures = length(features), Accuracy = max(mean_acc))
      )
      
      # Train full model 
      model_full <- randomForest(
        x = X_train[, features, drop = FALSE],
        y = y_train,
        ntree = ntree,
        mtry = opt_mtry,
        importance = TRUE
      )
      
      preds_test <- predict(model_full, X_test[, features, drop = FALSE])
      acc_ext <- mean(preds_test == y_test)
      
      external_acc_history <- rbind(
        external_acc_history,
        data.frame(Nfeatures = length(features), Accuracy = acc_ext)
      )
      
      if (max(mean_acc) > best_acc) {
        best_acc <- max(mean_acc)
        best_features <- features
        best_mtry <- opt_mtry
      }
      
      if (length(features) <= 10) break
      
      # Feature elimination
      imp <- importance(model_full, type = 2)[, 1]
      imp <- sort(imp)
      
      drop_n <- max(1, ceiling(length(features) * drop_rate))
      drop_feat <- names(imp)[1:drop_n]
      
      features <- setdiff(features, drop_feat)
    }
    
    # Final model
    final_model <- randomForest(
      x = X_train[, best_features, drop = FALSE],
      y = y_train,
      ntree = ntree,
      mtry = best_mtry
    )
    
    test_prob <- predict(final_model, X_test[, best_features, drop = FALSE], type = "prob")[, 2]
    test_class <- as.integer(test_prob > 0.5)
    
    predictions <- data.frame(
      SampleIndex = test_id,
      TrueLabel = y_test_num,
      PredProb = test_prob,
      PredClass = test_class
    )
    
    list(
      best_mtry = best_mtry,
      best_features = best_features,
      internal_acc_history = internal_acc_history,
      external_acc_history = external_acc_history,
      predictions = predictions
    )
  })
  
  # Global predictions
  preds_all <- do.call(rbind, lapply(results, function(z) z$predictions))
  overall_acc <- mean(preds_all$TrueLabel == preds_all$PredClass)
  
  return(list(
    results = results,
    predictions = preds_all,
    overall_accuracy = overall_acc
  ))
}




rf_nested_cv <- function(X, y, positive_class, negative_class,
                         outer_folds = 7, inner_folds = 6,
                         ntree = 500, mtry_grid = c(2, 4, 6),
                         strat_vars = NULL) {
  
  library(randomForest)
  y <- factor(y, levels = c(negative_class, positive_class))
  
  if (is.null(strat_vars)) {
    outer_folds_idx <- make_folds(y, k = outer_folds)
  } else {
    outer_folds_idx <- do.call(make_folds, c(list(y), strat_vars, k = outer_folds))
  }
  
  preds_all <- data.frame()
  
  # --- PARALLEL OUTER LOOP ---
  results <- bplapply(seq_along(outer_folds_idx), function(o) {
    
    test_id <- outer_folds_idx[[o]]
    
    train_id <- setdiff(seq_along(y), test_id)
    X_train <- X[train_id, , drop = FALSE]
    y_train <- y[train_id]
    X_test  <- X[test_id, , drop = FALSE]
    y_test  <- y[test_id]
    
    # Inner folds
    if (is.null(strat_vars)) {
      inner_idx <- make_folds(y_train, k = inner_folds)
    } else {
      strat_train <- lapply(strat_vars, function(v) v[train_id])
      inner_idx <- do.call(make_folds, c(list(y_train), strat_train, k = inner_folds))
    }
    
    inner_acc <- numeric(length(mtry_grid))
    
    for (m in seq_along(mtry_grid)) {
      acc_inner <- numeric(length(inner_idx))
      
      for (f in seq_along(inner_idx)) {
        inner_test <- inner_idx[[f]]
        inner_train <- setdiff(seq_along(y_train), inner_test)
        
        rf_inner <- randomForest(
          x = X_train[inner_train, , drop = FALSE],
          y = y_train[inner_train],
          ntree = ntree,
          mtry = mtry_grid[m]
        )
        
        pred_inner <- predict(rf_inner, newdata = X_train[inner_test, , drop = FALSE])
        acc_inner[f] <- mean(pred_inner == y_train[inner_test])
      }
      inner_acc[m] <- mean(acc_inner)
    }
    
    best_mtry <- mtry_grid[which.max(inner_acc)]
    
    # Final model
    rf_final <- randomForest(
      x = X_train,
      y = y_train,
      ntree = ntree,
      mtry = best_mtry
    )
    
    prob_mat <- predict(rf_final, newdata = X_test, type = "prob")
    pred_prob <- prob_mat[, positive_class]
    pred_class <- ifelse(pred_prob > 0.5, 1, 0)
    true_class <- ifelse(y_test == positive_class, 1, 0)
    
    preds_fold <- data.frame(
      SampleIndex = test_id,
      TrueLabel = true_class,
      PredProb = pred_prob,
      PredClass = pred_class
    )
    
    list(
      best_mtry = best_mtry,
      inner_acc = inner_acc,
      preds = preds_fold
    )
  })
  
  # --- Unify predictions ---
  preds_all <- do.call(rbind, lapply(results, function(z) z$preds))
  
  overall_acc <- mean(preds_all$TrueLabel == preds_all$PredClass)
  
  list(
    results = results,
    predictions = preds_all,
    overall_accuracy = overall_acc
  )
}


# SHOW RESULTS ---------------------------------------------------------------

show_results <- function(res_obj, title = "Model Results") {
  # --- Extract predictions and accuracy ---
  preds <- res_obj$predictions
  overall_acc <- res_obj$overall_accuracy
  
  # --- Confusion matrix ---
  conf_matrix <- table(True = preds$TrueLabel, Predicted = preds$PredClass)
  
  cat("###", title, "\n")
  cat("Overall Leave-One-Out Accuracy:", round(overall_acc, 3), "\n\n")
  
  # --- Show confusion matrix ---
  print(
    knitr::kable(
      conf_matrix,
      caption = "Global Confusion Matrix (Aggregated Across All LOO Iterations)."
    )
  )
  
  # --- ROC Curve ---
  library(pROC)
  library(ggplot2)
  
  roc_global <- roc(
    preds$TrueLabel,
    preds$PredProb,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  auc_val <- auc(roc_global)
  ci_val  <- ci.auc(roc_global)
  
  ggroc(roc_global, legacy.axes = TRUE, colour = "navy", linewidth = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    labs(
      title = paste0(
        "ROC Curve — AUC = ", round(auc_val, 3),
        " [", round(ci_val[1], 3), "–", round(ci_val[3], 3), "]"
      ),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    theme_minimal(base_size = 10)
  
}


