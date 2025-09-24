# XGBoost Model for snoRNA Classification
# Extreme Gradient Boosting implementation

# Load required libraries
suppressPackageStartupMessages({
  library(xgboost)
  library(caret)
  library(ggplot2)
  library(dplyr)
  library(pROC)
  library(Matrix)
})

# Function to train and evaluate XGBoost model
train_xgboost_model <- function(snorna_machine_learning, config = list()) {
  
  cat("Starting XGBoost model training...\n")
  
  # Default parameters
  default_config <- list(
    nrounds = 100,
    max_depth = 6,
    eta = 0.3,
    subsample = 1,
    colsample_bytree = 1,
    gamma = 0,
    min_child_weight = 1,
    reg_alpha = 0,
    reg_lambda = 1,
    test_size = 0.2,
    cv_folds = 5,
    early_stopping_rounds = 10,
    tune_hyperparameters = TRUE,
    objective = "binary:logistic",
    eval_metric = "auc"
  )
  
  # Merge with provided config
  config <- modifyList(default_config, config)
  
  # --------------------------------------
  # Step 1: Data preprocessing
  # --------------------------------------
  cat("Preprocessing data...\n")
  
  # Remove any rows with missing values
  snorna_machine_learning <- na.omit(snorna_machine_learning)
  
  # Normalize numerical features
  numeric_cols <- sapply(snorna_machine_learning, is.numeric)
  preproc <- preProcess(snorna_machine_learning[, numeric_cols], 
                        method = c("center", "scale"))
  data_norm <- snorna_machine_learning
  data_norm[, numeric_cols] <- predict(preproc, snorna_machine_learning[, numeric_cols])
  
  # Convert target to numeric (XGBoost requirement)
  data_norm$snoRNA <- factor(data_norm$snoRNA)
  target_levels <- levels(data_norm$snoRNA)
  data_norm$snoRNA_numeric <- as.numeric(data_norm$snoRNA) - 1
  
  # Handle multi-class vs binary classification
  n_classes <- length(target_levels)
  if (n_classes > 2) {
    config$objective <- "multi:softprob"
    config$eval_metric <- "mlogloss"
    config$num_class <- n_classes
  }
  
  cat(sprintf("Dataset dimensions: %d samples, %d features\n", 
              nrow(data_norm), ncol(data_norm) - 2))
  cat(sprintf("Class distribution: %s\n", 
              paste(table(data_norm$snoRNA), collapse = ", ")))
  
  # --------------------------------------
  # Step 2: Train-test split
  # --------------------------------------
  set.seed(123)
  train_idx <- createDataPartition(data_norm$snoRNA, p = 1 - config$test_size, list = FALSE)
  train_data <- data_norm[train_idx, ]
  test_data <- data_norm[-train_idx, ]
  
  # Prepare matrices for XGBoost
  feature_cols <- !names(train_data) %in% c("snoRNA", "snoRNA_numeric")
  train_matrix <- xgb.DMatrix(
    data = as.matrix(train_data[, feature_cols]),
    label = train_data$snoRNA_numeric
  )
  
  test_matrix <- xgb.DMatrix(
    data = as.matrix(test_data[, feature_cols]),
    label = test_data$snoRNA_numeric
  )
  
  cat(sprintf("Training set: %d samples\n", nrow(train_data)))
  cat(sprintf("Test set: %d samples\n", nrow(test_data)))
  
  # --------------------------------------
  # Step 3: Hyperparameter tuning (optional)
  # --------------------------------------
  if (config$tune_hyperparameters) {
    cat("Performing hyperparameter tuning with cross-validation...\n")
    
    # Define parameter grid
    param_grid <- list(
      max_depth = c(3, 6, 9),
      eta = c(0.1, 0.3, 0.5),
      subsample = c(0.8, 1.0),
      colsample_bytree = c(0.8, 1.0)
    )
    
    best_params <- list()
    best_score <- ifelse(config$objective == "binary:logistic", 0, Inf)
    
    # Manual grid search with CV
    for (max_depth in param_grid$max_depth) {
      for (eta in param_grid$eta) {
        for (subsample in param_grid$subsample) {
          for (colsample_bytree in param_grid$colsample_bytree) {
            
            params <- list(
              objective = config$objective,
              eval_metric = config$eval_metric,
              max_depth = max_depth,
              eta = eta,
              subsample = subsample,
              colsample_bytree = colsample_bytree,
              gamma = config$gamma,
              min_child_weight = config$min_child_weight,
              reg_alpha = config$reg_alpha,
              reg_lambda = config$reg_lambda
            )
            
            if (n_classes > 2) {
              params$num_class <- n_classes
            }
            
            # Cross-validation
            cv_result <- xgb.cv(
              params = params,
              data = train_matrix,
              nrounds = config$nrounds,
              nfold = config$cv_folds,
              early_stopping_rounds = config$early_stopping_rounds,
              verbose = FALSE,
              showsd = FALSE
            )
            
            # Get best score
            if (config$objective == "binary:logistic") {
              score <- max(cv_result$evaluation_log$test_auc_mean, na.rm = TRUE)
              if (score > best_score) {
                best_score <- score
                best_params <- params
                best_nrounds <- cv_result$best_iteration
              }
            } else {
              score <- min(cv_result$evaluation_log[[paste0("test_", config$eval_metric, "_mean")]], na.rm = TRUE)
              if (score < best_score) {
                best_score <- score
                best_params <- params
                best_nrounds <- cv_result$best_iteration
              }
            }
          }
        }
      }
    }
    
    cat(sprintf("Best CV score: %.4f\n", best_score))
    cat(sprintf("Best nrounds: %d\n", best_nrounds))
    
    # Train final model with best parameters
    model <- xgboost(
      data = train_matrix,
      params = best_params,
      nrounds = best_nrounds,
      verbose = 0
    )
    
  } else {
    cat("Training XGBoost with default parameters...\n")
    
    params <- list(
      objective = config$objective,
      eval_metric = config$eval_metric,
      max_depth = config$max_depth,
      eta = config$eta,
      subsample = config$subsample,
      colsample_bytree = config$colsample_bytree,
      gamma = config$gamma,
      min_child_weight = config$min_child_weight,
      reg_alpha = config$reg_alpha,
      reg_lambda = config$reg_lambda
    )
    
    if (n_classes > 2) {
      params$num_class <- n_classes
    }
    
    model <- xgboost(
      data = train_matrix,
      params = params,
      nrounds = config$nrounds,
      verbose = 0
    )
  }
  
  # --------------------------------------
  # Step 4: Model evaluation
  # --------------------------------------
  cat("Evaluating model performance...\n")
  
  # Predictions on test set
  pred_probs <- predict(model, test_matrix)
  
  if (n_classes == 2) {
    # Binary classification
    preds <- factor(ifelse(pred_probs > 0.5, target_levels[2], target_levels[1]), 
                    levels = target_levels)
    pred_prob_matrix <- cbind(1 - pred_probs, pred_probs)
    colnames(pred_prob_matrix) <- target_levels
  } else {
    # Multi-class classification
    pred_prob_matrix <- matrix(pred_probs, ncol = n_classes, byrow = TRUE)
    colnames(pred_prob_matrix) <- target_levels
    preds <- factor(target_levels[apply(pred_prob_matrix, 1, which.max)], 
                    levels = target_levels)
  }
  
  # Confusion matrix
  cm <- confusionMatrix(preds, test_data$snoRNA, positive = target_levels[2])
  print(cm)
  
  # ROC curve (for binary classification)
  roc_obj <- NULL
  if (n_classes == 2) {
    roc_obj <- roc(test_data$snoRNA, pred_probs)
    cat(sprintf("AUC: %.3f\n", auc(roc_obj)))
  }
  
  # --------------------------------------
  # Step 5: Feature importance
  # --------------------------------------
  cat("Computing feature importance...\n")
  
  # Get feature importance from XGBoost
  importance_matrix <- xgb.importance(
    feature_names = colnames(train_matrix),
    model = model
  )
  
  # Create importance dataframe
  importance_df <- data.frame(
    Feature = importance_matrix$Feature,
    Importance = importance_matrix$Gain
  ) %>% arrange(desc(Importance))
  
  # --------------------------------------
  # Step 6: Visualizations
  # --------------------------------------
  
  # Feature importance plot
  importance_plot <- ggplot(importance_df[1:min(20, nrow(importance_df)), ], 
                            aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "#27AE60") +
    coord_flip() +
    labs(title = "Feature Importance (XGBoost - Gain)",
         x = "Features",
         y = "Importance Score") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # ROC plot (for binary classification)
  roc_plot <- NULL
  if (!is.null(roc_obj)) {
    roc_plot <- ggroc(roc_obj) +
      labs(title = sprintf("ROC Curve (AUC = %.3f)", auc(roc_obj))) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  # Learning curves (if available)
  learning_plot <- NULL
  
  # --------------------------------------
  # Step 7: Return results
  # --------------------------------------
  results <- list(
    model = model,
    predictions = preds,
    probabilities = pred_prob_matrix,
    confusion_matrix = cm,
    importance_scores = setNames(importance_df$Importance, importance_df$Feature),
    importance_df = importance_df,
    preprocessing = preproc,
    plots = list(
      importance = importance_plot,
      roc = roc_plot,
      learning = learning_plot
    ),
    performance_metrics = list(
      accuracy = cm$overall["Accuracy"],
      sensitivity = if(n_classes == 2) cm$byClass["Sensitivity"] else NA,
      specificity = if(n_classes == 2) cm$byClass["Specificity"] else NA,
      auc = if(!is.null(roc_obj)) auc(roc_obj) else NA
    ),
    config = config,
    target_levels = target_levels
  )
  
  cat("XGBoost model training completed!\n")
  return(results)
}

# Example usage function
run_xgboost_example <- function(snorna_machine_learning) {
  
  # Configuration for XGBoost
  xgboost_config <- list(
    nrounds = 200,
    max_depth = 6,
    eta = 0.1,
    tune_hyperparameters = TRUE,
    cv_folds = 5,
    test_size = 0.2,
    early_stopping_rounds = 20
  )
  
  # Train model
  xgboost_results <- train_xgboost_model(snorna_machine_learning, xgboost_config)
  
  # Print results
  cat("\n=== XGBOOST MODEL RESULTS ===\n")
  print(xgboost_results$performance_metrics)
  
  # Display plots
  if (!is.null(xgboost_results$plots$importance)) {
    print(xgboost_results$plots$importance)
  }
  
  if (!is.null(xgboost_results$plots$roc)) {
    print(xgboost_results$plots$roc)
  }
  
  return(xgboost_results)
}

# If running as standalone script
if (!interactive() && exists("snorna_machine_learning")) {
  xgboost_results <- run_xgboost_example(snorna_machine_learning)
}