# SVM Model for snoRNA Classification
# Support Vector Machine implementation using e1071 package

# Load required libraries
suppressPackageStartupMessages({
  library(e1071)
  library(caret)
  library(ggplot2)
  library(dplyr)
  library(pROC)
  library(kernlab)
})

# Function to train and evaluate SVM model
train_svm_model <- function(snorna_machine_learning, config = list()) {
  
  cat("Starting SVM model training...\n")
  
  # Default parameters
  default_config <- list(
    kernel = "radial",  # radial, linear, polynomial, sigmoid
    cost = 1,
    gamma = "scale",
    test_size = 0.2,
    cv_folds = 5,
    tune_hyperparameters = TRUE
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
  
  # Ensure target variable is factor
  data_norm$snoRNA <- factor(data_norm$snoRNA)
  
  cat(sprintf("Dataset dimensions: %d samples, %d features\n", 
              nrow(data_norm), ncol(data_norm) - 1))
  cat(sprintf("Class distribution: %s\n", 
              paste(table(data_norm$snoRNA), collapse = ", ")))
  
  # --------------------------------------
  # Step 2: Train-test split
  # --------------------------------------
  set.seed(123)
  train_idx <- createDataPartition(data_norm$snoRNA, p = 1 - config$test_size, list = FALSE)
  train_data <- data_norm[train_idx, ]
  test_data <- data_norm[-train_idx, ]
  
  cat(sprintf("Training set: %d samples\n", nrow(train_data)))
  cat(sprintf("Test set: %d samples\n", nrow(test_data)))
  
  # --------------------------------------
  # Step 3: Hyperparameter tuning (optional)
  # --------------------------------------
  if (config$tune_hyperparameters) {
    cat("Performing hyperparameter tuning...\n")
    
    # Define parameter grid
    if (config$kernel == "radial") {
      tune_grid <- expand.grid(
        cost = c(0.1, 1, 10, 100),
        gamma = c(0.001, 0.01, 0.1, 1)
      )
    } else if (config$kernel == "linear") {
      tune_grid <- expand.grid(
        cost = c(0.1, 1, 10, 100),
        gamma = NA
      )
    } else {
      tune_grid <- expand.grid(
        cost = c(0.1, 1, 10, 100),
        gamma = c(0.001, 0.01, 0.1, 1)
      )
    }
    
    # Perform tuning
    svm_tune <- tune(svm, snoRNA ~ ., data = train_data,
                     kernel = config$kernel,
                     ranges = list(cost = tune_grid$cost, 
                                   gamma = if(config$kernel != "linear") tune_grid$gamma else NULL))
    
    best_params <- svm_tune$best.parameters
    cat(sprintf("Best parameters: Cost = %.3f", best_params$cost))
    if (!is.null(best_params$gamma)) {
      cat(sprintf(", Gamma = %.6f", best_params$gamma))
    }
    cat("\n")
    
    # Train final model with best parameters
    if (config$kernel == "linear") {
      model <- svm(snoRNA ~ ., data = train_data,
                   kernel = config$kernel,
                   cost = best_params$cost,
                   probability = TRUE)
    } else {
      model <- svm(snoRNA ~ ., data = train_data,
                   kernel = config$kernel,
                   cost = best_params$cost,
                   gamma = best_params$gamma,
                   probability = TRUE)
    }
  } else {
    cat("Training SVM with default parameters...\n")
    model <- svm(snoRNA ~ ., data = train_data,
                 kernel = config$kernel,
                 cost = config$cost,
                 gamma = config$gamma,
                 probability = TRUE)
  }
  
  # --------------------------------------
  # Step 4: Model evaluation
  # --------------------------------------
  cat("Evaluating model performance...\n")
  
  # Predictions on test set
  preds <- predict(model, test_data, probability = TRUE)
  pred_probs <- attr(preds, "probabilities")
  
  # Confusion matrix
  cm <- confusionMatrix(preds, test_data$snoRNA, positive = "TRUE")
  print(cm)
  
  # ROC curve (for binary classification)
  if (length(levels(test_data$snoRNA)) == 2) {
    if ("TRUE" %in% colnames(pred_probs)) {
      roc_obj <- roc(test_data$snoRNA, pred_probs[, "TRUE"])
      cat(sprintf("AUC: %.3f\n", auc(roc_obj)))
    }
  }
  
  # --------------------------------------
  # Step 5: Feature importance (approximation)
  # --------------------------------------
  cat("Computing feature importance...\n")
  
  # For SVM, we approximate feature importance using permutation importance
  feature_names <- names(train_data)[names(train_data) != "snoRNA"]
  importance_scores <- numeric(length(feature_names))
  names(importance_scores) <- feature_names
  
  # Baseline accuracy
  baseline_preds <- predict(model, test_data)
  baseline_acc <- mean(baseline_preds == test_data$snoRNA)
  
  # Permutation importance
  for (i in seq_along(feature_names)) {
    # Create permuted dataset
    test_permuted <- test_data
    test_permuted[[feature_names[i]]] <- sample(test_permuted[[feature_names[i]]])
    
    # Get predictions on permuted data
    permuted_preds <- predict(model, test_permuted)
    permuted_acc <- mean(permuted_preds == test_data$snoRNA)
    
    # Importance is the decrease in accuracy
    importance_scores[i] <- baseline_acc - permuted_acc
  }
  
  # Create importance dataframe
  importance_df <- data.frame(
    Feature = names(importance_scores),
    Importance = importance_scores
  ) %>% 
    arrange(desc(Importance)) %>%
    mutate(Importance = pmax(0, Importance))  # Ensure non-negative
  
  # --------------------------------------
  # Step 6: Visualizations
  # --------------------------------------
  
  # Feature importance plot
  importance_plot <- ggplot(importance_df[1:min(20, nrow(importance_df)), ], 
                            aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_bar(stat = "identity", fill = "#E74C3C") +
    coord_flip() +
    labs(title = "Feature Importance (SVM - Permutation Method)",
         x = "Features",
         y = "Importance Score") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # ROC plot (for binary classification)
  roc_plot <- NULL
  if (length(levels(test_data$snoRNA)) == 2 && exists("roc_obj")) {
    roc_plot <- ggroc(roc_obj) +
      labs(title = sprintf("ROC Curve (AUC = %.3f)", auc(roc_obj))) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  # --------------------------------------
  # Step 7: Return results
  # --------------------------------------
  results <- list(
    model = model,
    predictions = preds,
    probabilities = pred_probs,
    confusion_matrix = cm,
    importance_scores = importance_scores,
    importance_df = importance_df,
    preprocessing = preproc,
    plots = list(
      importance = importance_plot,
      roc = roc_plot
    ),
    performance_metrics = list(
      accuracy = cm$overall["Accuracy"],
      sensitivity = cm$byClass["Sensitivity"],
      specificity = cm$byClass["Specificity"],
      auc = if(exists("roc_obj")) auc(roc_obj) else NA
    ),
    config = config
  )
  
  cat("SVM model training completed!\n")
  return(results)
}

# Example usage function
run_svm_example <- function(snorna_machine_learning) {
  
  # Configuration for SVM
  svm_config <- list(
    kernel = "radial",
    tune_hyperparameters = TRUE,
    cv_folds = 5,
    test_size = 0.2
  )
  
  # Train model
  svm_results <- train_svm_model(snorna_machine_learning, svm_config)
  
  # Print results
  cat("\n=== SVM MODEL RESULTS ===\n")
  print(svm_results$performance_metrics)
  
  # Display plots
  if (!is.null(svm_results$plots$importance)) {
    print(svm_results$plots$importance)
  }
  
  if (!is.null(svm_results$plots$roc)) {
    print(svm_results$plots$roc)
  }
  
  return(svm_results)
}

# If running as standalone script
if (!interactive() && exists("snorna_machine_learning")) {
  svm_results <- run_svm_example(snorna_machine_learning)
}