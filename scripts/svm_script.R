# # SVM Model for snoRNA Classification
# # Support Vector Machine implementation using e1071 package
# 
# # Load required libraries
# suppressPackageStartupMessages({
#   library(e1071)
#   library(caret)
#   library(ggplot2)
#   library(dplyr)
#   library(pROC)
#   library(kernlab)
# })
# 
# # Function to train and evaluate SVM model
# train_svm_model <- function(snorna_machine_learning, config = list()) {
#   
#   cat("Starting SVM model training...\n")
#   
#   # Default parameters
#   default_config <- list(
#     kernel = "radial",  # radial, linear, polynomial, sigmoid
#     cost = 1,
#     gamma = "scale",
#     test_size = 0.2,
#     cv_folds = 5,
#     tune_hyperparameters = TRUE
#   )
#   
#   # Merge with provided config
#   config <- modifyList(default_config, config)
#   
#   # --------------------------------------
#   # Step 1: Data preprocessing
#   # --------------------------------------
#   cat("Preprocessing data...\n")
#   
#   # Remove any rows with missing values
#   snorna_machine_learning <- na.omit(snorna_machine_learning)
#   
#   # Normalize numerical features
#   numeric_cols <- sapply(snorna_machine_learning, is.numeric)
#   preproc <- preProcess(snorna_machine_learning[, numeric_cols], 
#                         method = c("center", "scale"))
#   data_norm <- snorna_machine_learning
#   data_norm[, numeric_cols] <- predict(preproc, snorna_machine_learning[, numeric_cols])
#   
#   # Ensure target variable is factor
#   data_norm$snoRNA <- factor(data_norm$snoRNA)
#   
#   cat(sprintf("Dataset dimensions: %d samples, %d features\n", 
#               nrow(data_norm), ncol(data_norm) - 1))
#   cat(sprintf("Class distribution: %s\n", 
#               paste(table(data_norm$snoRNA), collapse = ", ")))
#   
#   # --------------------------------------
#   # Step 2: Train-test split
#   # --------------------------------------
#   set.seed(123)
#   train_idx <- createDataPartition(data_norm$snoRNA, p = 1 - config$test_size, list = FALSE)
#   train_data <- data_norm[train_idx, ]
#   test_data <- data_norm[-train_idx, ]
#   
#   cat(sprintf("Training set: %d samples\n", nrow(train_data)))
#   cat(sprintf("Test set: %d samples\n", nrow(test_data)))
#   
#   # --------------------------------------
#   # Step 3: Hyperparameter tuning (optional)
#   # --------------------------------------
#   if (config$tune_hyperparameters) {
#     cat("Performing hyperparameter tuning...\n")
#     
#     # Define parameter grid
#     if (config$kernel == "radial") {
#       tune_grid <- expand.grid(
#         cost = c(0.1, 1, 10, 100),
#         gamma = c(0.001, 0.01, 0.1, 1)
#       )
#     } else if (config$kernel == "linear") {
#       tune_grid <- expand.grid(
#         cost = c(0.1, 1, 10, 100),
#         gamma = NA
#       )
#     } else {
#       tune_grid <- expand.grid(
#         cost = c(0.1, 1, 10, 100),
#         gamma = c(0.001, 0.01, 0.1, 1)
#       )
#     }
#     
#     # Perform tuning
#     svm_tune <- tune(svm, snoRNA ~ ., data = train_data,
#                      kernel = config$kernel,
#                      ranges = list(cost = tune_grid$cost, 
#                                    gamma = if(config$kernel != "linear") tune_grid$gamma else NULL))
#     
#     best_params <- svm_tune$best.parameters
#     cat(sprintf("Best parameters: Cost = %.3f", best_params$cost))
#     if (!is.null(best_params$gamma)) {
#       cat(sprintf(", Gamma = %.6f", best_params$gamma))
#     }
#     cat("\n")
#     
#     # Train final model with best parameters
#     if (config$kernel == "linear") {
#       model <- svm(snoRNA ~ ., data = train_data,
#                    kernel = config$kernel,
#                    cost = best_params$cost,
#                    probability = TRUE)
#     } else {
#       model <- svm(snoRNA ~ ., data = train_data,
#                    kernel = config$kernel,
#                    cost = best_params$cost,
#                    gamma = best_params$gamma,
#                    probability = TRUE)
#     }
#   } else {
#     cat("Training SVM with default parameters...\n")
#     model <- svm(snoRNA ~ ., data = train_data,
#                  kernel = config$kernel,
#                  cost = config$cost,
#                  gamma = config$gamma,
#                  probability = TRUE)
#   }
#   
#   # --------------------------------------
#   # Step 4: Model evaluation
#   # --------------------------------------
#   cat("Evaluating model performance...\n")
#   
#   # Predictions on test set
#   preds <- predict(model, test_data, probability = TRUE)
#   pred_probs <- attr(preds, "probabilities")
#   
#   # Confusion matrix
#   cm <- confusionMatrix(preds, test_data$snoRNA, positive = "TRUE")
#   print(cm)
#   
#   # ROC curve (for binary classification)
#   if (length(levels(test_data$snoRNA)) == 2) {
#     if ("TRUE" %in% colnames(pred_probs)) {
#       roc_obj <- roc(test_data$snoRNA, pred_probs[, "TRUE"])
#       cat(sprintf("AUC: %.3f\n", auc(roc_obj)))
#     }
#   }

# Load required libraries
library(e1071)  # For SVM
library(caret)
library(ggplot2)
library(dplyr)

# Check if running with Snakemake or in RStudio
# Need to clear the environment first -> TODO(FIX)
rm(snakemake, envir = .GlobalEnv)
if (!exists("snakemake")) {
  # need to change the static hardcoded setwd -> maybe i do not need one (script that creates folder and such)
  setwd("C:/Users/Marco/RProject/snoMatcher/")
  # Create mock snakemake object for testing in Rstudio
  snakemake <- list(
    input = list(
      snorna_machine_learning = "results/intermediate/snorna_machine_learning.RData"
    ),
    output = list(
      
      svm_model = "results/models/svm_model.RData"
    ),
    config = list(
      generate_plots = TRUE,
      export_tables = TRUE
    )
  )
  
  # Helper function for mock object
  get_input <- function(name) snakemake$input[[name]]
  get_output <- function(name) snakemake$output[[name]]
  get_config <- function(name) snakemake$config[[name]]
  # get_threads <- function() snakemake$threads
  ("Debug execution")
  
} else {
  # Helper functions for real snakemake object
  get_input <- function(name) snakemake@input[[name]]
  get_output <- function(name) snakemake@output[[name]]
  get_config <- function(name) snakemake@config[[name]]
  # get_threads <- function() snakemake@threads
  ("snakemake execution")
}

load(file = get_input("snorna_machine_learning"))

# Normalize numerical features (e.g., distances)
preproc <- preProcess(snorna_machine_learning, method = c("center", "scale"))
data_norm <- predict(preproc, snorna_machine_learning)

data_norm$snoRNA <- factor(data_norm$snoRNA)

# Split data into training and testing sets
set.seed(123)
train_idx <- createDataPartition(data_norm$snoRNA, p = 0.8, list = FALSE)
train_data <- data_norm[train_idx, ]
test_data <- data_norm[-train_idx, ]

# Train the SVM model
model <- svm(
  snoRNA ~ ., 
  data = train_data,
  kernel = "radial",        # RBF kernel (equivalent to Gaussian)
  cost = 1,                 # Regularization parameter
  gamma = 1/ncol(train_data-1),  # Kernel parameter (default: 1/number of features)
  probability = TRUE        # Enable probability estimates
)

#===
# Step 2: Make prediction on test data
#===
pred <- predict(model, test_data, probability = TRUE)
pred_probs <- attr(pred, "probabilities")

# --------------------------------------
# Step 3: Evaluate performance
# --------------------------------------
# For SVM, predictions are already factor levels, no need to extract from probabilities
preds <- pred
confusionMatrix(preds, test_data$snoRNA, positive = "TRUE")

# --------------------------------------
# Step 4: Hyperparameter tuning (optional but recommended for SVM)
# --------------------------------------
# SVM performance is highly dependent on hyperparameters
# You might want to tune cost and gamma parameters
tune_results <- tune(svm, 
                     snoRNA ~ ., 
                     data = train_data,
                     kernel = "radial",
                     ranges = list(
                       cost = c(0.1, 1, 10, 100),
                       gamma = c(0.001, 0.01, 0.1, 1)
                     ),
                     tunecontrol = tune.control(cross = 5)
)

# Get best model from tuning
best_model <- tune_results$best.model

# Make predictions with tuned model
pred_tuned <- predict(best_model, test_data)
confusionMatrix(pred_tuned, test_data$snoRNA, positive = "TRUE")

# --------------------------------------
# Step 5: Feature importance (SVM doesn't have direct feature importance)
# --------------------------------------
# SVM doesn't provide feature importance like Random Forest
# Alternative approaches:

# Option 1: Permutation importance (similar concept to RF)
permutation_importance <- function(model, test_data, n_repeats = 10) {
  baseline_accuracy <- mean(predict(model, test_data) == test_data$snoRNA)
  
  importance_scores <- c()
  feature_names <- names(test_data)[names(test_data) != "snoRNA"]
  
  for (feature in feature_names) {
    accuracies <- c()
    for (i in 1:n_repeats) {
      # Create permuted dataset
      test_permuted <- test_data
      test_permuted[[feature]] <- sample(test_permuted[[feature]])
      
      # Calculate accuracy with permuted feature
      permuted_accuracy <- mean(predict(model, test_permuted) == test_data$snoRNA)
      accuracies <- c(accuracies, permuted_accuracy)
    }
    # Importance = decrease in accuracy when feature is permuted
    importance_scores <- c(importance_scores, baseline_accuracy - mean(accuracies))
  }
  
  names(importance_scores) <- feature_names
  return(importance_scores)
}

# Calculate permutation importance
importance_scores <- permutation_importance(best_model, test_data)
importance_df <- data.frame(
  Feature = names(importance_scores),
  Importance = importance_scores
) %>% arrange(desc(Importance))

# Plot importance
ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Feature Importance (SVM - Permutation Method)",
       x = "Features",
       y = "Importance (Accuracy Drop)")

# Option 2: Weight vector analysis (only for linear SVM)
# If you use kernel = "linear", you can examine the weight vector:
# linear_model <- svm(snoRNA ~ ., data = train_data, kernel = "linear")
# weights <- t(linear_model$coefs) %*% linear_model$SV
# This gives you the linear coefficients for each feature

save(model, best_model, file = get_output("svm_model"))
