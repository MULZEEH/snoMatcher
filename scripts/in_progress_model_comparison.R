library(Biostrings)
library(data.table)
library(stringr)
library(writexl)
library(readxl)
library(ggseqlogo)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(tidyverse)
library(patchwork)
library(yaml)


compare_testing <- function(){
  preds <- pred
  confusionMatrix(preds, test_data$snoRNA, positive = "TRUE")
  
}
# Load your data
load("your_snorna_data.RData")

# Prepare your data (assuming it's already normalized)
set.seed(123)
train_idx <- createDataPartition(your_data$snoRNA, p = 0.8, list = FALSE)
train_data <- your_data[train_idx, ]
test_data <- your_data[-train_idx, ]

# Run the comparison
results <- train_and_evaluate_models(train_data, test_data, "snoRNA", model, model)
plots <- create_model_comparison_plots(results, test_data, "snoRNA", model, model)

# Display plots
print(plots$accuracy_comparison)
print(plots$detailed_metrics)
