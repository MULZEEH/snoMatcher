# Load required libraries
library(randomForest)
library(ranger)
library(caret)



# Normalize numerical features (e.g., distances)
preproc <- preProcess(snorna_machine_learning, method = c("center", "scale"))
data_norm <- predict(preproc, snorna_machine_learning)

data_norm$snoRNA <- factor(data_norm$snoRNA)


# Split data into training and testing sets
set.seed(123)
train_idx <- createDataPartition(data_norm$snoRNA, p = 0.8, list = FALSE)
train_data <- data_norm[train_idx, ]
test_data <- data_norm[-train_idx, ]

# Train the Random Forest model
model <- ranger(
  snoRNA ~ ., 
  data = train_data,
  importance = "permutation", # For feature importance
  num.trees = 500, # to check w
  mtry = sqrt(ncol(train_data) - 1),
  probability = T
)

# --------------------------------------
# Step 3: Evaluate performance
# --------------------------------------
preds <- colnames(pred$predictions)[max.col(pred$predictions, ties.method = "first")]
preds <- factor(preds, levels = levels(test_data$snoRNA))
confusionMatrix(preds, test_data$snoRNA, positive = "TRUE")


# --------------------------------------
# Step 4: Extract feature importance
# --------------------------------------
importance_scores <- importance(model)
importance_df <- data.frame(
  Feature = names(importance_scores),
  Importance = importance_scores
) %>% arrange(desc(Importance))


# Plot importance
ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Feature Importance (Random Forest)")
