# Load required libraries
library(randomForest)
library(ranger)
library(caret)



# Check if running with Snakemake or in RStudio
# Need to clear the environment first -> TODO(FIX)
rm(snakemake, envir = .GlobalEnv)
if (!exists("snakemake")) {
  # need to change the static hardcoded setwd
  setwd("C:/Users/Marco/RProject/snoMatcher/")
  # Create mock snakemake object for testing in Rstudio
  snakemake <- list(
    input = list(
      snorna_machine_learning = "results/intermediate/snorna_machine_learning.RData"
    ),
    output = list(
       model ="results/models/rf_model.RData",
      
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

# Train the Random Forest model
model <- ranger(
  snoRNA ~ ., 
  data = train_data,
  importance = "permutation", # For feature importance
  num.trees = 1000, # to check w -> add to config file
  mtry = sqrt(ncol(train_data) - 1),
  probability = T
)
#===
#Step 2: Make predition on test data
#===
pred <- predict(model, test_data)
# pred <- predict(model, train_data)

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

save(model, preproc, file = get_output("model"))