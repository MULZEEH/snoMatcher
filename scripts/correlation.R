setwd("C:/Users/Marco/RProject/snoMatcher/")

library(Biostrings)
library(GenomicRanges)
library(BSgenome)
library(seqinr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(corrplot)

load(file = "results/intermediate/processed_info_box.RData")

replace_minus_one_with_na <- function(df) {
  df[df == -1] <- NA
  df[df == 0] <- NA
  
  return(df)
}

pos_seq <- DNAStringSet(snodb_boxes[["DNA Sequence"]])
pos_seq

# looking for correlation between various length
length_seq <- as.data.frame(sapply(snodb_boxes["DNA Sequence"], nchar))

# library(dplyr)
# length_seq <- snodb_boxes %>%
#   mutate(length = nchar(`DNA Sequence`)) %>%
#   select(length)

length_c <- as.data.frame(sapply(snodb_boxes["c_seq"], nchar))
length_c_prime <- as.data.frame(sapply(snodb_boxes["c_prime_seq"], nchar))
length_d <- as.data.frame(sapply(snodb_boxes["d_seq"], nchar))
length_d_prime <- as.data.frame(sapply(snodb_boxes["d_prime_seq"], nchar))

pos_c <- as.data.frame(snodb_boxes["c_start"])
pos_c_prime <- as.data.frame(snodb_boxes["c_prime_start"])
pos_d <- as.data.frame(snodb_boxes["d_start"])
pos_d_prime <- as.data.frame(snodb_boxes["d_prime_start"])

dist_cpd <- as.data.frame(snodb_boxes["dist_c_prime_d"])
dist_cd <- as.data.frame(snodb_boxes["dist_c_d"])
dist_dpc <- as.data.frame(snodb_boxes["dist_c_d_prime"])
dist_cpdp <- as.data.frame(snodb_boxes["dist_d_prime_c_prime"])

guide1_length <- as.data.frame(sapply(snodb_boxes["guide1_seq"], nchar))
guide2_length <- as.data.frame(sapply(snodb_boxes["guide2_seq"], nchar))

guide1_start <- as.data.frame(snodb_boxes["guide1_start"])
guide2_start <- as.data.frame(snodb_boxes["guide2_start"])


master_list <- (c(length_seq, length_c, length_c_prime, length_d, length_d_prime,pos_c, pos_c_prime, pos_d, pos_d_prime, dist_cd, dist_cpd, dist_cpdp, dist_dpc, guide1_start, guide2_start, guide1_length, guide2_length))

final_df <- as.data.frame(do.call(cbind, master_list))

cleaned_df <- replace_minus_one_with_na(final_df)

# If you have non-numeric columns, select only numeric ones first
numeric_df <- cleaned_df[, sapply(cleaned_df, is.numeric)]
cor_matrix <- cor(numeric_df, use = "complete.obs")
print(cor_matrix)

# Melt the correlation matrix
melted_cor <- melt(cor_matrix)

# Create heatmap
ggplot(melted_cor, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), size = 3) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1)) +
  theme_minimal() +
  labs(title = "Correlation Matrix", x = "", y = "")


# synthetic data

load(file = "results/intermediate/processed_info_box.RData")


pos_seq <- DNAStringSet(snodb_boxes[["DNA Sequence"]])
pos_seq

generate_random_sequences <- function(positive_seqs, n_random = 1000, 
                                      output_file = "random_negatives.fasta") {
  cat("Generating random sequences matching positive set properties...\n")
  positive_seqs <- snodb_boxes
  # Analyze positive sequences
  seq_lengths <- sapply(positive_seqs["DNA Sequence"], nchar)
  gc_contents <- letterFrequency(DNAStringSet(positive_seqs[["DNA Sequence"]]), letters = "GC", as.prob = TRUE)
  
  # Get distributions
  length_mean <- mean(seq_lengths)
  length_sd <- sd(seq_lengths)
  gc_mean <- mean(gc_contents)
  gc_sd <- sd(gc_contents)
  
  # random_seqs <- DNAStringSet()
  new_sequences_list <- list()
  new_sequences_name <- character(50)
  
  for (i in 1:n_random) {
    # Sample length from normal distribution
    target_length <- round(rnorm(1, mean = length_mean, sd = length_sd))
    target_length <- max(50, min(500, target_length))  # Constrain to reasonable range
    
    # Sample GC content
    target_gc <- rnorm(1, mean = gc_mean, sd = gc_sd)
    target_gc <- max(0.2, min(0.8, target_gc))  # Constrain between 20-80%
    
    # Generate sequence with target GC content
    n_gc <- round(target_length * target_gc)
    n_at <- target_length - n_gc
    
    bases <- sample(c(rep("G", floor(n_gc/2)), 
                      rep("C", ceiling(n_gc/2)),
                      rep("A", floor(n_at/2)), 
                      rep("T", ceiling(n_at/2))))
    
    # seq <- DNAString(paste(bases, collapse = ""))
    # names(seq) <- sprintf("random_neg_%04d|len=%d|gc=%.2f", i, target_length, target_gc)
    # random_seqs <- c(random_seqs, seq)
    seq <- DNAString(paste(bases, collapse = ""))
    
    new_sequences_name[i] <- sprintf("random_neg_%04d|len=%d|gc=%.2f", i, target_length, target_gc)
    
    new_sequences_list[[i]] <- seq 
  }
  
  random_seqs <- DNAStringSet(new_sequences_list)
  
  names(random_seqs) <- new_sequences_name
  
  # Write to file
  # writeXStringSet(random_seqs, filepath = output_file)
  
  cat(sprintf("Generated %d random sequences -> %s\n", n_random, output_file))
  
  return(random_seqs)
}
negatives <- generate_random_sequences(snodb_boxes, 500, "data/generated.fasta")
negatives

# given the lengths of the transcript look for the possible sequences

# define c,d,c',d' positions using SVM
# Multi-Output Position Predictor Model using SVM
# Predicts pos1, pos2, pos3, pos4 given a length

library(e1071)
library(caret)

# Sample training data (replace with your actual data)
set.seed(123)
train_data <- data.frame(
  length = unlist(length_seq),
  pos1 = unlist(pos_c),
  pos2 = unlist(pos_c_prime),
  pos3 = unlist(pos_d),
  pos4 = unlist(pos_d_prime)
)

# Display training data
print("Training Data:")
print(train_data)

# Method 1: SVM with RBF kernel (default - good for non-linear relationships)
cat("\n=== Training SVM models (RBF kernel) for each position ===\n")

# Train individual SVM models
svm_pos1 <- svm(pos1 ~ length, data = train_data, kernel = "radial")
svm_pos2 <- svm(pos2 ~ length, data = train_data, kernel = "radial")
svm_pos3 <- svm(pos3 ~ length, data = train_data, kernel = "radial")
svm_pos4 <- svm(pos4 ~ length, data = train_data, kernel = "radial")

# Function to predict all positions
predict_positions_svm <- function(length_value) {
  new_data <- data.frame(length = length_value)
  
  predictions <- data.frame(
    length = length_value,
    pos1 = predict(svm_pos1, new_data),
    pos2 = predict(svm_pos2, new_data),
    pos3 = predict(svm_pos3, new_data),
    pos4 = predict(svm_pos4, new_data)
  )
  
  return(predictions)
}

# Test predictions
cat("\n=== Testing Predictions (RBF kernel) ===\n")
test_lengths <- c(90, 120, 76)
for (len in test_lengths) {
  result <- predict_positions_svm(len)
  cat(sprintf("\nLength: %.0f -> pos1: %.2f, pos2: %.2f, pos3: %.2f, pos4: %.2f\n",
              result$length, result$pos1, result$pos2, result$pos3, result$pos4))
}

# Method 2: SVM with linear kernel (simpler, for linear relationships)
cat("\n\n=== Alternative: SVM with Linear Kernel ===\n")

svm_lin_pos1 <- svm(pos1 ~ length, data = train_data, kernel = "linear")
svm_lin_pos2 <- svm(pos2 ~ length, data = train_data, kernel = "linear")
svm_lin_pos3 <- svm(pos3 ~ length, data = train_data, kernel = "linear")
svm_lin_pos4 <- svm(pos4 ~ length, data = train_data, kernel = "linear")

predict_positions_svm_linear <- function(length_value) {
  new_data <- data.frame(length = length_value)
  
  predictions <- data.frame(
    length = length_value,
    pos1 = predict(svm_lin_pos1, new_data),
    pos2 = predict(svm_lin_pos2, new_data),
    pos3 = predict(svm_lin_pos3, new_data),
    pos4 = predict(svm_lin_pos4, new_data)
  )
  
  return(predictions)
}

# Test linear kernel
cat("\n=== Testing Predictions (Linear kernel) ===\n")
for (len in test_lengths) {
  result <- predict_positions_svm_linear(len)
  cat(sprintf("\nLength: %.0f -> pos1: %.2f, pos2: %.2f, pos3: %.2f, pos4: %.2f\n",
              result$length, result$pos1, result$pos2, result$pos3, result$pos4))
}



# Method 3: Tuned SVM with cross-validation (best performance)
cat("\n\n=== Training Tuned SVM with Cross-Validation ===\n")

# Define tuning parameters
tune_control <- tune.control(cross = 5)  # 5-fold cross-validation

# Tune SVM for each position
cat("Tuning pos1...\n")
tuned_pos1 <- tune(svm, pos1 ~ length, data = train_data,
                   ranges = list(epsilon = seq(0.1, 0.5, 0.1),
                                 cost = c(1, 10, 100)),
                   tunecontrol = tune_control)

cat("Tuning pos2...\n")
tuned_pos2 <- tune(svm, pos2 ~ length, data = train_data,
                   ranges = list(epsilon = seq(0.1, 0.5, 0.1),
                                 cost = c(1, 10, 100)),
                   tunecontrol = tune_control)

cat("Tuning pos3...\n")
tuned_pos3 <- tune(svm, pos3 ~ length, data = train_data,
                   ranges = list(epsilon = seq(0.1, 0.5, 0.1),
                                 cost = c(1, 10, 100)),
                   tunecontrol = tune_control)

cat("Tuning pos4...\n")
tuned_pos4 <- tune(svm, pos4 ~ length, data = train_data,
                   ranges = list(epsilon = seq(0.1, 0.5, 0.1),
                                 cost = c(1, 10, 100)),
                   tunecontrol = tune_control)

# Extract best models
best_svm_pos1 <- tuned_pos1$best.model
best_svm_pos2 <- tuned_pos2$best.model
best_svm_pos3 <- tuned_pos3$best.model
best_svm_pos4 <- tuned_pos4$best.model

# Function with tuned models
predict_positions_tuned <- function(length_value) {
  new_data <- data.frame(length = length_value)
  
  predictions <- data.frame(
    length = length_value,
    pos1 = predict(best_svm_pos1, new_data),
    pos2 = predict(best_svm_pos2, new_data),
    pos3 = predict(best_svm_pos3, new_data),
    pos4 = predict(best_svm_pos4, new_data)
  )
  
  return(predictions)
}

# Show best parameters
cat("\n=== Best Parameters Found ===\n")
cat("pos1:", sprintf("cost=%.0f, epsilon=%.2f\n", 
                     tuned_pos1$best.parameters$cost,
                     tuned_pos1$best.parameters$epsilon))
cat("pos2:", sprintf("cost=%.0f, epsilon=%.2f\n", 
                     tuned_pos2$best.parameters$cost,
                     tuned_pos2$best.parameters$epsilon))
cat("pos3:", sprintf("cost=%.0f, epsilon=%.2f\n", 
                     tuned_pos3$best.parameters$cost,
                     tuned_pos3$best.parameters$epsilon))
cat("pos4:", sprintf("cost=%.0f, epsilon=%.2f\n", 
                     tuned_pos4$best.parameters$cost,
                     tuned_pos4$best.parameters$epsilon))

# Test tuned models
cat("\n=== Testing Predictions (Tuned SVM) ===\n")
for (len in test_lengths) {
  result <- predict_positions_tuned(len)
  cat(sprintf("\nLength: %.0f -> pos1: %.2f, pos2: %.2f, pos3: %.2f, pos4: %.2f\n",
              result$length, result$pos1, result$pos2, result$pos3, result$pos4))
}

# Save models (optional)
# saveRDS(list(pos1=best_svm_pos1, pos2=best_svm_pos2, 
#              pos3=best_svm_pos3, pos4=best_svm_pos4), 
#         "svm_position_models.rds")

cat("\n=== Complete! ===\n")
cat("Use predict_positions_svm(length) for basic RBF SVM\n")
cat("Use predict_positions_svm_linear(length) for linear SVM\n")
cat("Use predict_positions_tuned(length) for optimized SVM (recommended)\n")

# library(car)
# library(corrplot)
# multipleCorrelation_to_one <- model.matrix( ~ length_seq + d_start ,numeric_df )
# define c,d,c',d' positions using RF

negative_set <- data.frame(sequence = )







# Multi-Output Position Predictor Model
# Predicts pos1, pos2, pos3, pos4 given a length

library(caret)
library(randomForest)

# Sample training data (replace with your actual data) -> have to calculate each on a different stage? 3 model length ->  d -> 
set.seed(123)
train_data <- data.frame(
  length = sapply(snodb_boxes["DNA Sequence"], nchar)
  # posc = c(10, 15, 20, 25, 30, 35, 40, 45, 50),
  # poscp = c(30, 45, 60, 75, 90, 105, 120, 135, 150),
  # posd = c(60, 90, 120, 150, 180, 210, 240, 270, 300),
  # posdp = c(90, 135, 180, 225, 270, 315, 360, 405, 450)
)

# Display training data
print("Training Data:")
print(train_data)

# Method 1: Separate models for each position (recommended)
cat("\n=== Training separate models for each position ===\n")

# Train individual models
model_pos1 <- randomForest(pos1 ~ length, data = train_data, ntree = 100)
model_pos2 <- randomForest(pos2 ~ length, data = train_data, ntree = 100)
model_pos3 <- randomForest(pos3 ~ length, data = train_data, ntree = 100)
model_pos4 <- randomForest(pos4 ~ length, data = train_data, ntree = 100)

# Function to predict all positions
predict_positions <- function(length_value) {
  new_data <- data.frame(length = length_value)
  
  predictions <- data.frame(
    length = length_value,
    pos1 = predict(model_pos1, new_data),
    pos2 = predict(model_pos2, new_data),
    pos3 = predict(model_pos3, new_data),
    pos4 = predict(model_pos4, new_data)
  )
  
  return(predictions)
}

# Test predictions
cat("\n=== Testing Predictions ===\n")
test_lengths <- c(175, 275, 425)
for (len in test_lengths) {
  result <- predict_positions(len)
  cat(sprintf("\nLength: %.0f -> pos1: %.2f, pos2: %.2f, pos3: %.2f, pos4: %.2f\n",
              result$length, result$pos1, result$pos2, result$pos3, result$pos4))
}

# Method 2: Linear regression (simpler alternative)
cat("\n\n=== Alternative: Linear Regression Models ===\n")

lm_pos1 <- lm(pos1 ~ length, data = train_data)
lm_pos2 <- lm(pos2 ~ length, data = train_data)
lm_pos3 <- lm(pos3 ~ length, data = train_data)
lm_pos4 <- lm(pos4 ~ length, data = train_data)

# Function for linear model predictions
predict_positions_lm <- function(length_value) {
  new_data <- data.frame(length = length_value)
  
  predictions <- data.frame(
    length = length_value,
    pos1 = predict(lm_pos1, new_data),
    pos2 = predict(lm_pos2, new_data),
    pos3 = predict(lm_pos3, new_data),
    pos4 = predict(lm_pos4, new_data)
  )
  
  return(predictions)
}

# Model evaluation
cat("\n=== Model Performance (R-squared) ===\n")
cat(sprintf("pos1 R²: %.4f\n", summary(lm_pos1)$r.squared))
cat(sprintf("pos2 R²: %.4f\n", summary(lm_pos2)$r.squared))
cat(sprintf("pos3 R²: %.4f\n", summary(lm_pos3)$r.squared))
cat(sprintf("pos4 R²: %.4f\n", summary(lm_pos4)$r.squared))

# Save models (optional)
# saveRDS(list(pos1=model_pos1, pos2=model_pos2, pos3=model_pos3, pos4=model_pos4), 
#         "position_models.rds")

# To load later:
# models <- readRDS("position_models.rds")

cat("\n=== Complete! ===\n")
cat("Use predict_positions(length) for Random Forest predictions\n")
cat("Use predict_positions_lm(length) for Linear Model predictions\n")



#genome's intergenic and intronic regions
# i could use a model like a vae (the discriminator part) to analyse the transcript in batch of lengths =~2*max avrg length of the lonest snoRNAs, so that the start is at i and the next is i+i/2 si they overlap




