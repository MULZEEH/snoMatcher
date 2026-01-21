#==============================================================
#             GENERATING GUIDE SCORES FROM SAMPLE
#==============================================================
#this script aim to generate negative sample to give to the model.
# The type of sample used are in order:
# Firstly the hard negatives found by a first execution of the model, so the false positives ones
# then we have a set of scrambled generated ones
# then we have a set of defined by the dimension and maybe generated who knows

# In the end the script return the whole  

library(Biostrings)
library(data.table)
library(dplyr)
library(stringr)
library(writexl)
library(GenomicFeatures)
library(readxl)
library(stringi)
library(seqinr)
library(parallel)
library(ggseqlogo)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
# library(seqLogo)
library(ape)
library(tidyr)
library(patchwork)
library(tidyverse)


# Check if running with Snakemake or in RStudio
# Need to clear the environment first -> TODO(FIX)
rm(snakemake, envir = .GlobalEnv)
if (!exists("snakemake")) {
  # need to change the static hardcoded setwd
  setwd("C:/Users/Marco/RProject/snoMatcher/")
  # Create mock snakemake object for testing in Rstudio
  snakemake <- list(
    input = list(
      info_box = "results/intermediate/processed_info_box.RData",
      guide = "results/intermediate/guide.RData",
      negatives = "data/raw/gene_bodies_ENSG00000079974.19_site_6437_.csv",
      scores = "results/intermediate/scores.RData",
      mm_hits = "data/raw/snoRNA_hits_mm.csv",
      synthetic_threshold = 200
     
    ),
    output = list(
      
      snorna_machine_learning = "results/intermediate/snorna_machine_learning.RData"
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
# somewhat a graphic analysis of the negatives?
  # Helper functions for real snakemake object
  get_input <- function(name) snakemake@input[[name]]
  get_output <- function(name) snakemake@output[[name]]
  get_config <- function(name) snakemake@config[[name]]
  # get_threads <- function() snakemake@threads
  ("snakemake execution")
}

load(file = get_input("info_box"))
load(file = get_input("scores"))
load(file = get_input("guide"))

# mouse integration in the positive dataset

mm_snorna <- read.csv(get_input("mm_hits"))

# read the negatves that were computed previously
negatives <- read.csv(get_input("negatives"))

anno_function <- function (rna_seq){
  ("annotated")
}
# feature extraction function aims to return the dataframe with the correct feature to pass in the learning model
# as input it requires the 'rna_seq' df (autoesplicativo), 
# 'type': of the dataframe (boolean) that describe if its the positive value or negative one

feature_extraction <- function (rna_seq, anno){
  if(!anno){
    anno_function(rna_seq)
  }
  rna_seq <- rna_seq %>% dplyr::select(unique_id, seq_score, c_box_seq, d_prime_box_seq, c_prime_box_seq, d_box_seq,c_d_prime_dist, d_prime_c_prime_dist, c_prime_d_dist, c_d_dist)
  rna_seq$guide_score <- guides_scores[rna_seq$seq_score]
  rna_seq$c_box_score <- cbox_scores[rna_seq$c_box_seq]
  rna_seq$cprime_box_score <- c_prime_box_scores[rna_seq$c_prime_box_seq]
  rna_seq$d_box_score <- dbox_scores[rna_seq$d_box_seq]
  rna_seq$dprime_box_score <- d_prime_box_scores[rna_seq$d_prime_box_seq]
  rna_seq$snoRNA <- type
  rna_seq <- rna_seq %>% dplyr::rename(dist_d_prime_c_prime = d_prime_c_prime_dist, dist_c_prime_d = c_prime_d_dist, dist_c_d_prime = c_d_prime_dist, dist_c_d = c_d_dist)

  return rna_seq
    
}

mm_snorna <- feature_extraction(mm_snorna, T)
negatives <- feature_extraction(negatives, F)

#===
# ALTERNATIVE METHODS
#===

#function to scramble the sequence keeping a safe zone (the guide that we do NOT want to change)
scramble_sequence <- function(df) {
  
  # Function to scramble a single sequence
  scramble_one <- function(sequence, start, pace) {
    # Convert sequence to character vector
    seq_chars <- strsplit(sequence, "")[[1]]
    
    # Define the protected region (start to start+pace-1)
    protected_start <- start
    protected_end <- start + pace - 1 #?
    
    # Get indices to scramble (everything except the protected region)
    all_indices <- 1:length(seq_chars)
    protected_indices <- protected_start:protected_end
    scramble_indices <- setdiff(all_indices, protected_indices)
    
    # Scramble only the non-protected positions
    if (length(scramble_indices) > 0) {
      seq_chars[scramble_indices] <- sample(seq_chars[scramble_indices])
    }
    
    # Convert back to string
    return(paste(seq_chars, collapse = ""))
  }
  
  # Apply to each row of the dataframe
  df$scrambled_sequence <- mapply(scramble_one, 
                                  df$`DNA Sequence`, 
                                  df$guide_start, 
                                  df$guide_width)
  
  return(df)
}

# This function retrieves the 
get_sequence_interval <- function(seq, start_index, interval_length) {
  
  # The ending position is (start_index + length - 1)
  end_index <- start_index + interval_length - 1
  
  extracted_portion <- substr(seq, start = start_index, stop = end_index)
  
  return(extracted_portion)
}

# take a subset (some will automaticcally be excluded thanks to the NA of box scoring system)
# first i select only the sequence and the position to the boxes so that i can scamble without interfering with the guide seq

#need to first join (tbj)
guide_start_tbj <- snodb_boxes[, c('snoDB ID','guide1_start', 'guide1_width')]
guide_start_tbj <- guide_start_tbj[guide_start_tbj$guide1_start != -1,]

tmp_tbj <- snodb_data_guide[, c("snoDB ID",
                                "DNA Sequence",
                                "c_start",
                                "d_start",
                                "c_prime_start",
                                "d_prime_start",
                                "dist_c_d_prime",
                                "dist_c_prime_d",
                                "dist_d_prime_c_prime",
                                "dist_c_d",
                                "guide_score"
                                )]

tmp_joined <- merge(x = guide_start_tbj,
                    y = tmp_tbj,
                    by = "snoDB ID",
                    all = FALSE)

tmp <- tmp_joined %>% dplyr::rename(guide_start = guide1_start, guide_width = guide1_width)

#function that given a dna seq, it returns the same seq scrambled but in the range from start_pos to end_pos
# PROBLEM IF SCRAMBLED IT KEEPS THE POSITIONS SO DISTANCES AND SUCH ARE INVARIANT

# now with guide2:

guide_start_tbj <- snodb_boxes[, c('snoDB ID','guide2_start', 'guide2_width')]
guide_start_tbj <- guide_start_tbj[guide_start_tbj$guide2_start != -1,]

tmp_joined <- merge(x = guide_start_tbj,
                    y = tmp_tbj,
                    by = "snoDB ID",
                    all = FALSE)

tmp2 <- tmp_joined %>% dplyr::rename(guide_start = guide2_start, guide_width = guide2_width)


# tmp$`snoDB ID` <- paste(tmp$`snoDB ID`)
# maybe rename even if the name will be removed later

#since the analysis here is based on the scramble that introduce randomness i 
# will put a constraint to do the analysis n times until the number of synthetic 
# data are more then a defined threshold
alternative_negatives <- list()
repeat{
  alternative <- scramble_sequence(rbind(tmp, tmp2)) # <- here is to avoid getting the scrambled seq
  
  # Now i have to retrieve from the Scrambled Sequence the "new" boxes
  alternative$c_seq <- get_sequence_interval(alternative$scrambled_sequence, 
                                                       tmp$c_start, 7)
  alternative$c_prime_seq <- get_sequence_interval(alternative$scrambled_sequence, 
                                                             tmp$c_prime_start, 7)
  alternative$d_seq <- get_sequence_interval(alternative$scrambled_sequence, 
                                                       tmp$d_start, 4)
  alternative$d_prime_seq <- get_sequence_interval(alternative$scrambled_sequence, 
                                                             tmp$d_prime_start, 4)
  
  
  alternative <- feature_extraction(alternative, F)
  
  alternative_negatives <- unique(rbind(alternative, alternative_negatives))
  # print(nrow(unique(alternative_negatives)))
  number_of_synth <- nrow(alternative_negatives)
  print(number_of_synth)
  if( number_of_synth > get_input("synthetic_threshold") ){
    break
  }
}


#defining the list of dataframe that will be binded:
list_to_be_binded <- list(snodb_data_guide %>% dplyr::select(c_box_score, dprime_box_score, cprime_box_score, d_box_score, dist_c_d_prime, dist_d_prime_c_prime, dist_c_prime_d, dist_c_d, guide_score, snoRNA),
                          mm_snorna %>% dplyr::select(c_box_score, dprime_box_score, cprime_box_score, d_box_score, dist_c_d_prime, dist_d_prime_c_prime, dist_c_prime_d, dist_c_d, guide_score, snoRNA),
                          negatives %>% dplyr::select(c_box_score, dprime_box_score, cprime_box_score, d_box_score, dist_c_d_prime, dist_d_prime_c_prime, dist_c_prime_d, dist_c_d, guide_score, snoRNA),
                          alternative_negatives %>% dplyr::select(c_box_score, dprime_box_score, cprime_box_score, d_box_score, dist_c_d_prime, dist_d_prime_c_prime, dist_c_prime_d, dist_c_d, guide_score, snoRNA))

#binding the negatives with the positives ones
snorna_machine_learning <- rbindlist(list_to_be_binded)

list <- list(snodb_data_guide %>% dplyr::select(c_box_score, dprime_box_score, cprime_box_score, d_box_score, dist_c_d_prime, dist_d_prime_c_prime, dist_c_prime_d, dist_c_d, guide_score, snoRNA),
                          negatives %>% dplyr::select(c_box_score, dprime_box_score, cprime_box_score, d_box_score, dist_c_d_prime, dist_d_prime_c_prime, dist_c_prime_d, dist_c_d, guide_score, snoRNA)
)
#binding the negatives with the positives ones
small_snorna_ml <- rbindlist(list)
# save(snodb_data_guide, negatives, snorna_machine_learning, file = "snorna_machine_learning.RData")
snorna_ml_scr <- snorna_machine_learning[sample(1:nrow(snorna_machine_learning), nrow(snorna_machine_learning)),]

# a snorna has a shorter cbox, i need to remove it for now
snorna_machine_learning <- na.omit(snorna_machine_learning)
snorna_machine_learning$snoRNA <- as.factor(snorna_machine_learning$snoRNA)

save(snodb_data_guide, negatives, snorna_machine_learning, file = get_output("snorna_machine_learning"))

