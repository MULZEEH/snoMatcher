#==============================================================
#             GENERATING GUIDE SCORES FROM SAMPLE
#==============================================================
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
      scores = "results/intermediate/scores.RData"
     
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

negatives <- read.csv(get_input("negatives"))
# what is seq score?? rRNA su negatives e' sempre lo stesso???
negatives <- negatives %>% dplyr::select(unique_id, seq_score, c_box_seq, d_prime_box_seq, c_prime_box_seq, d_box_seq,c_d_prime_dist, d_prime_c_prime_dist, c_prime_d_dist, c_d_dist,)
negatives$guide_score <- guides_scores[negatives$seq_score]
negatives$c_box_score <- cbox_scores[negatives$c_box_seq]
negatives$cprime_box_score <- c_prime_box_scores[negatives$c_prime_box_seq]
negatives$d_box_score <- dbox_scores[negatives$d_box_seq]
negatives$dprime_box_score <- d_prime_box_scores[negatives$d_prime_box_seq]
negatives <- negatives[sample(1:nrow(negatives), 600),]
negatives$snoRNA <- F
negatives <- negatives %>% dplyr::rename(dist_d_prime_c_prime = d_prime_c_prime_dist, dist_c_prime_d = c_prime_d_dist, dist_c_d_prime = c_d_prime_dist, dist_c_d = c_d_dist)
snorna_machine_learning <- rbind(snodb_data_guide %>% dplyr::select(c_box_score, dprime_box_score, cprime_box_score, d_box_score, dist_c_d_prime, dist_d_prime_c_prime, dist_c_prime_d, dist_c_d, guide_score, snoRNA), negatives %>% dplyr::select(c_box_score, dprime_box_score, cprime_box_score, d_box_score, dist_c_d_prime, dist_d_prime_c_prime, dist_c_prime_d, dist_c_d, guide_score, snoRNA))

# save(snodb_data_guide, negatives, snorna_machine_learning, file = "snorna_machine_learning.RData")

# a snorna has a shorter cbox, i need to remove it for now -> what snorna?
snorna_machine_learning <- na.omit(snorna_machine_learning)
snorna_machine_learning$snoRNA <- as.factor(snorna_machine_learning$snoRNA)

save(snodb_data_guide, negatives, snorna_m