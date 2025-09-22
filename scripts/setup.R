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
library(seqLogo)
library(ape)
library(tidyr)

# Check if running in Snakemake or RStudio -> need to update to what i have
if (!exists("snakemake")) {
  # Create mock snakemake object for testing
  snakemake <- list(
    input = list(
      snodb = "data/snoDB_data.xlsx",
      interactions = "data/interactions.xlsx",
      boxes = "data/boxes.tsv"
    ),
    output = list(
      results = "results/analysis.rds",
      plots = "results/plots.pdf"
    ),
    config = list(
      threads = 4,
      significance_threshold = 0.05
    ),
    threads = 4
  )
}


snodb <- read_xlsx(snakemake@input[["sno_db"]])
met_sites <- read_xlsx(snakemake@input[["met_sites"]])
snodb_boxes <- fread(snakemake@input[["sno_boxes"]])

snodb_boxes <- snodb_boxes[snodb_boxes$c_start != -1 & snodb_boxes$c_prime_start != -1 & snodb_boxes$d_start != -1 & snodb_boxes$d_prime_start != -1,]
snodb_boxes <- merge(snodb[, c("snoDB ID", "Symbol", "DNA Sequence")], snodb_boxes, by.x = "snoDB ID", by.y = "snoDB_id")
snodb_boxes <- snodb_boxes[snodb_boxes$`snoDB ID` %in% met_sites$snoDB_id[met_sites$Status == "validated"],]

snoRNA_seq <- DNAStringSet(snodb_boxes$`DNA Sequence`)
names(snoRNA_seq) <- snodb_boxes$`snoDB ID`

met_sites <- met_sites[met_sites$snoDB_id %in% snodb_boxes$`snoDB ID` & met_sites$Status == "validated",]

snodb_boxes$guide1_width <- str_count(snodb_boxes$guide1_seq)
snodb_boxes$guide2_width <- str_count(snodb_boxes$guide2_seq)

# c <- c(snodb_boxes$c_seq[str_count(snodb_boxes$c_seq) == 7])
# 
# c_prime <- c(snodb_boxes$c_prime_seq)
# 
# d <- c(snodb_boxes$d_seq)
# 
# d_prime <- c(snodb_boxes$d_prime_seq)

save(met_sites, snodb_boxes, file = snakemake@output[["info_boxes"]])
