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
      info_box = "results/intermediate/processed_info_box.RData"
      #  add the fasta of the hs_rna mao
    ),
    output = list(
      
      guide = "results/intermediate/guide.RData"
      
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

load(file = get_input("info_box"))

#===
#PAIRING DISTRIBUTIONS -> guide distribution could also have its own file
#===

# why do i have 2 type(?) of sequences?
snodb_data <- snodb_boxes %>% dplyr::select(-c("guide1_start", "guide2_start")) %>% pivot_longer(c(guide1_seq, guide2_seq), names_to = "guide", values_to = "guide_seq") %>% filter(guide_seq != "") %>% mutate(guide_seq = str_sub(guide_seq, -10, -2))

rRNA_seq<- readDNAStringSet("data/raw/hs_rRNA_hebras_processed.fasta")
met_sites$rRNA_seq <- as.character(subseq(rRNA_seq[met_sites$rRNA], 
                                          met_sites$`Pos snoRNABase`-3, 
                                          met_sites$`Pos snoRNABase`+5))

# detection guide region 9 (+1 at the beginning) nt ???
snodb_data_guide <- merge(snodb_data, 
                          met_sites[, c("snoDB_id", "Symbol", 
                                        "rRNA", "Pos snoRNABase", 
                                        "rRNA_seq")], 
                          by.x = c("snoDB ID", "Symbol"), 
                          by.y = c("snoDB_id", "Symbol"))

# guide evaluation
snorna_match_seq <- as.data.frame(tstrsplit(as.character(reverse(snodb_data_guide$guide_seq)),
                                            "",
                                            fixed=T))

rRNA_match_seq <- as.data.frame(tstrsplit(as.character(snodb_data_guide$rRNA_seq),
                                          "",
                                          fixed=T))

for(j in 1:9){
  snorna_match_seq[,j] <- paste0(as.character(snorna_match_seq[,j]),
                                 as.character(rRNA_match_seq[,j]))
  snorna_match_seq[,j] <- str_replace_all(snorna_match_seq[,j],
                                          "AT|TA|GC|CG",
                                          "V")
  snorna_match_seq[,j] <- str_replace_all(snorna_match_seq[,j],
                                          "GT|TG",
                                          "U")
  
  non_matching <- !grepl("V|U", snorna_match_seq[,j])
  snorna_match_seq[non_matching,j] <- "M"
}

# snodb_data_guide$guide_match <- paste0(snorna_match_seq[,1], snorna_match_seq[,2], snorna_match_seq[,3],snorna_match_seq[,4],snorna_match_seq[,5], snorna_match_seq[,6], snorna_match_seq[,7], snorna_match_seq[,8],snorna_match_seq[,9])
snodb_data_guide$guide_match <- do.call(paste0, snorna_match_seq[, 1:9])
snodb_data_guide <- snodb_data_guide[str_split_fixed(snodb_data_guide$guide_match, "", 9)[,4] == "V",]
snodb_data_guide <- snodb_data_guide[str_count(snodb_data_guide$guide_match, "M") <2,]

snodb_data_guide <- snodb_data_guide %>% dplyr::select(-c(rRNA, `Pos snoRNABase`)) %>% distinct()

pos_pairing <-  as.data.frame(str_split_fixed(snodb_data_guide$guide_match, pattern = "", 9))
names(pos_pairing) <- c("i-3", "i-2", "i-1", "i", "i+1", "i+2", "i+3", "i+4", "i+5")


#===
#ACTUAL GENERATING
#===
# Sample data matrix
data_matrix <- as.matrix(pos_pairing)

# Find unique symbols and positions
unique_symbols <- unique(as.vector(data_matrix))
positions <- 1:9

# Calculate the frequency of each symbol in each position using apply()
frequency_matrix <- sapply(positions, function(pos) {
  col_data <- data_matrix[, pos]
  table(factor(col_data, levels = unique_symbols))
})

# Convert the result to a matrix, setting row and column names
frequency_matrix <- as.matrix(frequency_matrix)
frequency_matrix <- frequency_matrix/nrow(pos_pairing)

possible_guides <- utils::combn(rep(c("V", "U", "M"), 9),9)
possible_guides = unique(apply(possible_guides, 2, paste0, collapse=""))
dist.prova = sapply(possible_guides, function(x) adist(x, "VVVVVVVVV"))
possible_guides <- possible_guides[dist.prova < 4]
# save(possible_guides, file = "snoDB_analysis_validated/possible_guides.RData")

# Function to calculate the score of a string based on the frequency matrix
calculate_score <- function(input_string, frequency_matrix) {
  input_chars <- strsplit(input_string, "")[[1]]
  row_indices <- match(input_chars, rownames(frequency_matrix))
  col_indices <- seq_along(input_chars)
  
  score <- 0
  for (i in 1:length(input_chars)) {
    score <- score + frequency_matrix[row_indices[i], col_indices[i]]
  }
  return(score)
}

# Calculate the score for each element in the input vector
guides_scores <- sapply(possible_guides, calculate_score, frequency_matrix = frequency_matrix)
# guide
guides_scores_prop <- guides_scores/max(guides_scores)
guides_scores_prop <- setNames(guides_scores_prop, possible_guides)
names(guides_scores) <- str_split_fixed(names(guides_scores), "[.]",2)[,1]
# save(guides_scores,guides_scores_prop, file = "snoDB_analysis_validated/guides_scores.RData")

# guide scores of my results
snodb_data_guide$guide_score <- guides_scores[snodb_data_guide$guide_match]
# minimum is 0.69

snodb_data_guide$snoRNA <- T

save(guides_scores, guides_scores_prop, snodb_data_guide, file = get_output("guide"))

# plot con dati di analisi rRNA_seq on row 1295 (mold)
# ggplot(hits_final, aes(guide_score))+
#   geom_histogram(fill = "#434384")+
#   theme_bw()
# 
# ggsave("snoDB_analysis_validated/guide_score_dist.pdf", device = cairo_pdf, height = 3, width = 4)

#==============================================================
#                   HELIX SOMEHTING
#==============================================================
start_seq <- subseq(snoRNA_seq, pmax(1, snodb_boxes$c_start-6), snodb_boxes$c_start+7)
end_seq <- subseq(snoRNA_seq, snodb_boxes$d_start+1, pmin(snodb_boxes$d_start+11,width(snoRNA_seq)))

# start_seq <- subseq(snoRNA_seq, 1, snodb_boxes$c_start+7)
# end_seq <- subseq(snoRNA_seq, snodb_boxes$d_start+1, width(snoRNA_seq))




rev_end_seq <- reverse(end_seq)

# what ar eu doing?
score_system <- matrix(c(-3,-3,-3,+2,
                         -3,-3,+3,-3,
                         -3,+3,-3,+1,
                         +2,-3,+1,-3),
                       nrow=4, ncol = 4,
                       dimnames = list(c("A", "C", "G", "T"), c("A", "C","G","T")))

pwa <- pairwiseAlignment(start_seq, rev_end_seq, type = "local", substitutionMatrix = score_system, gapOpening = 4, gapExtension = 2)
align_df <- data.table("transcript"= names(snoRNA_seq), "start_align" = as.character(alignedPattern(pwa)), "end_align" = as.character(reverse(alignedSubject(pwa))), "score" = score(pwa))

#creo file senza -
align_df$start_align <- str_replace_all(align_df$start_align, "-", "")
align_df$end_align <- str_replace_all(align_df$end_align, "-", "")

table(str_count(align_df$start_align))
align_df$length_align <- str_count(align_df$start_align)

ggplot(align_df, aes(length_align))+
  geom_bar(fill = "#434384")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))+
  geom_text(stat="count", aes(x=length_align, label = stat(count)),position = position_dodge(width = 0.9), vjust = -0.3, size = 3)

ggsave("snoDB_analysis_validated/terminal_helix_with_boxes_align.pdf", device = cairo_pdf, width = 4, height = 3)


align_df <- align_df[str_count(align_df$start_align) >2]
# 215 su 228



# match pattern per trovare posizione allineamento
require(doParallel)
require(foreach)

ncores <- detectCores()-1
cluster_ext <- makeCluster(ncores, type = "SOCK")
registerDoParallel(cl = cluster_ext)

pos_align_df<-NULL

pos_align_df<-foreach(transcript = align_df$transcript,.combine='rbind',.packages = c("data.table","Biostrings")) %dopar% {
  data.table("transcript" = transcript, tail(as.data.table(vmatchPattern(align_df$start_align[align_df$transcript == transcript], start_seq[transcript])),1), as.data.table(vmatchPattern(align_df$end_align[align_df$transcript == transcript], end_seq[transcript]))[1,])
}
stopCluster(cluster_ext)


#chiamo le colonne nel modo giusto che se no mi sbaglio
names(pos_align_df)[c(4:6,9:11)] <- c("align_5_start", "align_5_end", "align_5_width", "align_3_start", "align_3_end", "align_3_width")

# aggiungo risultati match pattern a tabella principale ----
align_df <- merge(align_df, pos_align_df[,c(1, 4,5,6,9,10,11)], by = "transcript")

align_df <- merge(align_df, snodb_boxes[,c(1, 6,8)], by.x = "transcript", by.y = "snoDB ID")

# trovo posizione reale allineamento in 3'

align_df$align_3_start <- align_df$d_start +1 + align_df$align_3_start-1
align_df$align_3_end <- align_df$d_start +1 + align_df$align_3_end-1

# trovo posizione reale allineamento in 5'

align_df$align_5_start <- pmax(1, align_df$c_start-6) + align_df$align_5_start -1
align_df$align_5_end <- pmax(1, align_df$c_start-6) + align_df$align_5_end -1

# guardo distanza massima di inizio allineamento
align_df$dist_c <- align_df$c_start-align_df$align_5_end+1
align_df$dist_d <- align_df$align_3_start-(align_df$d_start+4)

write_xlsx(align_df, "snoDB_analysis_validated/alignments_terminal_helix.xlsx")

align_boxes <- align_df

#===
#SECOND METHOD (box excluded, on rfam no interaction between boxes is shown)
#===

start_seq <- subseq(snoRNA_seq, pmax(1, snodb_boxes$c_start-6), snodb_boxes$c_start)
end_seq <- subseq(snoRNA_seq, snodb_boxes$d_start+5, pmin(snodb_boxes$d_start+11,width(snoRNA_seq)))

# start_seq <- subseq(snoRNA_seq, 1, snodb_boxes$c_start)
# end_seq <- subseq(snoRNA_seq, snodb_boxes$d_start+5, width(snoRNA_seq))

rev_end_seq <- reverse(end_seq)


score_system <- matrix(c(-3,-3,-3,+2,
                         -3,-3,+3,-3,
                         -3,+3,-3,+1,
                         +2,-3,+1,-3),
                       nrow=4, ncol = 4,
                       dimnames = list(c("A", "C", "G", "T"), c("A", "C","G","T")))

pwa <- pairwiseAlignment(start_seq, rev_end_seq, type = "local", substitutionMatrix = score_system, gapOpening = 4, gapExtension = 2)
align_df <- data.table("transcript"= names(snoRNA_seq), "start_align" = as.character(alignedPattern(pwa)), "end_align" = as.character(reverse(alignedSubject(pwa))), "score" = score(pwa))

#creo file senza -
align_df$start_align <- str_replace_all(align_df$start_align, "-", "")
align_df$end_align <- str_replace_all(align_df$end_align, "-", "")

table(str_count(align_df$start_align))
align_df$length_align <- str_count(align_df$start_align)

ggplot(align_df, aes(length_align))+
  geom_bar(fill = "#434384")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_x_continuous(breaks = c(seq(0,20,2)))+
  geom_text(stat="count", aes(x=length_align, label = stat(count)),position = position_dodge(width = 0.9), vjust = -0.3, size = 3)

ggsave("snoDB_analysis_validated/terminal_helix_without_boxes_align.pdf", device = cairo_pdf, width = 4, height = 3)

ggplot(align_df, aes(score))+
  geom_bar(fill = "#434384")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  scale_x_continuous(breaks = c(seq(0,20,2)))+
  geom_text(stat="count", aes(x=score, label = stat(count)),position = position_dodge(width = 0.9), vjust = -0.3, size = 3)

ggsave("snoDB_analysis_validated/score_terminal_helix_without_boxes_align.pdf", device = cairo_pdf, width = 4, height = 3)


align_df <- align_df[str_count(align_df$start_align) >2]
# 187 su 228



# match pattern per trovare posizione allineamento
require(doParallel)
require(foreach)

ncores <- detectCores()-1
cluster_ext <- makeCluster(ncores, type = "SOCK")
registerDoParallel(cl = cluster_ext)

pos_align_df<-NULL

pos_align_df<-foreach(transcript = align_df$transcript,.combine='rbind',.packages = c("data.table","Biostrings")) %dopar% {
  data.table("transcript" = transcript, tail(as.data.table(vmatchPattern(align_df$start_align[align_df$transcript == transcript], start_seq[transcript])),1), as.data.table(vmatchPattern(align_df$end_align[align_df$transcript == transcript], end_seq[transcript]))[1,])
}
stopCluster(cluster_ext)


#chiamo le colonne nel modo giusto che se no mi sbaglio
names(pos_align_df)[c(4:6,9:11)] <- c("align_5_start", "align_5_end", "align_5_width", "align_3_start", "align_3_end", "align_3_width")

# aggiungo risultati match pattern a tabella principale ----
align_df <- merge(align_df, pos_align_df[,c(1, 4,5,6,9,10,11)], by = "transcript")

align_df <- merge(align_df, snodb_boxes[,c(1, 6,8)], by.x = "transcript", by.y = "snoDB ID")

# trovo posizione reale allineamento in 3'

align_df$align_3_start <- align_df$d_start +5 + align_df$align_3_start-1
align_df$align_3_end <- align_df$d_start +5 + align_df$align_3_end-1

# trovo posizione reale allineamento in 5'

align_df$align_5_start <- pmax(1, align_df$c_start-6) + align_df$align_5_start -1
align_df$align_5_end <- pmax(1, align_df$c_start-6) + align_df$align_5_end -1

# guardo distanza massima di inizio allineamento
align_df$dist_c <- align_df$c_start-align_df$align_5_end+1
align_df$dist_d <- align_df$align_3_start-(align_df$d_start+4)

write_xlsx(align_df, "snoDB_analysis_validated/alignments_terminal_helix.xlsx")

