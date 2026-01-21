"""
# In this script we compute the scores for each box (C, C', D, D') based on
# Position Frequency Matrices (PFMs) derived from snoDB sequences.
# We also calculate distances between boxes to generate relevant plots.

"""
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
library(fitdistrplus)


# Check if running with Snakemake or in RStudio
rm(snakemake, envir = .GlobalEnv)
if (!exists("snakemake")) {
  # need to change the static hardcoded setwd
  setwd("C:/Users/Marco/RProject/snoMatcher/")
  # Create mock snakemake object for testing in Rstudio
  snakemake <- list(
    input = list(
      info_box = "results/intermediate/info_box.RData"
    ),
    output = list(

      possible_box = "results/intermediate/possible_boxes.RData",
      scores = "results/intermediate/scores.RData",
      dist_fit = "results/intermediate/dist_fit.RData",
      
      #PLOTS
      box_mismatch_distribution = "results/plots/box_mismatch_distribution.pdf",
      #violin
      all_distance = "results/plots/all_distance.pdf",
      #"histogram" density
      all_distance_distributions = "results/plots/all_distance_distributions.pdf",
      motif_score_distribution = "results/plots/motif_score_snoDB_distribution.pdf",
      up_motif_score_distribution = "results/plots/up_motif_score_distribution.pdf",
      down_motif_score_distribution = "results/plots/down_motif_score_distribution.pdf",
      
      #TABLES
      pfm_cbox_snoDb = "results/tables/pfm_cbox_snoDb.csv",
      pfm_cbox_prime_snoDb = "results/tables/pfm_cbox_prime_snoDb.csv",
      pfm_dbox_snoDb = "results/tables/pfm_dbox_snoDb.csv",
      pfm_dbox_prime_snoDb = "results/tables/pfm_dbox_prime_snoDb.csv"
      
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

#=== LOADING DEPENDANCY DATA ===
load(file = get_input("info_box"))

# possibility of using not only the distance but the some form of correlation between the distance and the lenthg of the snoRNA
c <- c(snodb_boxes$c_seq[str_count(snodb_boxes$c_seq) == 7])
c_prime <- c(snodb_boxes$c_prime_seq)
d <- c(snodb_boxes$d_seq)
d_prime <- c(snodb_boxes$d_prime_seq)

#=== DISTANCE for SCOREs ===
snodb_boxes$dist_d_prime_c_prime <- snodb_boxes$c_prime_start - (snodb_boxes$d_prime_start+3)
snodb_boxes$dist_c_d_prime <- snodb_boxes$d_prime_start - 
  (snodb_boxes$c_start+6)
snodb_boxes$dist_c_prime_d <- snodb_boxes$d_start-
  (snodb_boxes$c_prime_start+6)
snodb_boxes$dist_c_d <- snodb_boxes$d_start - (snodb_boxes$c_start+6)


# FOR EACH OF THE BOX I COMPUTE THE SCORES -> TODO(find a way to loop)
#==============================================================
#                C-BOX (and C') SCOREs
#==============================================================

#=== BOX SCOREs ===
# computing C and C' box scores starting from PFMs
#==================

# generate all possible 7nt combinations with A,T,C,G
prova <-utils::combn(rep(c("C", "T", "G", "A"), 7),7)
prova.str = unique(apply(prova, 2, paste0, collapse=""))
dist.prova1 = sapply(prova.str, function(x) adist(x, "ATGATGA"))
dist.prova2 = sapply(prova.str, function(x) adist(x, "GTGATGA"))

possible_cboxes1 <- prova.str[dist.prova1 <= 4]
possible_cboxes2 <- prova.str[dist.prova2 <= 4]

possible_cboxes <- union(possible_cboxes1, possible_cboxes2)

#=== PFM ===
pfm_c <- consensusMatrix(c, as.prob = T)
pfm_c_prime <- consensusMatrix(c_prime, as.prob = T)

#=== ACTUAL SCOREs ===
cbox_scores <- sapply(possible_cboxes, PWMscoreStartingAt, pwm = as.matrix(pfm_c))
cbox_scores_prop <- cbox_scores/max(cbox_scores)

c_prime_box_scores <- sapply(possible_cboxes, PWMscoreStartingAt, pwm = as.matrix(pfm_c_prime))
c_prime_box_scores_prop <- c_prime_box_scores/max(c_prime_box_scores)


#=== MISMATCH SCOREs ===
snodb_boxes$cbox_mismatches <- pmin(dist.prova1[snodb_boxes$c_seq], 
                                    dist.prova2[snodb_boxes$c_seq])
snodb_boxes$cbox_prime_mismatches <- pmin(dist.prova1[snodb_boxes$c_prime_seq], 
                                          dist.prova2[snodb_boxes$c_prime_seq])


#load [Cs]scores in snodb main df snodb_boxes
snodb_boxes$c_box_score <- cbox_scores[snodb_boxes$c_seq]
snodb_boxes$cprime_box_score <- c_prime_box_scores[snodb_boxes$c_prime_seq]


#==============================================================
#                D-BOX (and D') SCORE
#==============================================================
prova <-utils::combn(rep(c("C", "T", "G", "A"), 4),4)
prova.str = unique(apply(prova, 2, paste0, collapse=""))
dist.prova = sapply(prova.str, function(x) adist(x, "CTGA"))
possible_dboxes <- prova.str[dist.prova <= 3]


#=== PFM ===
pfm_d <- consensusMatrix(d, as.prob = T)
pfm_d_prime <- consensusMatrix(d_prime, as.prob = T)

#=== ACTUAL SCOREs ===
dbox_scores <- sapply(possible_dboxes, PWMscoreStartingAt, 
                      pwm = as.matrix(pfm_d))
dbox_scores_prop <- dbox_scores/max(dbox_scores)

d_prime_box_scores <- sapply(possible_dboxes, PWMscoreStartingAt, 
                             pwm = as.matrix(pfm_d_prime))
d_prime_box_scores_prop <- d_prime_box_scores/max(d_prime_box_scores)



#=== MISMATCH SCOREs ===
snodb_boxes$dbox_mismatches <- dist.prova[snodb_boxes$d_seq]
snodb_boxes$dbox_prime_mismatches <- dist.prova[snodb_boxes$d_prime_seq]

#load [Ds]scores in snodb main df
snodb_boxes$dprime_box_score <- d_prime_box_scores[snodb_boxes$d_prime_seq]
snodb_boxes$d_box_score <- dbox_scores[snodb_boxes$d_seq]



# trying both pois and nbin and return best fitting
fit_function <- function(vect){
  fit_pois <- fitdist(vect, distr= "pois", method = "mle")
  fit_nbin <- fitdist(vect, distr = "nbinom", method = "mle")
  dist <- seq(min(vect),max(vect), 1)
  fit_pois <- dpois(dist, lambda = fit_pois$estimate)
  fit_nbin <- dnbinom(dist,  size= fit_nbin$estimate[1], mu = fit_nbin$estimate[2])
  df_fit_pois <- data.frame(dist, fit_pois)
  df_fit_nbin <- data.frame(dist, fit_nbin)
  vect_df <- as.data.table(table(vect))
  vect_df[, N_norm := N/sum(N)]
  vect_df <- vect_df[order(vect)]
  df_fit_pois <- merge(df_fit_pois, vect_df, by.x = "dist", by.y = "vect", all.x=T)
  df_fit_pois$N_norm[is.na(df_fit_pois$N_norm)] <-0 
  df_fit_nbin <- merge(df_fit_nbin, vect_df, by.x = "dist", by.y = "vect", all.x=T)
  df_fit_nbin$N_norm[is.na(df_fit_nbin$N_norm)] <-0 
  pois_cor <- cor(df_fit_pois$fit_pois, df_fit_pois$N_norm)
  nbin_cor <- cor(df_fit_nbin$fit_nbin, df_fit_nbin$N_norm)
  print(pois_cor)
  print(nbin_cor)
  if(pois_cor>nbin_cor){
    best_fit <- (df_fit_pois)
  } else{
    best_fit <-(df_fit_nbin)
  }
  best_fit$score <- c(best_fit[,2])/max(c(best_fit[,2]))
  return(best_fit)
}
snodb_boxes$motif_score <- cbox_scores_prop[snodb_boxes$c_seq]+
  c_prime_box_scores_prop[snodb_boxes$c_prime_seq]*0.5+
  dbox_scores_prop[snodb_boxes$d_seq]+
  d_prime_box_scores_prop[snodb_boxes$d_prime_seq]*0.5
snodb_boxes$up_motif_score <- cbox_scores_prop[snodb_boxes$c_seq]+
  d_prime_box_scores_prop[snodb_boxes$d_prime_seq]*0.5
snodb_boxes$down_motif_score <- c_prime_box_scores_prop[snodb_boxes$c_prime_seq]*0.5+
  dbox_scores_prop[snodb_boxes$d_seq]
c_d_dist_fit <- fit_function(snodb_boxes$dist_c_d)
v_c_d_dist_fit <- setNames(c_d_dist_fit$score, c_d_dist_fit$dist)
c_prime_d_prime_dist_fit <- fit_function(snodb_boxes$dist_d_prime_c_prime)
v_c_prime_d_prime_dist_fit <- setNames(c_prime_d_prime_dist_fit$score, c_prime_d_prime_dist_fit$dist)
c_prime_d_dist_fit <- fit_function(snodb_boxes$dist_c_prime_d)
v_c_prime_d_dist_fit <- setNames(c_prime_d_dist_fit$score, c_prime_d_dist_fit$dist)
c_d_prime_dist_fit <- fit_function(snodb_boxes$dist_c_d_prime)
v_c_d_prime_dist_fit <- setNames(c_d_prime_dist_fit$score, c_d_prime_dist_fit$dist)
snodb_boxes$dist_c_d_prime_score <- v_c_d_prime_dist_fit[as.character(snodb_boxes$dist_c_d_prime)]
snodb_boxes$dist_c_prime_d_score <- v_c_prime_d_dist_fit[as.character(snodb_boxes$dist_c_prime_d)]
snodb_boxes$dist_d_prime_c_prime_score <-v_c_prime_d_prime_dist_fit[as.character(snodb_boxes$dist_d_prime_c_prime)]
snodb_boxes$dist_c_d_score <- v_c_d_dist_fit[as.character(snodb_boxes$dist_c_d)]

snodb_boxes$dist_score <- v_c_d_prime_dist_fit[as.character(snodb_boxes$dist_c_d_prime)]+
  v_c_prime_d_dist_fit[as.character(snodb_boxes$dist_c_prime_d)]+
  v_c_prime_d_prime_dist_fit[as.character(snodb_boxes$dist_d_prime_c_prime)]+
  v_c_d_dist_fit[as.character(snodb_boxes$dist_c_d)]


#=== EXPORTING THE PROCESSED DATA ===
save(possible_cboxes, possible_dboxes, file = get_output("possible_box"))
save(cbox_scores, cbox_scores_prop,
     c_prime_box_scores, c_prime_box_scores_prop,
     d_prime_box_scores, d_prime_box_scores_prop,
     dbox_scores, dbox_scores_prop,
     file = get_output("scores"))
# posso ripassarlo nello stesso punto tanto computing_guide e negatives hanno bisogno anche degli score quindi non possono essere eseguiti prima di computing_scores nella pipeline
save(met_sites, snodb_boxes, file = get_input("info_box")) 
save(c_d_prime_dist_fit, v_c_d_prime_dist_fit,
     c_d_dist_fit, v_c_d_dist_fit,
     c_prime_d_prime_dist_fit, v_c_prime_d_prime_dist_fit,
     c_prime_d_dist_fit, v_c_prime_d_dist_fit,
     file = get_output("dist_fit"))

#=== PLOTTINI ===#
if(get_config("generate_plots")){
  
#==============================================================
#               MOTIVES PLOTS
#==============================================================
  # all motifs <- why is it commented ? chris what u hidin?
  snodb_boxes$motif_score <- cbox_scores_prop[snodb_boxes$c_seq]+
            c_prime_box_scores_prop[snodb_boxes$c_prime_seq]*0.5+
                             dbox_scores_prop[snodb_boxes$d_seq]+
            d_prime_box_scores_prop[snodb_boxes$d_prime_seq]*0.5
  
  ggplot(snodb_boxes, aes(motif_score))+
    geom_histogram(fill = "#434384")+
    theme_bw(base_size = 20)+
    xlab("snoRNA Boxes score")+
    ylab("# of events")
  # scale_x_continuous(breaks= c(2.8,3, 3.2,3.4,3.6,3.8,4))
  ggsave(get_output("motif_score_distribution"), device = cairo_pdf, width = 6, height = 5)
  # minimum is 2.17
  
  # the couples
  snodb_boxes$up_motif_score <- cbox_scores_prop[snodb_boxes$c_seq]+
    d_prime_box_scores_prop[snodb_boxes$d_prime_seq]*0.5
  
  ggplot(snodb_boxes, aes(up_motif_score))+
    geom_histogram(fill = "#434384")+
    theme_bw(base_size = 20)+
    xlab("snoRNA Boxes score")+
    ylab("# of events")+
    scale_x_continuous(breaks= c(1, 1.2, 1.4,1.6,1.8,2))
  ggsave(get_output("up_motif_score_distribution"), device = cairo_pdf, width = 6, height = 5)
  # minimum is 0.98
  
  snodb_boxes$down_motif_score <- c_prime_box_scores_prop[snodb_boxes$c_prime_seq]*0.5+
    dbox_scores_prop[snodb_boxes$d_seq]
  
  ggplot(snodb_boxes, aes(down_motif_score))+
    geom_histogram(fill = "#434384")+
    theme_bw(base_size = 20)+
    xlab("snoRNA Boxes score")+
    ylab("# of events")+
    scale_x_continuous(breaks= c(1, 1.2, 1.4,1.6,1.8,2))
  ggsave(get_output("down_motif_score_distribution"), device = cairo_pdf, width = 6, height = 5)
  
#==============================================================
#               MISMATCH PLOT
#==============================================================
  plot_df <- snodb_boxes %>% pivot_longer(c("cbox_mismatches", "cbox_prime_mismatches", "dbox_mismatches", "dbox_prime_mismatches"))
  plot_df$name <- factor(plot_df$name, levels = c("dbox_mismatches", "dbox_prime_mismatches","cbox_mismatches", "cbox_prime_mismatches"), labels = c("D-Box", "D'-Box","C-Box","C'-Box"))
  
  plot_df$box <- "D-Box"
  plot_df$box[plot_df$name %in% c("C-Box", "C'-Box")] <- "C-Box"
  
  plot_df$box <- factor(plot_df$box, levels = c("D-Box", "C-Box"))
  
  
  ggplot(plot_df, aes(value, color = name,fill=name))+
    geom_bar(position = position_dodge2(preserve = "single"))+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    xlab("# of mismatches from canonical box")+
    ylab("# of events")+
    geom_text(stat="count", aes(x=value, label = stat(count)),position = position_dodge2(width = 0.9), vjust = -0.3, size = 4, show.legend = FALSE)+
    scale_fill_manual(values = c("#d4aa00", "#ead580","#2e2e5b", "#7e81a9"))+
    scale_color_manual(values = c("#947600", "#d7b323","#434384", "#52547a"))+
    facet_grid(cols = vars(box), space = "free", scales = "free")+
    theme(
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )+
    guides(fill=guide_legend(title="Box type"), color = guide_legend(title="Box type"))
  
  ggsave(get_output("box_mismatch_distribution"), device = cairo_pdf, width=7, height = 5)
  
#==============================================================
#               DISTANCES PLOTS
#==============================================================  
  # pipe operator is just doin g pivot(snodb, c( ...)) -> instead of passing it a parameter is passing through pipe
  plot_df <- snodb_boxes %>% pivot_longer(c("dist_d_prime_c_prime", 
                                            "dist_c_prime_d", 
                                            "dist_c_d_prime", 
                                            "dist_c_d"))
  plot_df$name <- factor(plot_df$name, 
                         levels = c("dist_d_prime_c_prime", 
                                    "dist_c_prime_d", 
                                    "dist_c_d_prime", 
                                    "dist_c_d"), 
                         labels = c("C-D'", "D'-C'", "C'-D", "C-D"))
  
  ggplot(plot_df, aes(name,value))+
    geom_violin(alpha = 0.7,  color = "#2e2e5b", fill = "#434384", scale = "width")+
    geom_boxplot(fill = "white", width= 0.2, outlier.shape = NA)+
    theme_bw(base_size = 20)+
    xlab("Boxes")+
    ylab("Distances (nt)")+
    # scale_y_log10(breaks = c(10, 20,30,40,50,1000))
    scale_y_log10(breaks =c(1,5,10,20,50, 100, 200, 400, 800))+
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank())
  
  ggsave(get_output("all_distance"), device = cairo_pdf(), width = 6, height = 5)
 
  # Create individual plots to be merged
  plot1 <- ggplot(snodb_boxes, aes(dist_d_prime_c_prime, after_stat(count)/sum(after_stat(count)))) +
    geom_bar(fill = "#434384") +
    theme_bw() +
    xlab("C' - D' dist") +
    ylab("Density") +
    scale_x_continuous(breaks = c(0,10,20,40,60, 80, 100, 125, 150)) +
    theme(panel.grid.minor = element_blank())
  
  plot2 <- ggplot(snodb_boxes, aes(dist_c_prime_d, after_stat(count)/sum(after_stat(count)))) +
    geom_bar(fill = "#434384") +
    theme_bw() +
    xlab("D' - C dist") +
    ylab("Density") +
    scale_x_continuous(breaks = c(0,10,20,40,60, 80, 100, 125, 150)) +
    theme(panel.grid.minor = element_blank())
  
  plot3 <- ggplot(snodb_boxes, aes(dist_c_d_prime, after_stat(count)/sum(after_stat(count)))) +
    geom_bar(fill = "#434384") +
    theme_bw() +
    xlab("D - C' dist") +
    ylab("Density") +
    scale_x_continuous(breaks = c(0,10,20,40,60, 80, 100, 120, 140)) +
    theme(panel.grid.minor = element_blank())
  
  plot4 <- ggplot(snodb_boxes, aes(dist_c_d, after_stat(count)/sum(after_stat(count)))) +
    geom_bar(fill = "#434384") +
    theme_bw() +
    xlab("C - D dist") +
    ylab("Density") +
    scale_x_continuous(breaks = c(0,10,20,40,60, 80, 100, 120, 140, 160, 180, 200, 220)) +
    theme(panel.grid.minor = element_blank())
  
  # Combine in 2x2 layout
  combined_plot <- (plot1 + plot2) / (plot3 + plot4)
  
  # Add overall title (optional)
  combined_plot <- combined_plot + 
    plot_annotation(
      title = "snoRNA Box Distance Distributions",
      theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
    )
  
  # Save as single file -> ho visto che doposalva altre cose come distribution quindi devo trovare un altro nome
  ggsave(get_output("all_distance_distributions"), 
         combined_plot, 
         device = cairo_pdf, 
         width = 10, height = 8)

    ggplot(c_d_prime_dist_fit, aes(dist,N_norm))+
    geom_col(fill = "#434384")+
    geom_line(aes(dist, fit_nbin), color = "#d4aa00", linewidth = 1)+
    theme_bw()+
    xlab("Dist C-D' Box")+
    ylab("Density")+
    scale_color_manual(values = c(snoRNAs = "#434384", Nbinom = "#d4aa00"))
  ggsave("snoDB_analysis_validated/c_d_prime_dist_nbin_dist.pdf", device = cairo_pdf, width = 4, height = 3)
  
  ggplot(c_prime_d_dist_fit, aes(dist,N_norm))+
    geom_col(fill = "#434384")+
    geom_line(aes(dist, fit_nbin), color = "#d4aa00", linewidth = 1)+
    theme_bw()+
    xlab("Dist C'-D Box")+
    ylab("Density")+
    scale_color_manual(values = c(snoRNAs = "#434384", Nbinom = "#d4aa00"))
  ggsave("snoDB_analysis_validated/c_prime_d_dist_nbin_dist.pdf", device = cairo_pdf, width = 4, height = 3)
  
  ggplot(c_prime_d_prime_dist_fit, aes(dist,N_norm))+
    geom_col(fill = "#434384")+
    geom_line(aes(dist, fit_nbin), color = "#d4aa00", linewidth = 1)+
    theme_bw()+
    xlab("Dist C'-D' Box")+
    ylab("Density")+
    scale_color_manual(values = c(snoRNAs = "#434384", Nbinom = "#d4aa00"))
  ggsave("snoDB_analysis_validated/c_prime_d_prime_dist_fit_nbin_dist.pdf", device = cairo_pdf, width = 4, height = 3)
  
  ggplot(c_d_dist_fit, aes(dist,N_norm))+
    geom_col(fill = "#434384")+
    geom_line(aes(dist, fit_pois), color = "#d4aa00", linewidth = 1)+
    theme_bw()+
    xlab("Dist C-D Box")+
    ylab("Density")+
    scale_color_manual(values = c(snoRNAs = "#434384", Nbinom = "#d4aa00"))
  ggsave("snoDB_analysis_validated/c_d_dist_pois_dist.pdf", device = cairo_pdf, width = 4, height = 3)
  
  ggplot(snodb_boxes, aes(dist_score, after_stat(count/sum(count))))+
    geom_histogram(fill = "#434384")+
    theme_bw()+
    scale_x_continuous(breaks = c(0.5,1,1.5,2,2.5,3,3.5,4))+
    ylab("Density")
  ggsave("snoDB_analysis_validated/dist_score_dist.pdf", device = cairo_pdf, width = 4, height = 3)
  # min is 0.69
  
}

#=== EXPORTING PFM TABLES ===#
if(get_config("export_tables")){
  
  write.csv(pfm_c, get_output("pfm_cbox_snoDb"))
  write.csv(pfm_c_prime, get_output("pfm_cbox_prime_snoDb"))
  write.csv(pfm_d, get_output("pfm_dbox_snoDb"))
  write.csv(pfm_d_prime, get_output("pfm_dbox_prime_snoDb"))
  
  write_xlsx(snodb_boxes, "snoDB_analysis_validated/snodb_analysis.xlsx")
}
