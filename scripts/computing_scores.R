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
      info_box = "results/intermediate/info_box.RData"
    ),
    output = list(

      possible_box = "results/intermediate/possible_boxes.RData",
      scores = "results/intermediate/scores.RData",
      #trovare un altro nome o solo sovrascrivere quello precedente
      processed_info_box = "results/intermediate/processed_info_box.RData",
      
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

load(file = get_input("info_box"))

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
# cbox_scores_prop <- cbox_scores/max(cbox_scores)

c_prime_box_scores <- sapply(possible_cboxes, PWMscoreStartingAt, pwm = as.matrix(pfm_c_prime))
# c_prime_box_scores_prop <- c_prime_box_scores/max(c_prime_box_scores)


#=== MISMATCH SCOREs ===
snodb_boxes$cbox_mismatches <- pmin(dist.prova1[snodb_boxes$c_seq], 
                                    dist.prova2[snodb_boxes$c_seq])
snodb_boxes$cbox_prime_mismatches <- pmin(dist.prova1[snodb_boxes$c_prime_seq], 
                                          dist.prova2[snodb_boxes$c_prime_seq])


#load [Cs]scores in snodb main df snodb_boxes
snodb_boxes$c_box_score <- cbox_scores[snodb_boxes$c_seq]
snodb_boxes$cprime_box_score <- c_prime_box_scores[snodb_boxes$c_prime_seq]

# CHECK
# c[dist.prova1[c]>3]
# c[dist.prova2[c]>3]

# CHECK
# c_prime[dist.prova1[c_prime]>3]
# c_prime[dist.prova2[c_prime]>3]


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
# dbox_scores_prop <- dbox_scores/max(dbox_scores)

d_prime_box_scores <- sapply(possible_dboxes, PWMscoreStartingAt, 
                             pwm = as.matrix(pfm_d_prime))
# d_prime_box_scores_prop <- d_prime_box_scores/max(d_prime_box_scores)



#=== MISMATCH SCOREs ===
snodb_boxes$dbox_mismatches <- dist.prova[snodb_boxes$d_seq]
snodb_boxes$dbox_prime_mismatches <- dist.prova[snodb_boxes$d_prime_seq]

#load [Ds]scores in snodb main df
snodb_boxes$dprime_box_score <- d_prime_box_scores[snodb_boxes$d_prime_seq]
snodb_boxes$d_box_score <- dbox_scores[snodb_boxes$d_seq]

# CHECK
# d[dist.prova[d]>2]
# d_prime[dist.prova[d_prime]>2]


#=== EXPORTING THE PROCESSED DATA ===
save(possible_cboxes, possible_dboxes, file = get_output("possible_box"))
save(cbox_scores, cbox_scores_prop,
     c_prime_box_scores, c_prime_box_scores_prop,
     d_prime_box_scores, d_prime_box_scores_prop,
     dbox_scores, dbox_scores_prop,
     file = get_output("box_scores"))

# DO I REALLY NEED TO SAVE THE PROP? 
# save(cbox_scores, c_prime_box_scores,
#      dbox_scores, d_prime_box_scores, 
#      file = get_output("box_scores"))

#resaving the input (since i added something inthe snodb_boxes)
save(met_sites, snodb_boxes, file = get_output("processed_info_box"))

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
  
  # # Reshape data for easier plotting
  # plot_data <- snodb_boxes %>%
  #   select(dist_c_prime_d_prime, dist_d_prime_c, dist_d_c_prime, dist_c_d) %>%
  #   pivot_longer(everything(), names_to = "distance_type", values_to = "distance") %>%
  #   mutate(
  #     distance_type = case_when(
  #       distance_type == "dist_c_prime_d_prime" ~ "C' - D' dist",
  #       distance_type == "dist_d_prime_c" ~ "D' - C dist", 
  #       distance_type == "dist_d_c_prime" ~ "D - C' dist",
  #       distance_type == "dist_c_d" ~ "C - D dist"
  #     ),
  #     distance_type = factor(distance_type, 
  #                            levels = c("C' - D' dist", "D' - C dist", 
  #                                       "D - C' dist", "C - D dist"))
  #   )
  # 
  # # Create single plot with facets
  # combined_plot <- ggplot(plot_data, aes(distance, after_stat(count)/sum(after_stat(count)))) +
  #   geom_bar(fill = "#434384") +
  #   theme_bw() +
  #   facet_wrap(~distance_type, scales = "free_x", ncol = 2) +
  #   xlab("Distance") +
  #   ylab("Density") +
  #   theme(
  #     panel.grid.minor = element_blank(),
  #     strip.text = element_text(size = 12),
  #     plot.title = element_text(size = 16, hjust = 0.5)
  #   ) +
  #   ggtitle("snoRNA Box Distance Distributions")
  # 
  # # Save
  # ggsave("snoDB_analysis_validated/all_distance_distributions.pdf", 
  #        combined_plot, 
  #        device = cairo_pdf, 
  #        width = 10, height = 8)
}

#=== EXPORTING PFM TABLES ===#
if(get_config("export_tables")){
  
  write.csv(pfm_c, get_output("pfm_cbox_snoDb"))
  write.csv(pfm_c_prime, get_output("pfm_cbox_prime_snoDb"))
  write.csv(pfm_d, get_output("pfm_dbox_snoDb"))
  write.csv(pfm_d_prime, get_output("pfm_dbox_prime_snoDb"))
}
