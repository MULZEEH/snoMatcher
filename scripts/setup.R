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

#
# Pulire il Global Environment prima di eseguire l'if condition
# Al momento controlla se esiste snakemake (come avviene nella
# normale procedura con snakemake) in fase di debug pero' non 
# deve esistere quindi se si e' eseguito l'if una volta la volta
# dopo questo non considera piu di essere in Debug mode
#
#
# Per eseguire lo script mediante snakemake 
# |-> snakemake --cores N results/intermediate/info_box.Rdata
#



# Check if running with Snakemake or in RStudio
# Need to clear the environment first -> TODO(FIX)
rm(snakemake, envir = .GlobalEnv)
if (!exists("snakemake")) {
  setwd("C:/Users/Marco/RProject/snoMatcher/")
  # Create mock snakemake object for testing in Rstudio
  snakemake <- list(
    input = list(
      sno_db = "data/raw/snoDB_data.xlsx",
      met_sites = "data/raw/snoDB_rRNA_interactions.xlsx",
      sno_boxes = "data/raw/cd_boxes.tsv"
    ),
    output = list(
      info_box = "results/intermediate/info_box.RData",
      all_logos = "results/plots/all_logos.pdf",
      guide_width = "results/plots/guide_width.pdf"
      ),
    config = list(
      generate_plots = TRUE,
      export_tables = TRUE
    )
  )
  
# ORRENDO -> TROVA SOLUZIONE NICER
  
  # Helper function for mock object
  get_input <- function(name) snakemake$input[[name]]
  get_output <- function(name) snakemake$output[[name]]
  get_config <- function(name) snakemake$config[[name]]
  # get_threads <- function() snakemake$threads
  cat("Debug mode...")
  
} else {
  # Helper functions for real snakemake object
  get_input <- function(name) snakemake@input[[name]]
  get_output <- function(name) snakemake@output[[name]]
  get_config <- function(name) snakemake@config[[name]]
  # get_threads <- function() snakemake@threads
  cat("snakemake execution...")
}

snodb <- read_xlsx(get_input("sno_db"))
met_sites <- read_xlsx(get_input("met_sites"))
snodb_boxes <- fread(get_input("sno_boxes"))

snodb_boxes <- snodb_boxes[snodb_boxes$c_start != -1 & snodb_boxes$c_prime_start != -1 & snodb_boxes$d_start != -1 & snodb_boxes$d_prime_start != -1,]
snodb_boxes <- merge(snodb[, c("snoDB ID", "Symbol", "DNA Sequence")], snodb_boxes, by.x = "snoDB ID", by.y = "snoDB_id")
snodb_boxes <- snodb_boxes[snodb_boxes$`snoDB ID` %in% met_sites$snoDB_id[met_sites$Status == "validated"],]

snoRNA_seq <- DNAStringSet(snodb_boxes$`DNA Sequence`)
names(snoRNA_seq) <- snodb_boxes$`snoDB ID`

met_sites <- met_sites[met_sites$snoDB_id %in% snodb_boxes$`snoDB ID` & met_sites$Status == "validated",]

snodb_boxes$guide1_width <- str_count(snodb_boxes$guide1_seq)
snodb_boxes$guide2_width <- str_count(snodb_boxes$guide2_seq)

# (c <- c(snodb_boxes$c_seq[str_count(snodb_boxes$c_seq) == 7]))
# (c_prime <- c(snodb_boxes$c_prime_seq))
# (d <- c(snodb_boxes$d_seq))
# (d_prime <- c(snodb_boxes$d_prime_seq))


# GENERATE PLOTS -> Generate a stranage WARNING 
if (get_config("generate_plots")) {
  
#---GUIDE WIDTH PLOT---#
  plot_df <- snodb_boxes %>% pivot_longer(c(guide1_width, guide2_width))
  
  ggplot(plot_df[plot_df$value != 0,], aes(value-1))+
    geom_bar(fill="#434384")+
    theme_bw()+
    geom_text(stat="count", aes(x=value-1, label = stat(count)),position = position_dodge(width = 0.9), vjust = -0.3, size = 3)+
    xlab("Guide width")+
    ylab("# of events")+
    scale_x_continuous(breaks = seq(7,18,1))+
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
  ggsave(get_output("guide_width"), device = cairo_pdf, width=4, height = 3)
  
  
  # Define box types and sequences
  box_data <- list(
    "C box" = snodb_boxes$c_seq[str_count(snodb_boxes$c_seq) == 7],
    "C' box" = snodb_boxes$c_prime_seq,
    "D box" = snodb_boxes$d_seq,
    "D' box" = snodb_boxes$d_prime_seq
  )
  
  # Create plots -> TOTEST!()
  plots <- map(names(box_data), function(box_name) {
    sequences <- box_data[[box_name]]
    ggseqlogo(gsub(x = sequences, pattern = "T", replacement = "U")) + 
      ggtitle(box_name) +
      theme(plot.title = element_text(hjust = 0.5, size = 14))
  })
  
  # Combine in 2x2 grid
  combined_plot <- wrap_plots(plots, ncol = 2, nrow = 2)
  
  ggsave(get_output("all_logos"), combined_plot, 
         device = cairo_pdf, width = 14, height = 10)
  
}

# EXPORT THE SNORNADB + THE METILATION SITES -> info_boxes.pdf
save(met_sites, snodb_boxes, file = get_output("info_box"))
