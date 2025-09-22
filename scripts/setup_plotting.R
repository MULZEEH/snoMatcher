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
library(patchwork)

# Check if running in Snakemake or RStudio
if (!exists("snakemake")) {
  #TODO!
}

load(file = snakemake@input[["info_boxes"]], verbose = TRUE)

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
ggsave(snakemake@output[["guide_width"]], device = cairo_pdf, width=4, height = 3)


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

ggsave(snakemake@output[["all_logos"]], combined_plot, 
       device = cairo_pdf, width = 14, height = 10)