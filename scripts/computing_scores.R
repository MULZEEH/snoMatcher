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

# Check if running in Snakemake or RStudio
if (!exists("snakemake")) {
  #TODO!
}

load(file = snakemake@input[["info_boxes"]], verbose = TRUE)

c <- c(snodb_boxes$c_seq[str_count(snodb_boxes$c_seq) == 7])
c_prime <- c(snodb_boxes$c_prime_seq)
d <- c(snodb_boxes$d_seq)
d_prime <- c(snodb_boxes$d_prime_seq)

# FOR EACH OF THE BOX I COMPUTE THE SCORES (can i loop it?)
prova <-utils::combn(rep(c("C", "T", "G", "A"), 7),7)
prova.str = unique(apply(prova, 2, paste0, collapse=""))
dist.prova1 = sapply(prova.str, function(x) adist(x, "ATGATGA"))
dist.prova2 = sapply(prova.str, function(x) adist(x, "GTGATGA"))


possible_cboxes1 <- prova.str[dist.prova1 <= 4]
possible_cboxes2 <- prova.str[dist.prova2 <= 4]

possible_cboxes <- union(possible_cboxes1, possible_cboxes2)

# save(possible_cboxes, file = "snoDB_analysis_validated/possible_cboxes.RData")

pfm_c <- consensusMatrix(c, as.prob = T)

# if statement for the print makes much more sense
write.csv(pfm_c, "snoDB_analysis_validated/pfm_cbox_snoDb_snoRNA.csv")

# c[dist.prova1[c]>3]
# c[dist.prova2[c]>3]

cbox_scores <- sapply(possible_cboxes, PWMscoreStartingAt, pwm = as.matrix(pfm_c))
cbox_scores_prop <- cbox_scores/max(cbox_scores)
# save(cbox_scores, cbox_scores_prop, file = "snoDB_analysis_validated/cbox_scores.RData")

pfm_c_prime <- consensusMatrix(c_prime, as.prob = T)
write.csv(pfm_c_prime, "snoDB_analysis_validated/pfm_c_prime_box_snoDb_snoRNA.csv")

# c_prime[dist.prova1[c_prime]>3]
# c_prime[dist.prova2[c_prime]>3]

c_prime_box_scores <- sapply(possible_cboxes, PWMscoreStartingAt, pwm = as.matrix(pfm_c_prime))
c_prime_box_scores_prop <- c_prime_box_scores/max(c_prime_box_scores)
# save(c_prime_box_scores, c_prime_box_scores_prop, file = "snoDB_analysis_validated/c_prime_box_scores.RData")

prova <-utils::combn(rep(c("C", "T", "G", "A"), 4),4)
prova.str = unique(apply(prova, 2, paste0, collapse=""))
dist.prova = sapply(prova.str, function(x) adist(x, "CTGA"))
possible_dboxes <- prova.str[dist.prova <= 3]
save(possible_dboxes, file = "snoDB_analysis_validated/possible_dboxes.RData")

pfm_d <- consensusMatrix(d, as.prob = T)
write.csv(pfm_d, "snoDB_analysis_validated/pfm_Dbox_snoDb_snoRNA.csv")

# d[dist.prova[d]>2]

dbox_scores <- sapply(possible_dboxes, PWMscoreStartingAt, pwm = as.matrix(pfm_d))
dbox_scores_prop <- dbox_scores/max(dbox_scores)
save(dbox_scores, dbox_scores_prop, file = "snoDB_analysis_validated/dbox_scores.RData")

pfm_d_prime <- consensusMatrix(d_prime, as.prob = T)
write.csv(pfm_d_prime, "snoDB_analysis_validated/pfm_D_prime_box_snoDb_snoRNA.csv")

# d_prime[dist.prova[d_prime]>2]

d_prime_box_scores <- sapply(possible_dboxes, PWMscoreStartingAt, pwm = as.matrix(pfm_d_prime))
d_prime_box_scores_prop <- d_prime_box_scores/max(d_prime_box_scores)
save(d_prime_box_scores, d_prime_box_scores_prop, file = "snoDB_analysis_validated/d_prime_box_scores.RData")


snodb_boxes$c_box_score <- cbox_scores[snodb_boxes$c_seq]
snodb_boxes$dprime_box_score <- d_prime_box_scores[snodb_boxes$d_prime_seq]
snodb_boxes$cprime_box_score <- c_prime_box_scores[snodb_boxes$c_prime_seq]
snodb_boxes$d_box_score <- dbox_scores[snodb_boxes$d_seq]


# all motifs ->> this needs to get executed for the sequents ggplot
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
ggsave("snoDB_analysis_validated/motif_score_snoDB_distribution.pdf", device = cairo_pdf, width = 6, height = 5)
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
ggsave("snoDB_analysis_validated/up_motif_score_snoDB_distribution.pdf", device = cairo_pdf, width = 6, height = 5)
# minimum is 0.98

snodb_boxes$down_motif_score <- c_prime_box_scores_prop[snodb_boxes$c_prime_seq]*0.5+
  dbox_scores_prop[snodb_boxes$d_seq]

ggplot(snodb_boxes, aes(down_motif_score))+
  geom_histogram(fill = "#434384")+
  theme_bw(base_size = 20)+
  xlab("snoRNA Boxes score")+
  ylab("# of events")+
  scale_x_continuous(breaks= c(1, 1.2, 1.4,1.6,1.8,2))
ggsave("snoDB_analysis_validated/down_motif_score_snoDB_distribution.pdf", device = cairo_pdf, width = 6, height = 5)
# minimum is 0.76

prova <-utils::combn(rep(c("C", "T", "G", "A"), 7),7)
prova.str = unique(apply(prova, 2, paste0, collapse=""))
dist.prova1 = sapply(prova.str, function(x) adist(x, "ATGATGA"))
dist.prova2 = sapply(prova.str, function(x) adist(x, "GTGATGA"))
snodb_boxes$cbox_mismatches <- pmin(dist.prova1[snodb_boxes$c_seq], dist.prova2[snodb_boxes$c_seq])
snodb_boxes$cbox_prime_mismatches <- pmin(dist.prova1[snodb_boxes$c_prime_seq], dist.prova2[snodb_boxes$c_prime_seq])

prova <-utils::combn(rep(c("C", "T", "G", "A"), 4),4)
prova.str = unique(apply(prova, 2, paste0, collapse=""))
dist.prova = sapply(prova.str, function(x) adist(x, "CTGA"))
snodb_boxes$dbox_mismatches <- dist.prova[snodb_boxes$d_seq]
snodb_boxes$dbox_prime_mismatches <- dist.prova[snodb_boxes$d_prime_seq]

plot_df <- snodb_boxes %>% pivot_longer(c("cbox_mismatches", "cbox_prime_mismatches", "dbox_mismatches", "dbox_prime_mismatches"))
plot_df$name <- factor(plot_df$name, levels = c("dbox_mismatches", "dbox_prime_mismatches","cbox_mismatches", "cbox_prime_mismatches"), labels = c("D-Box", "D'-Box","C-Box", "C'-Box"))
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

ggsave("snoDB_analysis_validated/box_mismatch_distribution.pdf", device = cairo_pdf, width=7, height = 5)

# removing snord10, 31, 111 - 4 mismatches in c'-box
# removing snord19C - 3 mismatches in d'-box

snodb_boxes$dist_d_prime_c_prime <- snodb_boxes$c_prime_start - (snodb_boxes$d_prime_start+3)
# 1-51
# removing SNORD17


snodb_boxes$dist_c_d_prime <- snodb_boxes$d_prime_start - 
  (snodb_boxes$c_start+6)
# 6-50

snodb_boxes$dist_c_prime_d <- snodb_boxes$d_start-
  (snodb_boxes$c_prime_start+6)
# 8-63

# total distance C - D
snodb_boxes$dist_c_d <- snodb_boxes$d_start - (snodb_boxes$c_start+6)
# 40-126
# removing SNORD17

# snodb_boxes$dist_d_prime_d <- snodb_boxes$d_start-(snodb_boxes$d_prime_start+3)
# # 22-109
# # removing SNORD17
# 
# snodb_boxes$dist_c_c_prime <- snodb_boxes$c_prime_start-(snodb_boxes$c_start+6)
# # 17-98

# median(snodb_boxes$dist_d_prime_c)
# # 17
# 
# median(snodb_boxes$dist_c_prime_d_prime)
# # 8
# 
# median(snodb_boxes$dist_d_c_prime)
# # 17
# 
# median(snodb_boxes$dist_c_d)
# # 56

plot_df <- snodb_boxes %>% pivot_longer(c("dist_c_prime_d_prime", "dist_d_prime_c", "dist_d_c_prime", "dist_c_d"))
plot_df$name <- factor(plot_df$name, levels = c("dist_d_prime_c", "dist_c_prime_d_prime", "dist_d_c_prime", "dist_c_d"), labels = c("C-D'", "D'-C'", "C'-D", "C-D"))

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

ggsave("snoDB_analysis_validated/all_distances.pdf", device = cairo_pdf(), width = 6, height = 5)


ggplot(snodb_boxes, aes(dist_c_prime_d_prime, after_stat(count)/sum(after_stat(count))))+
  geom_bar(fill = "#434384")+
  theme_bw()+
  xlab("C' - D' dist")+
  ylab("Density")+
  scale_x_continuous(breaks =c(0,10,20,40,60, 80, 100, 125, 150))+
  theme(panel.grid.minor = element_blank())
ggsave("snoDB_analysis_validated/c_prime_d_prime_dist.pdf", device = cairo_pdf, width = 4, height = 3)

ggplot(snodb_boxes, aes(dist_d_prime_c, after_stat(count)/sum(after_stat(count))))+
  geom_bar(fill = "#434384")+
  theme_bw()+
  xlab("C - D' dist")+
  ylab("Density")+
  scale_x_continuous(breaks =c(0,10,20,40,60, 80, 100, 125, 150))+
  theme(panel.grid.minor = element_blank())
ggsave("snoDB_analysis_validated/d_prime_c_dist.pdf", device = cairo_pdf, width = 4, height = 3)

ggplot(snodb_boxes, aes(dist_d_c_prime, after_stat(count)/sum(after_stat(count))))+
  geom_bar(fill = "#434384")+
  theme_bw()+
  xlab("C' - D dist")+
  ylab("Density")+
  scale_x_continuous(breaks =c(0,10,20,40,60, 80, 100, 120, 140))+
  theme(panel.grid.minor = element_blank())
ggsave("snoDB_analysis_validated/c_prime_d_dist.pdf", device = cairo_pdf, width = 4, height = 3)

ggplot(snodb_boxes, aes(dist_c_d, after_stat(count)/sum(after_stat(count))))+
  geom_bar(fill = "#434384")+
  theme_bw()+
  xlab("C - D dist")+
  ylab("Density")+
  scale_x_continuous(breaks =c(0,10,20,40,60, 80, 100, 120, 140, 160, 180, 200, 220))+
  theme(panel.grid.minor = element_blank())
ggsave("snoDB_analysis_validated/c_d_dist.pdf", device = cairo_pdf, width = 4, height = 3)