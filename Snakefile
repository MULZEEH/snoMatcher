
# snoRNA Feature Analysis Pipeline

# snoMatcher/
# |-Snakefile
# |-config.yaml
# |-scripts/
# | |-somescripts
# |-results/
# | |-plots/
# | |-tables/
# | |-intermediate or processed data/?
# | |-final ?
# |-data/
# | |-raw/
# | |-.../
# |-env/
# | |-sssss.yaml or something
#


#import pandas as pd
#samples = pd.read_csv("sample.tsv", sep='\t').set_index("sample", drop= False)

# POSSIBLE IMPROVEMENT:|- add some additional configuration for the plottingsystem

configfile: "config.yaml"

# Output is defined based on the flag chosen
def get_final_output():
  outs = [
    # DATA THAT WILL ALWAYS BE PRESENT
    "results/noidea.somedataRdata" 
    # Maybe to be added in the final folder (results/final/...)
  ]
  
  
  if config.get("generate_plots", True):
    outs.extend([
      #actually are tons of these
      "results/plots/someplots.pdf"
    ])
  
  if config.get("export_tables", True):
    outs.extend([
      "results/tables/sometables.csv"
    ])
    
  return outs


# Final target -> what we want to have in the end (i guess some validation report and Rdata?)
rule all:
    input:
      get_final_output()

###=====================================================
###=====================================================
###========= CORE DATA PROCESSING  =====================
###=====================================================
###=====================================================

#=====================================================
# SETUP - Raw Data Loading and getting the firsts PLOTS
#
# POSSIBLE IMPROVEMENT:|- maybe "parallelizing" the computational of c/d/c'/d' since they are independent
#                       |-> adding the cores/threads filed
#=====================================================
rule setup : 
	 input: 
	   # dir_setup = "results/.dir_setup_complete",
	   # sno_db = config["input_files"]["snodb"] -> in yaml: snodb: "path"
	   sno_db = "data/raw/snoDB_data.xlsx",
	   met_sites = "data/raw/snoDB_rRNA_interactions.xlsx",
	   sno_boxes = "data/raw/cd_boxes.tsv"
	   
	 output: 
	   info_box = "results/intermediate/info_box.RData",
	   #if i add this i cannot run snakemake wih --config generate_plots = false (nespole)
	   all_logos = "results/plots/all_logos.pdf",
	   guide_width = "results/plots/guide_width.pdf"
	   
	   
	 script: 
	   "scripts/setup.R"

#=====================================================
# SCORING SYSTEM
#=====================================================
rule computing_scores :
	 input:
	   info_box = "results/intermediate/info_box.Rdata"
	   
	 output: 
# 	   possible_box = "results/intermediate/possible_boxes.RData",
#      scores = "results/intermediate/scores.RData",
#      box_mismatch_distribution = "results/plots/box_mismatch_distribution.pdf",
#       #violin
#      all_distance = "results/plots/all_distance.pdf",
#       #"histogram" density
#      relative_dist = "results/plots/relative_dist.pdf",
#      motif_score_distribution = "results/plots/motif_score_snoDB_distribution.pdf",
#      up_motif_score_distribution = "results/plots/up_motif_score_distribution.pdf",
#      down_motif_score_distribution = "results/plots/down_motif_score_distribution.pdf",
#       
#       #TABLES
#       
#      pfm_c_box_snoDb_snoRNA = "results/tables/pfm_c_box_snoDb_snoRNA.csv",
#      pfm_c_prime_box_snoDb_snoRNA = "results/tables/pfm_c_prime_box_snoDb_snoRNA.csv",
#      pfm_d_box_snoDb_snoRNA = "results/tables/pfm_d_box_snoDb_snoRNA.csv",
#      pfm_d_prime_box_snoDb_snoRNA = "results/tables/pfm_d_prime_box_snoDb_snoRNA.csv"
	   
	 script: 
	   "scripts/computing_scores.R"

#=====================================================
# PAIRING DISTRIBUTIONS AND GUIDE -> TODO ctrl+shift+c is ma friend
#=====================================================
 
# rule pairing_distributions :
# 	 input:
# 	   snodb_boxes = ...
# 	   metsites = (done in the setup.R)
# 	   human_genome = .fasta
# 	 output: 
# 	 script: 

#=====================================================
#
#=====================================================
 
# rule other :
# 	 input: 
# 	 output: 
# 	 script: 
	   

#=====================================================
# UTILITY RULES
#=====================================================

rule clean:
  output:
    ".dir_clean_ready"
  shell:
    """
        rm -rf results/
        touch {output}
        """

# rule clean_plot:
#   shell:
#     "rm -rf results."


#not working
rule setup_directory:
  input:
    clean_setup = ".dir_clean_ready"
  output: 
    dir_setup = "results/.dir_setup_complete"
  # message: "setup Directories"
  shell:
        """
        mkdir  results
        mkdir  results/intermediate
        mkdir  results/plots
        mkdir  results/tables
        touch {output}
        """



# 
# ###=====================================================
# ###================TO DELETE============================
# ###============= OPTIONAL DATA  ========================
# ###============ (tables & plots) =======================
# ###=====================================================
# # IT MAKEES NO SENSE TO DO SOMETHING LIKE THIS, JUST LOOK AT A snakemake@config var that tells me if i have to do a condition or not
# if config.get("generate_plots", True):
#   rule setup_plots:
#     input:
#       info_boxes = "results/intermediate/info_boxes.Rdata"
#       
#     output:
#       all_logos = "results/plots/all_logos.pdf",
#       guide_width = "results/plots/guide_width.pdf"
#       
#     script:
#       "scripts/setup_plotting.R"
#       
#   rule scores_plots:
#     input:
#       info_boxes = "results/intermediate/info_boxes.Rdata"
#       # 
#     output:
#        box_mismatch_distribution = "results/plots/box_mismatch_distribution.pdf",
#   	   all_distance = "results/plots/all_distance.pdf",
#   	   c_prime_d_prime_dist = "results/plots/c_prime_d_prime_dist.pdf",
#   	   #...and such 
#   	   motif_score_distribution = "results/plots/motif_score_snoDB_distribution.pdf",
#   	   up_motif_score_distribution = "results/plots/up_motif_score_distribution.pdf",
#   	   down_motif_score_distribution = "results/plots/down_motif_score_distribution.pdf"
#   	   
#   	 script:
#   	   "scripts/scores_plotting.R"
#   	   
#   if config.get("export_tables", True):
#     rule scores_tables:
#       input: 
#         info_boxes = "results/intermediate/info_boxes.Rdata"
#         
#       output:
#         pfm_c_box_snoDb_snoRNA = "results/tables/pfm_c_box_snoDb_snoRNA.csv",
#         pfm_c_prime_box_snoDb_snoRNA = "results/tables/pfm_c_prime_box_snoDb_snoRNA.csv",
#         pfm_d_box_snoDb_snoRNA = "results/tables/pfm_d_box_snoDb_snoRNA.csv",
#         pfm_d_prime_box_snoDb_snoRNA = "results/tables/pfm_d_prime_box_snoDb_snoRNA.csv"
#         
#       script:
#         "scripts/score_tables.R"
#         
#   
