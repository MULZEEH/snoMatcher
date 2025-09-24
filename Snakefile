
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

# TODO 24/09:
# @use the config as actual config also for the Snakefile
#   LIKE
#   rule setup:
#     """Load raw data and perform initial processing"""
#     input:
#         sno_db = config["input_files"]["snodb"],
#         met_sites = config["input_files"]["met_sites"], 
#         sno_boxes = config["input_files"]["sno_boxes"],
#         dirs = rules.setup_directories.output
#     output:
#         info_box = "results/intermediate/info_box.RData",
#         all_logos = "results/plots/all_logos.pdf" if config.get("generate_plots", True) else [],
#         guide_width = "results/plots/guide_width.pdf" if config.get("generate_plots", True) else []
#     threads: config.get("threads", 1)
#     resources:
#         mem_mb = config.get("memory", 4000)
#     script:
#         "scripts/setup.R"
# @use expand and if condition inside the rule to expand or not the output based on the flags
# @set threads but implement it later mdear
# @take a look at the NEW configuration in config.yaml (in the end commented)



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
	   # sno_db = config["input_files"]["snodb"] -> in yaml: snodb: "path" -> will always switch to this
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
	   # NOT WORKING PROPERLY something like this should work:
	   # plots = expand("results/plots/{plot}.pdf", 
                      # plot=["box_mismatch_distribution", "all_distance", "motif_score_distribution"]) if config.get("generate_plots", True) else [],
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

rule model_comparing:
  input: 
    rf_model = "models/rf_model.RData",
    svm_model = "models/svm_model.RData",
    xgboost_model = "models/xgboost_model.RData",
    cnn_model = "models/cnn_model.RData"
  output:
    prediction_comparison = "results/final/models_comp_plot.pdf"
  
	   

#=====================================================
# UTILITY RULES
#=====================================================

rule clean:
    """Remove all results and temporary files"""
    shell:
        "rm -rf results/ .snakemake/"

rule clean_plots:
    """Remove only plot files"""
    shell:
        "rm -rf results/plots/"

rule clean_tables:
    """Remove only table files"""  
    shell:
        "rm -rf results/tables/"



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
