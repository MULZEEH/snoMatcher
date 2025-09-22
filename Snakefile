# snoRNA Feature Analysis Pipeline
configfile: "config.yaml"

# Final target -> what we want to have in the end (i guess some validation report and Rdata?)
rule all:
    input:
      "results/noidea.Rdata"

#=====================================================
# SETUP - Raw Data Loading and getting the firsts PLOTS
#=====================================================
rule setup :
	 input: 
	   sno_db = "data/raw/snoDB_data.xlsx"
	   met_sites = "data/raw/snoDB_rRNA_interactions.xlsx"
	   # ok different tsv
	   sno_boxes = "data/raw/cd_boxes.tsv"
	 output: 
	   guide_width_plot = "results/plots/guide_width.pdf"
	   # I think it would be better to have only 1 box_logo.pdf that incorporates all the other 4 <-
	   Cbox_logo = "results/plots/Cbox_logo.pdf"
	   Dbox_logo = "results/plots/Dbox_logo.pdf"
	   C_prime_box_logo = "results/plots/C_prime_logo.pdf"
	   D_prime_box_logo = "results/plots/D_prime_logo.pdf"
	 script: 
	   "scripts/setup.R"

#=====================================================
#
#=====================================================

 
rule matching :
	 input: 
	 output: 
	 script: 

#=====================================================
#
#=====================================================

 
rule boxes_scores :
	 input: 
	 output: 
	 script: 

#=====================================================
#
#=====================================================

 
rule snodb_boxes_scores :
	 input: 
	 output: 
	 script: 

#=====================================================
#
#=====================================================
 
rule mismatch_boxes :
	 input: 
	 output: 
	 script: 

#=====================================================
#
#=====================================================

 
rule distance :
	 input: 
	 output: 
	 script: 

#=====================================================
#
#=====================================================
 
rule pairing_distributions :
	 input: 
	 output: 
	 script: 

#=====================================================
#
#=====================================================
 
rule generated_guide_scores_my_resul :
	 input: 
	 output: 
	 script: 

#=====================================================
#
#=====================================================
 
rule negative_examples :
	 input: 
	 output: 
	 script: 

#=====================================================
#
#=====================================================

 
rule learning_random_forest :
	 input: 
	 output: 
	 script: 

#=====================================================
#
#=====================================================

 
rule other_learning :
	 input: 
	 output: 
	 script: 
	   
#=====================================================
#
#=====================================================

 
rule testingvalidating :
	 input: 
	 output: 
	 script: 

#=====================================================
#
#=====================================================

 
rule fitting_distance_with_distribution :
	 input: 
	 output: 
	 script: 

#=====================================================
#
#=====================================================
 
rule other :
	 input: 
	 output: 
	 script: 

#=====================================================
#
#=====================================================

 
