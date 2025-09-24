~~setup~~
~~matching~~
~~boxes_scores~~
~~snodb_boxes_scores~~
~~mismatch_boxes~~
~~distance~~
pairing_distributions
generated_guide_scores_my_resul
negative_examples
learning_random_forest
other_learning
testingvalidating
fitting_distance_with_distribution
other

ampliare dati(e negativi)-> guarda altre pubblicazioni (snorna guida pochi principalmente predizioni)# [ mancanza di validazioni ] TOPO -> fbk metilaz\tion sites
# ampliare dati(e negativi, altre sequenze !snoRNA allinterno degli introni)-> guarda altre pubblicazioni (snorna guida pochi principalmente predizioni)
# [ aggiungere altre specie -> tabelle un po cosi da sistemare ]
# analisi strutture secondarie?
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
# @wrap bunch of code that do shit as a single function

# # --- Dynamic Configuration Loading For Correct Use Of Yamz ---
#
# # 1. Load the entire YAML file into an R list
# # Adjust the path below to where your config file is located
# mock_config <- yaml::read_yaml("workflow_config.yaml")
#
# # 2. (Optional but recommended) Set working directory
# # Only do this IF you absolutely need to and IF the path is in the config file
# # setwd(mock_config$project_root) # <-- Remove this if possible (see below)
#
# # 3. Create mock snakemake object using the loaded config
# snakemake <- list(
#   input = mock_config$mock_input,
#   output = mock_config$mock_output,
#   config = mock_config$parameters
# )
# metodi ml nuovi |->  random forest comparison

