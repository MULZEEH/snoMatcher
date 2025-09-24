# snoMatcher: snoRNA Feature Analysis Pipeline

A bioinformatics pipeline/tool idk for identifying and analyzing small nucleolar RNA (snoRNA) features in the human genome using machine learning techniques.

## Overview

snoMatcher is a Snakemake-based pipeline that processes snoRNA sequence data to:
- Analyze C/D box snoRNA structural features
- Score and classify snoRNA sequences
- Generate visualizations of sequence motifs and distributions
- Export analysis results in multiple formats


## Installation

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)
- Git

### Quick Start

1. **Clone the repository:**
```bash
git clone https://github.com/MULZEEH/snoMatcher.git
cd snoMatcher
```

2. **Create the conda environment:**
```bash
conda env create -f env/snomatcher.yaml
conda activate snomatcher
```

3. **Prepare your data:**
Place your input files in the `data/raw/` directory:
- `snoDB_data.xlsx`: snoRNA database
- `snoDB_rRNA_interactions.xlsx`: Methylation sites data
- `cd_boxes.tsv`: C/D box sequences
- `human_genome.fasta`: Human genome reference (NEED TO BE ADDED) 

4. **Configure the pipeline:**
Edit `config.yaml` to match your data paths and analysis preferences.

5. **Run the pipeline:**
```bash
# Dry run to check the workflow
snakemake --dry-run

# Run with 4 cores
snakemake --cores 4

# Run with custom configuration
snakemake --cores 4 --config generate_plots=true export_tables=false
```

## Project Structure

```
snoMatcher/
├── Snakefile                 # Main workflow definition
├── config.yaml               # Configuration parameters
├── README.md                 # This file
├── scripts/                  # Analysis scripts
│   ├── setup.R               # Data loading and preprocessing
│   ├── computing_scores.R    # Scoring system implementation
│   └── other.R               # Others
├── data/                     # Data directory
│   └── raw/                  # Raw input files
├── results/                  # Output directory (created by pipeline)
│   ├── intermediate/         # Intermediate data files
│   ├── plots/                # Generated visualizations
│   ├── tables/               # Exported data tables
│   └── end_report.html/pdf   #
└── env/                      # Conda environment
    └── r_analysis.yaml       # Environment specification
```

## Usage Examples

### Basic Analysis
```bash
# Run complete pipeline
snakemake --cores 4

# Generate only plots
snakemake --cores 4 --config export_tables=false

# Run without plots (faster)
snakemake --cores 4 --config generate_plots=false
```

### Partial Runs
```bash
# Run only data setup
snakemake --cores 4 results/intermediate/info_box.RData

# Generate specific plots
snakemake --cores 4 results/plots/all_logos.pdf

# Clean and restart
snakemake clean
snakemake --cores 4
```

## Output Files

### Intermediate Data
- `results/intermediate/info_box.RData`: Processed snoRNA data
- `results/intermediate/scores.RData`: Computed feature scores
- `results/intermediate/pairing_distributions.RData`: Pairing analysis

### Visualizations (if enabled)
- `results/plots/all_logos.pdf`: Sequence logos for C/D boxes
- `results/plots/guide_width.pdf`: Guide sequence length distribution
- `results/plots/box_mismatch_distribution.pdf`: Mismatch analysis
- `results/plots/motif_score_distribution.pdf`: Score distributions
- Various other plots 

### Data Tables (if enabled)
- `results/tables/pfm_*_box_snoDb_snoRNA.csv`: Position frequency matrices
- Various other analysis tables


### Debug Mode

To run R scripts interactively for debugging: (NEED TO REMOVE HARDCODED PART)
```r
# In RStudio, set working directory (for now):
source("scripts/setup.R")
# When debugging remember to always run:
rm(snakemake, envir = .GlobalEnv)
# before the if(!exists("snakemake")) statement
```
