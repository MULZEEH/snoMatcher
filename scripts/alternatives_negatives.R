#!/usr/bin/env Rscript
# snoRNA Negative Examples Generator for ML Training
# Generates three types of negative sequences:
# 1) Random sequences matching length/GC-content (using oligoApp-like approach)
# 2) Real genomic intergenic sequences
# 3) Decoy sequences with motifs but lacking structure

library(Biostrings)
library(GenomicRanges)
library(BSgenome)
library(seqinr)
library(dplyr)

# ============================================================================
# PART 1: Random Sequences Matching Length and GC-Content Distribution
# ============================================================================

load(file = "results/intermediate/processed_info_box.RData")


pos_seq <- DNAStringSet(snodb_boxes[["DNA Sequence"]])
pos_seq

generate_random_sequences <- function(positive_seqs, n_random = 1000, 
                                      output_file = "random_negatives.fasta") {
  cat("Generating random sequences matching positive set properties...\n")
  positive_seqs <- snodb_boxes
  # Analyze positive sequences
  seq_lengths <- sapply(positive_seqs["DNA Sequence"], nchar)
  gc_contents <- letterFrequency(DNAStringSet(positive_seqs[["DNA Sequence"]]), letters = "GC", as.prob = TRUE)
  
  # Get distributions
  length_mean <- mean(seq_lengths)
  length_sd <- sd(seq_lengths)
  gc_mean <- mean(gc_contents)
  gc_sd <- sd(gc_contents)
  
  # random_seqs <- DNAStringSet()
  new_sequences_list <- list()
  new_sequences_name <- character(50)
  
  for (i in 1:n_random) {
    # Sample length from normal distribution
    target_length <- round(rnorm(1, mean = length_mean, sd = length_sd))
    target_length <- max(50, min(500, target_length))  # Constrain to reasonable range
    
    # Sample GC content
    target_gc <- rnorm(1, mean = gc_mean, sd = gc_sd)
    target_gc <- max(0.2, min(0.8, target_gc))  # Constrain between 20-80%
    
    # Generate sequence with target GC content
    n_gc <- round(target_length * target_gc)
    n_at <- target_length - n_gc
    
    bases <- sample(c(rep("G", floor(n_gc/2)), 
                      rep("C", ceiling(n_gc/2)),
                      rep("A", floor(n_at/2)), 
                      rep("T", ceiling(n_at/2))))
    
    # seq <- DNAString(paste(bases, collapse = ""))
    # names(seq) <- sprintf("random_neg_%04d|len=%d|gc=%.2f", i, target_length, target_gc)
    # random_seqs <- c(random_seqs, seq)
    seq <- DNAString(paste(bases, collapse = ""))
    
    new_sequences_name[i] <- sprintf("random_neg_%04d|len=%d|gc=%.2f", i, target_length, target_gc)
    
    new_sequences_list[[i]] <- seq 
  }
  
  random_seqs <- DNAStringSet(new_sequences_list)
  
  names(random_seqs) <- new_sequences_name
  
  # Write to file
  # writeXStringSet(random_seqs, filepath = output_file)
  
  cat(sprintf("Generated %d random sequences -> %s\n", n_random, output_file))
  
  return(random_seqs)
}
negatives <- generate_random_sequences(snodb_boxes, 500, "data/generated.fasta")
negatives

# synth <- as.character(negatives)
# pos <- as.character(pos_seq)
# common <- intersect(synth, pos)
# length(common)

#===
#1.5 synthetic data but with a grammar
#===



# ============================================================================
# PART 2: Real Genomic Intergenic Sequences -> non coding sequences? regions
# ============================================================================

extract_intergenic_sequences <- function(genome_bsgenome, 
                                         annotation_file,
                                         n_sequences = 1000,
                                         min_length = 50,
                                         max_length = 500,
                                         output_file = "intergenic_negatives.fasta") {
  cat("Extracting intergenic sequences from genome...\n")
  
  # Load genome (e.g., BSgenome.Hsapiens.UCSC.hg38)
  # genome <- getBSgenome(genome_bsgenome)
  
  # Load gene annotations
  # This assumes GTF/GFF format - adjust as needed
  if (file.exists(annotation_file)) {
    genes <- rtracklayer::import(annotation_file)
    genes <- genes[genes$type == "gene"]
  } else {
    stop("Annotation file not found. Provide GTF/GFF file.")
  }
  
  # Get intergenic regions
  cat("Identifying intergenic regions...\n")
  chromosomes <- unique(seqnames(genes))
  intergenic_ranges <- GRangesList()
  
  for (chr in chromosomes) {
    chr_genes <- genes[seqnames(genes) == chr]
    chr_genes <- reduce(chr_genes)  # Merge overlapping genes
    
    # Get gaps between genes
    gaps <- gaps(chr_genes)
    
    # Filter by length
    gaps <- gaps[width(gaps) >= min_length & width(gaps) <= max_length * 2]
    intergenic_ranges[[as.character(chr)]] <- gaps
  }
  
  intergenic_ranges <- unlist(intergenic_ranges)
  
  # Sample random intergenic regions
  if (length(intergenic_ranges) < n_sequences) {
    cat(sprintf("Warning: Only %d intergenic regions available\n", 
                length(intergenic_ranges)))
    n_sequences <- length(intergenic_ranges)
  }
  
  sampled_ranges <- sample(intergenic_ranges, n_sequences)
  
  # Extract sequences
  cat("Extracting sequences...\n")
  intergenic_seqs <- getSeq(genome, sampled_ranges)
  
  # Trim to target length if needed
  for (i in seq_along(intergenic_seqs)) {
    if (width(intergenic_seqs[i]) > max_length) {
      start_pos <- sample(1:(width(intergenic_seqs[i]) - max_length + 1), 1)
      intergenic_seqs[i] <- subseq(intergenic_seqs[i], start_pos, start_pos + max_length - 1)
    }
  }
  
  names(intergenic_seqs) <- sprintf("intergenic_neg_%04d|%s:%d-%d", 
                                    seq_along(intergenic_seqs),
                                    seqnames(sampled_ranges),
                                    start(sampled_ranges),
                                    end(sampled_ranges))
  
  # Write to file
  writeXStringSet(intergenic_seqs, filepath = output_file)
  cat(sprintf("Extracted %d intergenic sequences -> %s\n", 
              length(intergenic_seqs), output_file))
  
  return(intergenic_seqs)
}


# ============================================================================
# MAIN EXECUTION FUNCTION
# ============================================================================

main <- function(positive_fasta, 
                 genome_bsgenome = "BSgenome.Hsapiens.UCSC.hg38",
                 annotation_file = NULL,
                 n_random = 1000,
                 n_intergenic = 1000,
                 n_decoy = 500,
                 output_dir = "negative_examples") {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load positive sequences
  cat("Loading positive snoRNA sequences...\n")
  positive_seqs <- readDNAStringSet(positive_fasta)
  cat(sprintf("Loaded %d positive sequences\n", length(positive_seqs)))
  
  # Generate all three types of negatives
  cat("\n=== Type 1: Random Sequences ===\n")
  random_file <- file.path(output_dir, "random_negatives.fasta")
  random_seqs <- generate_random_sequences(positive_seqs, n_random, random_file)
  
  cat("\n=== Type 2: Intergenic Sequences ===\n")
  if (!is.null(annotation_file) && file.exists(annotation_file)) {
    intergenic_file <- file.path(output_dir, "intergenic_negatives.fasta")
    
    # Load genome
    if (requireNamespace(genome_bsgenome, quietly = TRUE)) {
      genome <- getBSgenome(genome_bsgenome)
      intergenic_seqs <- extract_intergenic_sequences(genome, annotation_file, 
                                                      n_intergenic, 
                                                      output_file = intergenic_file)
    } else {
      cat(sprintf("Warning: %s not installed. Skipping intergenic extraction.\n", 
                  genome_bsgenome))
    }
  } else {
    cat("Warning: No annotation file provided. Skipping intergenic extraction.\n")
    cat("Provide GTF/GFF file with annotation_file parameter.\n")
  }
  
  cat("\n=== Type 3: Decoy Sequences ===\n")
  decoy_file <- file.path(output_dir, "decoy_negatives.fasta")
  decoy_seqs <- generate_decoy_sequences(positive_seqs, n_decoy, decoy_file)
  
  # Summary statistics
  cat("\n=== Summary ===\n")
  cat(sprintf("Total negative examples generated: %d\n", 
              n_random + n_decoy))
  cat(sprintf("Output directory: %s\n", output_dir))
  
  # Create combined negative set
  all_negatives <- c(random_seqs, decoy_seqs)
  combined_file <- file.path(output_dir, "all_negatives.fasta")
  writeXStringSet(all_negatives, filepath = combined_file)
  cat(sprintf("Combined negatives file: %s\n", combined_file))
  
  return(list(
    random = random_seqs,
    decoy = decoy_seqs
  ))
}


# ============================================================================
# USAGE EXAMPLE
# ============================================================================

# Example usage:
# negatives <- main(
#   positive_fasta = "positive_snornas.fasta",
#   genome_bsgenome = "BSgenome.Hsapiens.UCSC.hg38",
#   annotation_file = "gencode.v38.annotation.gtf",
#   n_random = 1000,
#   n_intergenic = 1000,
#   n_decoy = 500,
#   output_dir = "negative_examples"
# )

# For testing without genome data:
# negatives <- main(
#   positive_fasta = "positive_snornas.fasta",
#   n_random = 1000,
#   n_decoy = 500,
#   output_dir = "negative_examples"
# )

cat("Script loaded. Run main() function with your parameters.\n")
cat("Example: negatives <- main('positive_snornas.fasta', output_dir='negatives')\n")