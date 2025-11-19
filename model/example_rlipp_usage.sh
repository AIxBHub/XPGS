#!/bin/bash

# Example script showing how to calculate RLIPP scores
# This is a template - update the paths to match your actual files

# ==============================================================================
# EXAMPLE 1: Basic RLIPP calculation with minimal parameters
# ==============================================================================

python model/calculate_rlipp.py \
    -model MODEL/my_trained_model.pt \
    -onto data/ontology.txt \
    -gene2id data/gene2ind.txt \
    -test data/test_data.pt \
    -output results/rlipp_scores.tsv

# ==============================================================================
# EXAMPLE 2: Full parameter specification
# ==============================================================================

python model/calculate_rlipp.py \
    -model MODEL/202410141234_model.pt \
    -onto PGS001990/PGS001990_GO_term_genes_NodeNorm.txt \
    -gene2id PGS001990/gene2ind_PGS001990.txt \
    -test PGS001990/test_data.pt \
    -output results/rlipp_scores_full.tsv \
    -genotype_hiddens 50 \
    -min_dropout_layer 2 \
    -dropout_fraction 0.1 \
    -batchsize 64 \
    -predictor_epochs 100 \
    -cuda 0

# ==============================================================================
# EXAMPLE 3: CPU-only mode (no GPU)
# ==============================================================================

python model/calculate_rlipp.py \
    -model MODEL/my_model.pt \
    -onto data/ontology.txt \
    -gene2id data/gene2ind.txt \
    -test data/test_data.pt \
    -output rlipp_cpu.tsv \
    -cuda -1  # Use CPU

# ==============================================================================
# EXAMPLE 4: Quick calculation with fewer epochs (faster but less accurate)
# ==============================================================================

python model/calculate_rlipp.py \
    -model MODEL/my_model.pt \
    -onto data/ontology.txt \
    -gene2id data/gene2ind.txt \
    -test data/test_data.pt \
    -output rlipp_quick.tsv \
    -predictor_epochs 50  # Faster calculation

# ==============================================================================
# EXAMPLE 5: After calculation, analyze results
# ==============================================================================

# View top 20 terms
echo "Top 20 terms by RLIPP score:"
head -n 21 results/rlipp_scores.tsv | column -t

# Count positive vs negative RLIPP
echo ""
echo "Distribution of RLIPP scores:"
awk 'NR>1 {if ($2>0) pos++; else if ($2<0) neg++; else zero++} END {print "Positive:", pos, "Negative:", neg, "Zero:", zero}' results/rlipp_scores.tsv

# Find GO terms with RLIPP > 0.5
echo ""
echo "Highly important terms (RLIPP > 0.5):"
awk 'NR>1 && $2>0.5 {print $1, $2}' results/rlipp_scores.tsv | column -t
