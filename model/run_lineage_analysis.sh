#!/bin/bash

# Script to analyze RLIPP lineages
#
# This script identifies lineages (paths from leaf terms to root) in the ontology
# that have consistently high positive or negative RLIPP scores.
#
# Usage: ./run_lineage_analysis.sh
#
# Before running, update the paths below to match your files.

# ==============================================================================
# CONFIGURATION - Update these paths
# ==============================================================================

# Path to RLIPP scores file (output from calculate_rlipp.py)
RLIPP_FILE="rlipp_scores.tsv"

# Path to ontology file (same as used in training)
ONTOLOGY_FILE="data/ontology.txt"

# Output prefix for results
OUTPUT_PREFIX="lineage_analysis"

# ==============================================================================
# ANALYSIS PARAMETERS
# ==============================================================================

# Direction to analyze: positive, negative, or both
DIRECTION="both"

# Number of top lineages to report
TOP_N=20

# Minimum consistency (fraction of terms with expected sign, 0-1)
# 0.5 = at least 50% of terms must have the expected sign
MIN_CONSISTENCY=0.5

# Minimum lineage length (number of terms in path)
MIN_LENGTH=3

# Number of top lineages to visualize individually
VISUALIZE_TOP=5

# ==============================================================================
# RUN ANALYSIS
# ==============================================================================

echo "========================================"
echo "RLIPP Lineage Analysis"
echo "========================================"
echo "RLIPP file: $RLIPP_FILE"
echo "Ontology: $ONTOLOGY_FILE"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Direction: $DIRECTION"
echo "========================================"

python model/analyze_rlipp_lineages.py \
    -rlipp $RLIPP_FILE \
    -onto $ONTOLOGY_FILE \
    -output $OUTPUT_PREFIX \
    -direction $DIRECTION \
    -top_n $TOP_N \
    -min_consistency $MIN_CONSISTENCY \
    -min_length $MIN_LENGTH \
    -visualize_top $VISUALIZE_TOP

echo ""
echo "========================================"
echo "Analysis complete!"
echo "Check output files: ${OUTPUT_PREFIX}_*.tsv and ${OUTPUT_PREFIX}_*.png"
echo "========================================"
