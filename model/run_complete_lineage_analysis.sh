#!/bin/bash

# Complete Lineage Analysis Workflow
#
# This script runs the full analysis pipeline:
# 1. Calculate RLIPP scores
# 2. Identify high-RLIPP lineages
# 3. Create network visualizations
# 4. Analyze sample and gene drivers
#
# Usage:
#   bash run_complete_lineage_analysis.sh

set -e  # Exit on error

echo "======================================================================="
echo "COMPLETE LINEAGE ANALYSIS WORKFLOW"
echo "======================================================================="
echo ""

# ============================================================================
# CONFIGURATION - Modify these paths for your data
# ============================================================================

# Input files (REQUIRED - update these paths!)
MODEL="results/model.pt"
ONTOLOGY="data/ontology.txt"
GENE2ID="data/gene2id.txt"
TEST_DATA="data/test.pt"
OBO_FILE="data/go-basic.obo"  # Optional, for readable term names

# Model parameters (must match training configuration)
GENOTYPE_HIDDENS=6

# Output directory
OUTPUT_DIR="results/lineage_analysis"
mkdir -p $OUTPUT_DIR

# Analysis parameters
MIN_CONSISTENCY=0.5      # Minimum fraction of terms with expected RLIPP sign
TOP_LINEAGES_TO_FIND=20  # Number of top lineages to identify
TOP_LINEAGES_TO_ANALYZE=5  # Number to do detailed driver analysis on
TOP_SAMPLES=20           # Number of top samples to analyze per lineage
NETWORK_VIZ_TOP=3        # Number of individual network diagrams to create
NETWORK_VIZ_COMBINED=5   # Number of lineages in combined network

# ============================================================================
# STEP 1: Calculate RLIPP Scores
# ============================================================================

echo "Step 1: Calculating RLIPP scores..."
echo "-------------------------------------------------------------------"

RLIPP_OUTPUT="$OUTPUT_DIR/rlipp_scores.tsv"

if [ -f "$RLIPP_OUTPUT" ]; then
    echo "RLIPP scores already exist at: $RLIPP_OUTPUT"
    echo "Skipping calculation (delete file to recalculate)"
else
    python calculate_rlipp.py \
        -model "$MODEL" \
        -onto "$ONTOLOGY" \
        -gene2id "$GENE2ID" \
        -test "$TEST_DATA" \
        -genotype_hiddens $GENOTYPE_HIDDENS \
        -output "$RLIPP_OUTPUT" \
        -predictor_epochs 100 \
        -batch_size 64

    echo "✓ RLIPP scores saved to: $RLIPP_OUTPUT"
fi

echo ""

# ============================================================================
# STEP 2: Identify High-RLIPP Lineages
# ============================================================================

echo "Step 2: Identifying high-RLIPP lineages..."
echo "-------------------------------------------------------------------"

LINEAGES_PREFIX="$OUTPUT_DIR/lineages"

python analyze_rlipp_lineages.py \
    -rlipp "$RLIPP_OUTPUT" \
    -onto "$ONTOLOGY" \
    -output "$LINEAGES_PREFIX" \
    -direction both \
    -top_n $TOP_LINEAGES_TO_FIND \
    -min_consistency $MIN_CONSISTENCY \
    -min_length 2 \
    -visualize_top 5

echo "✓ Lineage analysis complete"
echo "  - Positive lineages: ${LINEAGES_PREFIX}_positive.tsv"
echo "  - Negative lineages: ${LINEAGES_PREFIX}_negative.tsv"
echo ""

# ============================================================================
# STEP 3: Create Network Visualizations
# ============================================================================

echo "Step 3: Creating network-style visualizations..."
echo "-------------------------------------------------------------------"

# Check if OBO file exists
OBO_PARAM=""
if [ -f "$OBO_FILE" ]; then
    OBO_PARAM="-obo $OBO_FILE"
    echo "Using GO term names from: $OBO_FILE"
else
    echo "Warning: OBO file not found, will use GO IDs"
fi

# Positive lineages
if [ -f "${LINEAGES_PREFIX}_positive.tsv" ]; then
    echo "Creating network diagrams for POSITIVE lineages..."
    python visualize_lineage_network.py \
        -lineages "${LINEAGES_PREFIX}_positive.tsv" \
        -rlipp "$RLIPP_OUTPUT" \
        -onto "$ONTOLOGY" \
        -output "$OUTPUT_DIR/network_positive" \
        -direction positive \
        -top_n $NETWORK_VIZ_TOP \
        -combined $NETWORK_VIZ_COMBINED \
        $OBO_PARAM

    echo "✓ Positive lineage networks created"
fi

# Negative lineages
if [ -f "${LINEAGES_PREFIX}_negative.tsv" ]; then
    echo "Creating network diagrams for NEGATIVE lineages..."
    python visualize_lineage_network.py \
        -lineages "${LINEAGES_PREFIX}_negative.tsv" \
        -rlipp "$RLIPP_OUTPUT" \
        -onto "$ONTOLOGY" \
        -output "$OUTPUT_DIR/network_negative" \
        -direction negative \
        -top_n $NETWORK_VIZ_TOP \
        -combined $NETWORK_VIZ_COMBINED \
        $OBO_PARAM

    echo "✓ Negative lineage networks created"
fi

echo ""

# ============================================================================
# STEP 4: Analyze Sample and Gene Drivers
# ============================================================================

echo "Step 4: Analyzing sample and gene drivers..."
echo "-------------------------------------------------------------------"

# Positive lineages
if [ -f "${LINEAGES_PREFIX}_positive.tsv" ]; then
    echo "Analyzing drivers for POSITIVE lineages..."
    python analyze_lineage_drivers.py \
        -model "$MODEL" \
        -onto "$ONTOLOGY" \
        -gene2id "$GENE2ID" \
        -test "$TEST_DATA" \
        -lineages "${LINEAGES_PREFIX}_positive.tsv" \
        -genotype_hiddens $GENOTYPE_HIDDENS \
        -output "$OUTPUT_DIR/drivers_positive" \
        -top_lineages $TOP_LINEAGES_TO_ANALYZE \
        -top_samples $TOP_SAMPLES

    echo "✓ Positive lineage drivers analyzed"
fi

# Negative lineages
if [ -f "${LINEAGES_PREFIX}_negative.tsv" ]; then
    echo "Analyzing drivers for NEGATIVE lineages..."
    python analyze_lineage_drivers.py \
        -model "$MODEL" \
        -onto "$ONTOLOGY" \
        -gene2id "$GENE2ID" \
        -test "$TEST_DATA" \
        -lineages "${LINEAGES_PREFIX}_negative.tsv" \
        -genotype_hiddens $GENOTYPE_HIDDENS \
        -output "$OUTPUT_DIR/drivers_negative" \
        -top_lineages $TOP_LINEAGES_TO_ANALYZE \
        -top_samples $TOP_SAMPLES

    echo "✓ Negative lineage drivers analyzed"
fi

echo ""

# ============================================================================
# SUMMARY
# ============================================================================

echo "======================================================================="
echo "ANALYSIS COMPLETE!"
echo "======================================================================="
echo ""
echo "All results saved to: $OUTPUT_DIR"
echo ""
echo "Output files:"
echo "-------------"
echo ""
echo "RLIPP Scores:"
echo "  $RLIPP_OUTPUT"
echo ""
echo "Lineage Analysis:"
echo "  ${LINEAGES_PREFIX}_positive.tsv - Top positive lineages"
echo "  ${LINEAGES_PREFIX}_negative.tsv - Top negative lineages"
echo "  ${LINEAGES_PREFIX}_positive_summary.png - Summary plots"
echo "  ${LINEAGES_PREFIX}_positive_lineage_*.png - Individual lineage plots"
echo ""
echo "Network Visualizations:"
echo "  $OUTPUT_DIR/network_positive_network_*.png - Individual networks"
echo "  $OUTPUT_DIR/network_positive_network_combined.png - Combined view"
echo ""
echo "Driver Analysis:"
echo "  $OUTPUT_DIR/drivers_positive_lineage*_samples.tsv - Sample contributions"
echo "  $OUTPUT_DIR/drivers_positive_lineage*_genes.tsv - Gene contributions"
echo "  $OUTPUT_DIR/drivers_positive_lineage*_samples.png - Sample plots"
echo "  $OUTPUT_DIR/drivers_positive_lineage*_genes.png - Gene plots"
echo ""
echo "Next Steps:"
echo "-----------"
echo "1. Review lineage TSV files to see top pathways"
echo "2. Examine network diagrams for pathway structure"
echo "3. Check sample TSV files to identify key samples"
echo "4. Review gene TSV files for important genes"
echo ""
echo "For interpretation help, see: LINEAGE_DRIVER_ANALYSIS.md"
echo "======================================================================="
