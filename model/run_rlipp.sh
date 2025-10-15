#!/bin/bash

# Script to calculate RLIPP scores from a trained VNN model
#
# Usage: ./run_rlipp.sh
#
# Before running:
# 1. Make sure you have a trained model file (*.pt)
# 2. Update the paths below to match your setup
# 3. Make sure the parameters match your training configuration

# Load Python module (adjust if needed)
module load python/3.12
source .venv/bin/activate

# ==============================================================================
# CONFIGURATION - Update these paths to match your setup
# ==============================================================================

# Path to your trained model file
MODEL_PATH="MODEL/202512345678_model.pt"  # UPDATE THIS

# Path to your ontology file (same as used in training)
ONTOLOGY_FILE="PGS001990/PGS001990_GO_term_genes_NodeNorm.txt"  # UPDATE THIS

# Path to gene2id mapping file (same as used in training)
GENE2ID_FILE="PGS001990/gene2ind_PGS001990.txt"  # UPDATE THIS

# Path to test dataset (.pt file)
TEST_DATA="PGS001990/test_data.pt"  # UPDATE THIS

# Output file for RLIPP scores
OUTPUT_FILE="rlipp_scores.tsv"

# ==============================================================================
# TRAINING PARAMETERS - These must match your trained model
# ==============================================================================

GENOTYPE_HIDDENS=50      # Number of hidden neurons per term
MIN_DROPOUT_LAYER=2      # Minimum layer for dropout
DROPOUT_FRACTION=0.1     # Dropout fraction
BATCHSIZE=64            # Batch size for inference

# ==============================================================================
# RLIPP CALCULATION PARAMETERS
# ==============================================================================

PREDICTOR_EPOCHS=100    # Epochs to train predictor for each term
CUDA=0                  # GPU device ID (use 0 if you have GPU, ignored if no GPU)

# ==============================================================================
# RUN RLIPP CALCULATION
# ==============================================================================

echo "========================================"
echo "RLIPP Score Calculation"
echo "========================================"
echo "Model: $MODEL_PATH"
echo "Ontology: $ONTOLOGY_FILE"
echo "Test Data: $TEST_DATA"
echo "Output: $OUTPUT_FILE"
echo "========================================"

python model/calculate_rlipp.py \
    -model $MODEL_PATH \
    -onto $ONTOLOGY_FILE \
    -gene2id $GENE2ID_FILE \
    -test $TEST_DATA \
    -output $OUTPUT_FILE \
    -genotype_hiddens $GENOTYPE_HIDDENS \
    -min_dropout_layer $MIN_DROPOUT_LAYER \
    -dropout_fraction $DROPOUT_FRACTION \
    -batchsize $BATCHSIZE \
    -predictor_epochs $PREDICTOR_EPOCHS \
    -cuda $CUDA

echo ""
echo "========================================"
echo "RLIPP calculation complete!"
echo "Results saved to: $OUTPUT_FILE"
echo "========================================"
