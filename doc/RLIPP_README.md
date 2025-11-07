# RLIPP Score Calculation for VNN Models

This guide explains how to calculate RLIPP (Relative Local Improvement in Predictive Power) scores for trained VNN (Visible Neural Network) models.

## What is RLIPP?

RLIPP measures how much additional predictive information each subsystem/term in your hierarchical neural network provides beyond its child subsystems.

- **Positive RLIPP**: The term improves prediction relative to its children (important subsystem)
- **Negative RLIPP**: The term decreases predictive power (may be redundant)
- **Near-zero RLIPP**: The term provides similar information to its children

This metric helps identify which biological subsystems (e.g., GO terms) are most important for predicting your phenotype of interest.

## Reference

Based on the method described in:
> Ma et al. (2018). "Using deep learning to model the hierarchical structure and function of a cell." Nature Methods. PMID: 29505029

## Quick Start

### 1. Prerequisites

You need:
- A trained VNN model file (`.pt` file)
- The same ontology file used during training
- The same gene2id mapping file used during training
- A test dataset (`.pt` file)

### 2. Configure and Run

Edit [`run_rlipp.sh`](run_rlipp.sh) and update the following paths:

```bash
# Required paths
MODEL_PATH="MODEL/your_model.pt"              # Your trained model
ONTOLOGY_FILE="path/to/ontology.txt"          # Same as training
GENE2ID_FILE="path/to/gene2ind.txt"           # Same as training
TEST_DATA="path/to/test_data.pt"              # Test dataset

# Output
OUTPUT_FILE="rlipp_scores.tsv"                # Where to save results
```

**Important**: Make sure the model parameters match your training configuration:
- `GENOTYPE_HIDDENS` (default: 50)
- `MIN_DROPOUT_LAYER` (default: 2)
- `DROPOUT_FRACTION` (default: 0.1)

Then run:
```bash
./model/run_rlipp.sh
```

### 3. Command Line Usage

You can also run the Python script directly:

```bash
python model/calculate_rlipp.py \
    -model MODEL/your_model.pt \
    -onto ontology.txt \
    -gene2id gene2ind.txt \
    -test test_data.pt \
    -output rlipp_scores.tsv \
    -genotype_hiddens 50 \
    -min_dropout_layer 2 \
    -dropout_fraction 0.1 \
    -batchsize 64 \
    -predictor_epochs 100 \
    -cuda 0
```

## Output

The script generates a TSV file (`rlipp_scores.tsv` by default) with the following columns:

| Column | Description |
|--------|-------------|
| `term` | Name of the subsystem/GO term |
| `rlipp` | RLIPP score (relative improvement) |
| `term_power` | Predictive power using the term's embedding |
| `children_power` | Predictive power using only children's embeddings |
| `num_children` | Number of child terms |
| `layer` | Layer index in the hierarchy (0 = leaves) |

Results are sorted by RLIPP score (highest to lowest).

### Example Output

```
term                    rlipp    term_power  children_power  num_children  layer
GO:0006412             0.4523   0.7821      0.4532          15            3
GO:0019538             0.3891   0.7234      0.4892          8             2
GO:0006396            -0.1234   0.5123      0.5987          12            3
```

## Interpretation

### High RLIPP Scores (> 0.2)
- These subsystems are **important** for prediction
- They capture information not fully represented by their children
- Focus on these for biological interpretation

### Low/Negative RLIPP Scores (< 0)
- These subsystems may be **redundant**
- Their children already capture most of the predictive information
- May be candidates for pruning in model simplification

### Distribution Analysis
The script prints summary statistics:
- Number of terms with positive vs. negative RLIPP
- Top and bottom ranked terms
- Overall distribution statistics

## How RLIPP is Calculated

For each term in the hierarchy:

1. **Extract embeddings**: Get the hidden layer representation for the term and its children from the trained model

2. **Train term predictor**: Train a simple linear model to predict the phenotype using the term's embedding
   - Calculate correlation: `term_power`

3. **Train children predictor**: Train a linear model using concatenated children embeddings
   - Calculate correlation: `children_power`

4. **Calculate relative improvement**:
   ```
   RLIPP = (term_power - children_power) / max(|term_power|, |children_power|)
   ```

5. **Rank terms**: Sort all terms by RLIPP score to identify most important subsystems

## Parameters

### Required Parameters
- `-model`: Path to trained model (`.pt` file)
- `-onto`: Ontology file used in training
- `-gene2id`: Gene to ID mapping file
- `-test`: Test dataset (`.pt` file)

### Model Architecture Parameters (Must Match Training)
- `-genotype_hiddens`: Number of neurons per term (default: 50)
- `-min_dropout_layer`: Starting layer for dropout (default: 2)
- `-dropout_fraction`: Dropout probability (default: 0.1)

### Calculation Parameters
- `-batchsize`: Batch size for inference (default: 64)
- `-predictor_epochs`: Training epochs for predictive power calculation (default: 100)
- `-output`: Output filename (default: `rlipp_scores.tsv`)
- `-cuda`: GPU device ID (default: 0)

## Tips

### For Faster Calculation
- Increase `-batchsize` (if you have enough memory)
- Decrease `-predictor_epochs` (may be less accurate)
- Use GPU if available

### For More Accurate Results
- Increase `-predictor_epochs` (100-200 is usually sufficient)
- Ensure your test set is large and representative

### Memory Issues
- Decrease `-batchsize`
- Process in smaller chunks if needed

## Troubleshooting

### "Model parameter mismatch"
Make sure `-genotype_hiddens`, `-min_dropout_layer`, and `-dropout_fraction` match your training configuration.

### "CUDA out of memory"
Reduce `-batchsize` or use CPU (script will automatically fall back to CPU if CUDA is unavailable).

### "File not found"
Check that all paths are correct and files exist. Use absolute paths if relative paths cause issues.

## Citation

If you use this RLIPP implementation, please cite:

1. The original DCell/VNN paper:
   > Ma et al. (2018). "Using deep learning to model the hierarchical structure and function of a cell." Nature Methods, 15(4), 290-298. PMID: 29505029

2. Your own work implementing the VNN model

## Contact

For questions or issues, please contact the repository maintainer.
