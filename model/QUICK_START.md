# VNN Training Quick Start Guide

## Prerequisites

- Python 3.7+
- PyTorch
- NetworkX
- Required data files:
  - Ontology file (e.g., GO terms)
  - Gene-to-ID mapping file
  - Training data (.pt file)
  - Test data (.pt file)

## Basic Usage

### 1. View All Options

```bash
python model/run_vnn.py --help
```

### 2. Run Training

```bash
python model/run_vnn.py \
    -onto path/to/ontology.txt \
    -gene2id path/to/gene2ind.txt \
    -train path/to/train.pt \
    -test path/to/test.pt \
    -modeldir MODEL/ \
    -epoch 50
```

### 3. Monitor Training

**In console:**
- Real-time progress updates
- Epoch metrics (loss, correlation)
- Model saving notifications

**In log file:**
```bash
tail -f MODEL/logs/*_training.log
```

### 4. Check Results

```bash
# View configuration
cat MODEL/logs/*_config.json

# View metrics
cat MODEL/*_metrics_output.tsv

# Plot training curves
python -c "
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('MODEL/*_metrics_output.tsv', sep='\t')
df.plot(x='epoch', y=['train_corr', 'val_corr'])
plt.savefig('training_curves.png')
"
```

## Common Options

| Option | Description | Default |
|--------|-------------|---------|
| `-onto` | Ontology file (required) | - |
| `-gene2id` | Gene mapping file (required) | - |
| `-train` | Training data file (required) | - |
| `-test` | Test data file | None |
| `-epoch` | Number of epochs | 50 |
| `-lr` | Learning rate | 0.0002 |
| `-batchsize` | Batch size | 64 |
| `-genotype_hiddens` | Hidden neurons per term | 50 |
| `-alpha` | Auxiliary loss weight | 0.1 |
| `-delta` | Improvement threshold | 0.01 |
| `-modeldir` | Output directory | MODEL/ |

## File Structure

After training, you'll have:

```
MODEL/
├── logs/
│   ├── 202410141530_training.log    # Detailed training log
│   └── 202410141530_config.json     # Configuration
├── 202410141530_metrics_output.tsv  # Epoch metrics
└── 202410141530_model.pt            # Trained model
```

## Example: Full Training Run

```bash
# Create output directory
mkdir -p MODEL/

# Run training with custom parameters
python model/run_vnn.py \
    -onto data/GO_terms.txt \
    -gene2id data/gene2ind.txt \
    -train data/train_data.pt \
    -test data/test_data.pt \
    -modeldir MODEL/ \
    -genotype_hiddens 100 \
    -lr 0.0001 \
    -wd 0.0001 \
    -alpha 0.2 \
    -epoch 100 \
    -batchsize 128 \
    -delta 0.005

# Training will show progress like:
# INFO: Starting VNN Training Pipeline
# INFO: Using device: cuda
# INFO: Ontology loaded: 2345 terms, root=GO:0008150
# INFO: Model architecture created
# INFO:   Total parameters: 1,234,567
# ...
# INFO: Epoch 1/100
# INFO:   Train Loss: 0.3456 | Train Corr: 0.4567 | Val Loss: 0.3789 | Val Corr: 0.4123 | Time: 5.2s
# INFO:   *** Model saved (initial, val_corr=0.4123)
# ...
```

## Calculate RLIPP Scores (Post-Training)

After training, calculate RLIPP scores:

```bash
# Edit run_rlipp.sh with your model path
vim model/run_rlipp.sh

# Run RLIPP calculation
./model/run_rlipp.sh

# Visualize results
python model/visualize_rlipp.py -input rlipp_scores.tsv -output rlipp_viz
```

See [RLIPP_README.md](RLIPP_README.md) for details.

## Troubleshooting

### "CUDA out of memory"
- Reduce `-batchsize`
- Reduce `-genotype_hiddens`

### "Ontology has multiple roots"
- Check your ontology file
- Ensure single root term

### "No improvement in validation"
- Try lower learning rate (`-lr`)
- Try different alpha (`-alpha`)
- Check your data quality

### Training is slow
- Increase `-batchsize` (if memory allows)
- Use GPU if available
- Reduce `-genotype_hiddens`

## Getting Help

1. Check annotations in source files (all functions documented)
2. See [ANNOTATION_LOGGING_SUMMARY.md](ANNOTATION_LOGGING_SUMMARY.md) for details
3. Check log files for detailed debug information
4. See [RLIPP_README.md](RLIPP_README.md) for RLIPP analysis

## Next Steps

1. **Analyze Results**: Use metrics file to plot training curves
2. **Calculate RLIPP**: Identify important biological subsystems
3. **Interpret Model**: Examine term-level predictions
4. **Refine Parameters**: Adjust based on performance

## Citation

If using this code, please cite:

> Ma et al. (2018). "Using deep learning to model the hierarchical structure and function of a cell."
> Nature Methods, 15(4), 290-298. PMID: 29505029
