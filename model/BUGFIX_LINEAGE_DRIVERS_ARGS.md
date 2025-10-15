# Bug Fix: MinimalArgs Missing Attributes

## Date
2025-10-15

## Problem

When running `analyze_lineage_drivers.py`, got error:
```
AttributeError: 'MinimalArgs' object has no attribute 'lr'
```

## Root Cause

The `MinimalArgs` class only provided a subset of attributes required by `TrainingDataWrapper.__init__()`:

**Required by TrainingDataWrapper:**
- Model architecture: `genotype_hiddens`, `min_dropout_layer`, `dropout_fraction`
- Training params: `lr`, `wd`, `alpha`, `epoch`, `batchsize`
- Output config: `modeldir`, `delta`, `metric_output`
- Dataset config: `train`, `test`, `testsetratio`
- Optimization: `optimize`

**What MinimalArgs had:**
- Only: `onto`, `gene2id`, `genotype_hiddens`, `min_dropout_layer`, `dropout_fraction`

**Missing:**
- `lr`, `wd`, `alpha`, `epoch`, `batchsize`
- `modeldir`, `delta`, `metric_output`
- `train`, `test`, `testsetratio`
- `optimize`

## Solution

Added all required attributes to `MinimalArgs` with dummy values (since they're not used for inference):

```python
class MinimalArgs:
    def __init__(self):
        # Required files
        self.onto = args.onto
        self.gene2id = args.gene2id

        # Model architecture
        self.genotype_hiddens = args.genotype_hiddens
        self.min_dropout_layer = 2
        self.dropout_fraction = 0.3

        # Training hyperparameters (dummy values, not used for inference)
        self.lr = 0.001
        self.wd = 0.001
        self.alpha = 0.3
        self.epoch = 100
        self.batchsize = 64

        # Output configuration (dummy values, not used)
        self.modeldir = '.'
        self.delta = 0.001
        self.metric_output = 'metrics.tsv'

        # Dataset configuration (dummy values, not used)
        self.train = None
        self.test = None
        self.testsetratio = 0.2

        # Optimization (not used for inference)
        self.optimize = 0
```

## Why Dummy Values Are OK

The `analyze_lineage_drivers.py` script only uses `TrainingDataWrapper` for:
1. Loading ontology structure (via `self.load_ontology()`)
2. Loading gene-to-ID mapping (via `utils.load_mapping()`)

The training hyperparameters and output configuration are **NOT** used because:
- No training is performed (model already trained)
- No optimization is performed (inference only)
- Output files are managed by the driver analysis script itself

## Files Modified

- **[analyze_lineage_drivers.py](analyze_lineage_drivers.py)** (lines 378-409)

## Verification

After fix, the script should run successfully:
```bash
python analyze_lineage_drivers.py \
    -model model.pt \
    -onto ontology.txt \
    -gene2id gene2id.txt \
    -test test.pt \
    -lineages lineages_positive.tsv \
    -output drivers
```

Expected output:
```
======================================================================
LINEAGE DRIVER ANALYSIS
======================================================================

Initializing data wrapper...
Total number of genes = 574
Ontology loaded: 1234 terms, root=ROOT
Loading model from: model.pt
...
```

No more AttributeError!

## Alternative Solution (Future)

A cleaner solution would be to refactor `TrainingDataWrapper` into two classes:
1. **OntologyLoader** - Only loads ontology and gene mappings (no training params)
2. **TrainingDataWrapper** - Extends OntologyLoader with training configuration

Then `analyze_lineage_drivers.py` would use `OntologyLoader` instead of creating dummy args.

However, this requires modifying the core VNN codebase, so the current fix with dummy values is simpler and safer.
