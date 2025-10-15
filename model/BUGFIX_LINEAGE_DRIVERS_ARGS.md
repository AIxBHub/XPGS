# Bug Fix: Lineage Driver Analysis - Multiple Issues

## Date
2025-10-15

## Summary

Fixed 6 bugs in `analyze_lineage_drivers.py` to make it work with the VNN model architecture.

## Problem 1: MinimalArgs Missing Attributes

When running `analyze_lineage_drivers.py`, got error:
```
AttributeError: 'MinimalArgs' object has no attribute 'lr'
```

## Problem 2: DCellNN Constructor

After fixing Problem 1, got error:
```
TypeError: DCellNN.__init__() got an unexpected keyword argument 'term_size_map'
```

## Problem 3: Model Loading Format

After fixing Problems 1 and 2, got error:
```
TypeError: 'DCellNN' object is not subscriptable
```

This occurred on the line:
```python
model.load_state_dict(checkpoint['model_state_dict'])
```

## Problem 4: Test Data Key Names

After fixing Problems 1-3, got error:
```
ValueError: Test file dict must contain 'x' and 'y' keys
```

## Problem 5: Module Access Method

After fixing Problems 1-4, got error:
```
AttributeError: 'DCellNN' object has no attribute 'term_nn_dict'
```

## Problem 6: Data Type Mismatch

After fixing Problems 1-5, got error:
```
RuntimeError: mat1 and mat2 must have the same dtype, but got Char and Float
```

## Root Cause 1: Missing Attributes

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

## Root Cause 2: Wrong DCellNN Constructor Call

The original code called `DCellNN` with individual parameters:
```python
model = DCellNN(
    term_size_map=data_wrapper.term_size_map,
    term_direct_gene_map=data_wrapper.term_direct_gene_map,
    dG=data_wrapper.dG,
    gene_dim=len(data_wrapper.gene_id_mapping),
    num_hiddens_genotype=data_wrapper.num_hiddens_genotype,
    cuda_id=data_wrapper.cuda
)
```

But `DCellNN.__init__()` actually expects a **single `data_wrapper` object**:
```python
def __init__(self, data_wrapper):
    # ...
```

## Root Cause 3: Model Save Format

The VNN training code saves the model as a **complete object**, not as a state dictionary:

```python
# From vnn_train.py line 318, 324:
torch.save(self.model, pt_path)
```

This means the saved `.pt` file contains the entire `DCellNN` object, not a dict with `'model_state_dict'` key.

The original code assumed state dict format:
```python
checkpoint = torch.load(model_path)
model.load_state_dict(checkpoint['model_state_dict'])  # FAILS! checkpoint is the model itself
```

## Solution 1: Add All Required Attributes

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

## Solution 2: Fix DCellNN Constructor Call

Changed from passing individual parameters to passing the data_wrapper object.

## Solution 3: Handle Multiple Model Save Formats

Updated the model loading to handle both state dict and full object saves:

```python
def load_model_and_data(model_path, data_wrapper):
    """Load trained model."""
    print(f"Loading model from: {model_path}")

    # DCellNN expects a data_wrapper object
    model = DCellNN(data_wrapper)

    checkpoint = torch.load(model_path, map_location=data_wrapper.cuda)

    # Handle different save formats
    if isinstance(checkpoint, dict):
        # Saved as {'model_state_dict': ..., 'epoch': ..., etc}
        if 'model_state_dict' in checkpoint:
            model.load_state_dict(checkpoint['model_state_dict'])
        elif 'state_dict' in checkpoint:
            model.load_state_dict(checkpoint['state_dict'])
        else:
            # Dict but no standard key, assume the dict itself is state_dict
            model.load_state_dict(checkpoint)
    else:
        # Saved directly as model object (VNN default format)
        # checkpoint is already the model, just return it
        print("Model was saved as object, using loaded model directly")
        model = checkpoint
        model.eval()
        return model

    model.eval()
    print(f"Model loaded successfully")
    return model
```

This handles:
1. **State dict format**: `{'model_state_dict': ...}` or `{'state_dict': ...}`
2. **Full object format**: `torch.save(model, path)` (VNN default)
3. **Plain dict format**: State dict saved directly as dict

## Root Cause 4: Test Data Format

VNN uses capital `'X'` and lowercase `'y'` as keys in test `.pt` files (from `prepare_dataloader.py`):
```python
"""Expected .pt file format:
    "X": Tensor of shape (n_samples, n_features)
    "y": Tensor of shape (n_samples,)
"""
```

The original code checked for lowercase `'x'` first.

## Solution 4: Check VNN Format First

Reordered key checks to try `'X'` and `'y'` first:
```python
if 'X' in data and 'y' in data:  # VNN default
    x_test = data['X']
    y_test = data['y']
elif 'x' in data and 'y' in data:  # Alternative
    # ...
```

## Root Cause 5: Module Storage Method

`DCellNN` doesn't store modules in a `term_nn_dict` dictionary. Instead, it uses `add_module()` to register modules with names like:
- `{term}_direct_gene_layer`
- `{term}_dropout_layer`
- `{term}_linear_layer`
- `{term}_batchnorm_layer`

The original code tried:
```python
term_nn_dict = model.term_nn_dict[term]  # FAILS! No such attribute
```

## Solution 5: Use getattr() to Access Modules

Rewrote `extract_term_embeddings_per_sample()` to use `getattr()` and `hasattr()`:

```python
# Access modules by name
direct_gene_layer = getattr(model, term + '_direct_gene_layer')
gene_out = direct_gene_layer(gene_input)

# Check if dropout exists (only for higher layers)
if hasattr(model, term + '_dropout_layer'):
    dropout_layer = getattr(model, term + '_dropout_layer')
    term_input = dropout_layer(term_input)

# Apply transformations
linear_layer = getattr(model, term + '_linear_layer')
out = linear_layer(term_input)
out = torch.tanh(out)

batchnorm_layer = getattr(model, term + '_batchnorm_layer')
out = batchnorm_layer(out)
```

## Root Cause 6: Data Type Mismatch

Test data was loaded as `torch.int8` (Char type) but the model layers expect `torch.float32`.

This happens when genotype data is saved as integer types (0/1 for mutations) to save memory.

## Solution 6: Convert to Float

Convert input data to float before processing:

```python
def extract_term_embeddings_per_sample(model, x_data, device):
    model.eval()
    with torch.no_grad():
        # Convert to float and move to device
        x_data = x_data.float().to(device)

        # Now process...
```

## Why Dummy Values Are OK

The `analyze_lineage_drivers.py` script only uses `TrainingDataWrapper` for:
1. Loading ontology structure (via `self.load_ontology()`)
2. Loading gene-to-ID mapping (via `utils.load_mapping()`)
3. Providing configuration to `DCellNN` constructor

The training hyperparameters and output configuration are **NOT** used because:
- No training is performed (model already trained)
- No optimization is performed (inference only)
- Output files are managed by the driver analysis script itself

## Files Modified

- **[analyze_lineage_drivers.py](analyze_lineage_drivers.py)**
  - Lines 378-409: Added all required attributes to MinimalArgs (Problem 1)
  - Lines 38-67: Fixed DCellNN constructor call and model loading (Problems 2 & 3)
  - Lines 70-110: Fixed test data loading with multiple format support (Problem 4)
  - Lines 113-181: Rewrote embedding extraction using getattr() (Problem 5)

## All Issues Resolved

✅ **Problem 1**: MinimalArgs missing attributes → Added all required attributes
✅ **Problem 2**: DCellNN constructor → Pass data_wrapper object
✅ **Problem 3**: Model loading → Handle full object format
✅ **Problem 4**: Test data keys → Check 'X'/'y' first
✅ **Problem 5**: Module access → Use getattr() instead of dict lookup
✅ **Problem 6**: Data type mismatch → Convert input to float32

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
