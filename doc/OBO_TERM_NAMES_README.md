# Using GO Term Names in Visualizations

## Overview

All RLIPP visualization scripts now support using human-readable GO term names (e.g., "mitochondrion inheritance") instead of GO IDs (e.g., "GO:0000001") by providing a GO OBO file.

## What is an OBO File?

OBO (Open Biomedical Ontologies) is a standard format for representing biological ontologies. The Gene Ontology Consortium distributes GO terms in OBO format.

**Download the latest GO OBO file:**
```bash
# Download go-basic.obo (recommended - filtered, acyclic version)
wget http://purl.obolibrary.org/obo/go/go-basic.obo

# Or the full version
wget http://purl.obolibrary.org/obo/go.obo
```

## Updated Scripts

Three scripts now support OBO files:

1. **[visualize_rlipp.py](visualize_rlipp.py)** - RLIPP visualization
2. **[analyze_rlipp_lineages.py](analyze_rlipp_lineages.py)** - Lineage analysis
3. **[obo_parser.py](obo_parser.py)** - OBO file parser (new utility module)

## Usage

### Basic Usage (Without OBO - Uses GO IDs)

```bash
# Old way - shows GO IDs like "GO:0000001"
python model/visualize_rlipp.py -input rlipp_scores.tsv -output viz
```

### Enhanced Usage (With OBO - Uses Term Names)

```bash
# New way - shows names like "mitochondrion inheritance"
python model/visualize_rlipp.py \
    -input rlipp_scores.tsv \
    -output viz \
    -obo go-basic.obo
```

## Examples

### visualize_rlipp.py

```bash
# Visualize RLIPP scores with term names
python model/visualize_rlipp.py \
    -input rlipp_scores.tsv \
    -output rlipp_viz \
    -top_n 20 \
    -obo data/go-basic.obo
```

**Output:** Plots will show "DNA repair" instead of "GO:0006281"

### analyze_rlipp_lineages.py

```bash
# Analyze lineages with term names
python model/analyze_rlipp_lineages.py \
    -rlipp rlipp_scores.tsv \
    -onto ontology.txt \
    -output lineages \
    -direction both \
    -obo data/go-basic.obo
```

**Output:** Lineage visualizations will show readable pathway names

## What Gets Updated

### In visualize_rlipp.py

**Top Terms Plot:**
- **Before:** `GO:0006281`, `GO:0006974`, `GO:0051716`
- **After:** `DNA repair`, `DNA damage response`, `cellular response to stress`

**All plots with term labels** now use names instead of IDs

### In analyze_rlipp_lineages.py

**Individual Lineage Plots:**
- **Before:**
  ```
  GO:0006281
  GO:0006974
  GO:0051716
  GO:0008150
  ```

- **After:**
  ```
  DNA repair
  DNA damage response
  cellular response to stress
  biological process
  ```

**Summary Plots:**
- Bar chart labels use term names
- More readable and interpretable

## OBO Parser Module

The new `obo_parser.py` module provides utilities:

```python
from obo_parser import create_term_name_dict, get_term_name, format_term_label

# Load term names from OBO file
name_dict = create_term_name_dict('go-basic.obo')

# Get name for a GO ID
name = get_term_name('GO:0006281', name_dict)
print(name)  # "DNA repair"

# Format label with optional ID
label = format_term_label('GO:0006281', name_dict, include_id=True)
print(label)  # "GO:0006281: DNA repair"

label = format_term_label('GO:0006281', name_dict, include_id=False)
print(label)  # "DNA repair"
```

## Features

### Automatic Fallback

If OBO file is not provided or cannot be loaded, scripts automatically fall back to using GO IDs:

```
Warning: OBO file provided but could not be loaded: [error message]
Continuing with GO IDs...
```

### Smart Truncation

Long term names are automatically truncated to fit in visualizations:

```python
# Names longer than max_length are truncated
"very long term name that would not fit" → "very long term name th..."
```

### Handles Missing Terms

If a GO ID is not found in the OBO file, the script uses the ID:

```python
get_term_name('GO:9999999', name_dict)  # Returns 'GO:9999999'
```

### Skips Obsolete Terms

The parser automatically skips obsolete GO terms:

```
[Term]
id: GO:0000002
name: obsolete mitochondrial genome maintenance
is_obsolete: true
```

This term will not be included in the name dictionary.

## Complete Workflow Example

### Step 1: Download OBO File

```bash
cd /path/to/your/project
wget http://purl.obolibrary.org/obo/go/go-basic.obo
```

### Step 2: Calculate RLIPP Scores

```bash
python model/calculate_rlipp.py \
    -model MODEL/trained_model.pt \
    -onto data/ontology.txt \
    -gene2id data/gene2ind.txt \
    -test data/test.pt \
    -output rlipp_scores.tsv
```

### Step 3: Visualize with Term Names

```bash
# Create visualizations with readable names
python model/visualize_rlipp.py \
    -input rlipp_scores.tsv \
    -output results/rlipp_viz \
    -top_n 30 \
    -obo go-basic.obo
```

### Step 4: Analyze Lineages with Term Names

```bash
# Identify important pathways with readable names
python model/analyze_rlipp_lineages.py \
    -rlipp rlipp_scores.tsv \
    -onto data/ontology.txt \
    -output results/lineages \
    -direction positive \
    -min_consistency 0.7 \
    -obo go-basic.obo
```

### Step 5: View Results

```bash
# Open visualizations - now with readable names!
open results/rlipp_viz_top_terms.png
open results/lineages_positive_lineage_1.png
```

## Comparison: With vs Without OBO

### Without OBO (GO IDs Only)

**Top 5 Terms:**
```
GO:0006281    0.623
GO:0006974    0.545
GO:0051716    0.401
GO:0050789    0.234
GO:0008150    0.156
```
❌ **Hard to interpret** - Need to look up each ID

### With OBO (Term Names)

**Top 5 Terms:**
```
DNA repair                          0.623
DNA damage response                 0.545
cellular response to stress         0.401
regulation of biological process    0.234
biological process                  0.156
```
✅ **Easy to interpret** - Immediately understand biology

## Benefits

### For Presentations

- **More professional** - Audience doesn't see cryptic IDs
- **Self-explanatory** - No need to explain what GO:0006281 means
- **Better figures** - Publication-ready without manual editing

### For Analysis

- **Faster interpretation** - Spot important pathways immediately
- **Better communication** - Share results with non-bioinformaticians
- **Easier validation** - Quickly check if results make biological sense

### For Publications

- **Clearer figures** - Reviewers can understand at a glance
- **Less annotation needed** - Figures are self-documenting
- **Professional appearance** - Standard practice in the field

## Troubleshooting

### "Warning: obo_parser not found"

**Problem:** The obo_parser module is not in your Python path

**Solution:**
```bash
# Make sure obo_parser.py is in the same directory as visualization scripts
ls model/obo_parser.py  # Should exist

# Or add model directory to Python path
export PYTHONPATH="/path/to/model:$PYTHONPATH"
```

### "Could not load OBO file"

**Problem:** OBO file path is incorrect or file is corrupted

**Solution:**
```bash
# Check file exists
ls go-basic.obo

# Try re-downloading
wget http://purl.obolibrary.org/obo/go/go-basic.obo

# Use absolute path
python model/visualize_rlipp.py -input rlipp.tsv -obo /full/path/to/go-basic.obo
```

### "Term names not showing up"

**Problem:** OBO parameter not provided

**Solution:**
```bash
# Make sure to include -obo parameter
python model/visualize_rlipp.py -input rlipp.tsv -output viz -obo go-basic.obo
                                                             # ^^^^ Don't forget!
```

### Names are truncated

**Problem:** Long term names don't fit in visualization

**Solution:** This is intentional to prevent overlapping labels. The truncation preserves readability. You can:
- View full names in the TSV output files
- Adjust `max_length` parameter in the code if needed
- Increase figure size for more space

## Testing OBO Parser

Test the parser directly:

```bash
# Test OBO parser
python model/obo_parser.py go-basic.obo

# Output:
# Parsing go-basic.obo...
# Found 45123 non-obsolete terms
#
# First 5 terms:
#
# GO:0000001:
#   Name: mitochondrion inheritance
#   Namespace: biological_process
#   Parents: ['GO:0048308', 'GO:0048311']
# ...
```

## Advanced: Custom Term Labels

If you want to customize how term names appear, edit `obo_parser.py`:

```python
def format_term_label(term_id, name_dict, include_id=True, max_length=50):
    """Customize this function to change label format"""
    name = get_term_name(term_id, name_dict, max_length=max_length)

    # Option 1: ID only (original behavior)
    return term_id

    # Option 2: Name only (current default)
    return name

    # Option 3: Both ID and name
    if include_id and term_id != name:
        return f"{term_id}: {name}"
    else:
        return name

    # Option 4: Custom format
    return f"[{term_id}] {name}"
```

## File Formats

### OBO File Structure

```
[Term]
id: GO:0000001
name: mitochondrion inheritance
namespace: biological_process
def: "The distribution of mitochondria..." [references]
is_a: GO:0048308 ! organelle inheritance
is_a: GO:0048311 ! mitochondrion distribution
```

### Name Dictionary Structure

```python
{
    'GO:0000001': 'mitochondrion inheritance',
    'GO:0006281': 'DNA repair',
    'GO:0006974': 'DNA damage response',
    ...
}
```

## Summary

✅ **Download go-basic.obo once**
✅ **Add `-obo go-basic.obo` to visualization commands**
✅ **Enjoy readable, professional visualizations!**

No other changes needed - scripts automatically use term names when OBO file is provided.
