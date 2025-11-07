# Network Visualization Improvements

## Changes Made

Updated `visualize_lineage_network.py` to improve readability:

### 1. Show Full Term Names (No Truncation)

**Before:**
```python
# Terms were truncated
label = get_term_name(node, name_dict, max_length=20)  # "chromatin remod..."
label = get_term_name(node, name_dict, max_length=15)  # "positive reg..."
```

**After:**
```python
# Show complete term names
label = get_term_name(node, name_dict, max_length=None)  # "chromatin remodeling"
```

### 2. Increased Font Sizes

**Individual Network Plots:**
- Font size: **9 → 14** (55% larger)

**Combined Network Plot:**
- Font size: **8 → 12** (50% larger)
- Figure size: **(16, 12) → (20, 16)** (25% larger)

### 3. Files Modified

- **[visualize_lineage_network.py](visualize_lineage_network.py)**
  - Lines 128-133: Show full names in individual plots
  - Lines 173-179: Increased font size to 14
  - Lines 244-252: Show full names and larger figure for combined plot
  - Lines 301-307: Increased font size to 12

## Results

### Before:
- Truncated labels: "chromatin re...", "positive reg...", "sensory perc..."
- Small font (hard to read)
- Cramped layout

### After:
- Full labels: "chromatin remodeling", "positive regulation", "sensory perception"
- Large font (size 12-14, easy to read)
- Spacious layout (20×16 for combined plots)

## Re-generate Visualizations

To create new plots with these improvements:

```bash
python visualize_lineage_network.py \
    -lineages results/lineages_positive.tsv \
    -rlipp results/rlipp_scores.tsv \
    -onto data/ontology.txt \
    -output results/network_improved \
    -direction positive \
    -top_n 3 \
    -combined 5 \
    -obo data/go-basic.obo
```

This will create:
- `network_improved_positive_network_1.png` - Full term names, larger font
- `network_improved_positive_network_2.png` - Full term names, larger font
- `network_improved_positive_network_3.png` - Full term names, larger font
- `network_improved_positive_network_combined.png` - Full term names, larger figure

All labels will be readable without truncation!
