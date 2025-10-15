# Bug Fix: OBO File Not Used in Lineage Analysis

## Date
2025-10-15

## Problem

When running the complete pipeline via `run_complete_lineage_analysis.sh`, the lineage visualization plots showed GO IDs instead of readable term names, even though an OBO file was specified.

**Example:**
- Y-axis labels showed: `GO:0045597`, `GO:0045595`, `GO:0050793`, etc.
- Should have shown: "regulation of...", "positive regulation of...", etc.

## Root Cause

The `run_complete_lineage_analysis.sh` script had the OBO file configured at the top:
```bash
OBO_FILE="data/go-basic.obo"
```

But **Step 2** (Identify High-RLIPP Lineages) didn't pass the OBO file to `analyze_rlipp_lineages.py`:

```bash
# MISSING -obo parameter!
python analyze_rlipp_lineages.py \
    -rlipp "$RLIPP_OUTPUT" \
    -onto "$ONTOLOGY" \
    -output "$LINEAGES_PREFIX" \
    -direction both \
    -top_n $TOP_LINEAGES_TO_FIND \
    -min_consistency $MIN_CONSISTENCY \
    -min_length 2 \
    -visualize_top 5
    # No -obo flag here!
```

Meanwhile, **Step 3** (Network Visualizations) correctly included the OBO parameter.

## Solution

Added OBO parameter to Step 2:

```bash
# Check if OBO file exists for term names
OBO_PARAM=""
if [ -f "$OBO_FILE" ]; then
    OBO_PARAM="-obo $OBO_FILE"
    echo "Using GO term names from: $OBO_FILE"
else
    echo "Warning: OBO file not found, will use GO IDs"
fi

python analyze_rlipp_lineages.py \
    -rlipp "$RLIPP_OUTPUT" \
    -onto "$ONTOLOGY" \
    -output "$LINEAGES_PREFIX" \
    -direction both \
    -top_n $TOP_LINEAGES_TO_FIND \
    -min_consistency $MIN_CONSISTENCY \
    -min_length 2 \
    -visualize_top 5 \
    $OBO_PARAM  # ← Now includes OBO file!
```

## Files Modified

- **[run_complete_lineage_analysis.sh](run_complete_lineage_analysis.sh)** (lines 84-102)

## Verification

After this fix, when you run the complete pipeline:

```bash
bash run_complete_lineage_analysis.sh
```

You should see:
```
Step 2: Identifying high-RLIPP lineages...
-------------------------------------------------------------------
Using GO term names from: data/go-basic.obo

...

✓ Lineage analysis complete
```

And the lineage plots will show:
- **Before (wrong)**: Y-axis labels = `GO:0045597`, `GO:0045595`
- **After (correct)**: Y-axis labels = "regulation of cell differentiation", "positive regulation of cell differentiation"

## Which Files Are Affected

This fix ensures term names appear in:
1. `lineages_positive_lineage_1.png` through `lineages_positive_lineage_5.png`
2. `lineages_negative_lineage_1.png` through `lineages_negative_lineage_5.png`
3. `lineages_positive_summary.png`
4. `lineages_negative_summary.png`

All individual lineage plots and summary plots will now use readable term names instead of GO IDs!

## Re-run to Get Updated Plots

To regenerate with term names:

```bash
# Delete old plots
rm -f results/lineage_analysis/lineages_*.png

# Re-run the complete pipeline
bash run_complete_lineage_analysis.sh
```

The new plots will have readable term names on all axes.
