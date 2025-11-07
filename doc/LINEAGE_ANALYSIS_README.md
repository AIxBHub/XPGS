# RLIPP Lineage Analysis

## Overview

This tool identifies **lineages** (paths from leaf terms to root) in your ontology hierarchy that have consistently high positive or negative RLIPP scores. This helps you understand which biological pathways or processes are:

- **Important for prediction** (high positive RLIPP lineages)
- **Redundant** (high negative RLIPP lineages)

## What is a Lineage?

A **lineage** is a path through the ontology hierarchy from a leaf term (bottom) to the root term (top). For example:

```
Leaf Term (GO:0001234)
    ↓
Intermediate Term 1 (GO:0002345)
    ↓
Intermediate Term 2 (GO:0003456)
    ↓
Root Term (GO:0008150)
```

Each term along this path has an RLIPP score. A lineage with consistently high scores indicates a biologically important pathway.

## Key Concepts

### RLIPP Score Interpretation

- **Positive RLIPP** (>0): Term adds predictive power beyond its children → Important
- **Negative RLIPP** (<0): Term is redundant with its children → Less important
- **Zero RLIPP** (≈0): Term provides similar information to children

### Lineage Scoring Metrics

For each lineage, we calculate:

1. **Consistency**: Fraction of terms with the expected sign (positive or negative)
   - 1.0 = All terms have the expected sign
   - 0.5 = Half the terms have the expected sign

2. **Mean Score**: Average absolute RLIPP score for terms with expected sign

3. **Magnitude**: Sum of RLIPP scores in the expected direction

4. **Combined Score**: `Consistency × Magnitude × Mean Score`
   - This metric balances all three factors
   - Higher values = better lineages

## Quick Start

### 1. Basic Usage

```bash
# Edit the configuration
vim model/run_lineage_analysis.sh

# Run the analysis
./model/run_lineage_analysis.sh
```

### 2. Command Line Usage

```bash
python model/analyze_rlipp_lineages.py \
    -rlipp rlipp_scores.tsv \
    -onto data/ontology.txt \
    -output lineage_results \
    -direction both \
    -top_n 20 \
    -min_consistency 0.5 \
    -min_length 3
```

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-rlipp` | RLIPP scores TSV file (from calculate_rlipp.py) | Required |
| `-onto` | Ontology file (same as used in training) | Required |
| `-output` | Output prefix for results files | `lineages` |
| `-direction` | Direction to analyze: positive/negative/both | `both` |
| `-top_n` | Number of top lineages to report | 20 |
| `-min_consistency` | Minimum consistency score (0-1) | 0.5 |
| `-min_length` | Minimum lineage length (number of terms) | 2 |
| `-visualize_top` | Number of top lineages to visualize | 5 |

### Parameter Guidelines

**`-min_consistency`**:
- `0.8-1.0`: Very strict (almost all terms must have expected sign)
- `0.5-0.7`: Moderate (majority of terms must have expected sign)
- `0.3-0.5`: Lenient (allows more mixed lineages)

**`-min_length`**:
- `2`: Include even short paths
- `3-5`: Focus on moderate-length pathways
- `>5`: Only analyze deep hierarchical pathways

## Output Files

### 1. Lineage Data Files (TSV)

**File**: `{output}_positive.tsv` and `{output}_negative.tsv`

Contains all identified lineages with metrics:

| Column | Description |
|--------|-------------|
| `rank` | Rank by combined score |
| `leaf` | Starting leaf term |
| `root` | Ending root term |
| `length` | Number of terms in lineage |
| `consistency` | Fraction of terms with expected sign |
| `mean_score` | Average absolute RLIPP score |
| `median_score` | Median absolute RLIPP score |
| `magnitude` | Sum of RLIPP scores in expected direction |
| `combined_score` | Overall lineage score |
| `min_score` | Minimum RLIPP score in lineage |
| `max_score` | Maximum RLIPP score in lineage |
| `path` | Full path from leaf to root |
| `rlipp_scores` | RLIPP scores for each term |

### 2. Summary Plots

**File**: `{output}_positive_summary.png` and `{output}_negative_summary.png`

Four-panel summary showing:
- **Top left**: Distribution of lineage lengths
- **Top right**: Distribution of consistency scores
- **Bottom left**: Mean score vs consistency (colored by combined score)
- **Bottom right**: Top 10 lineages by combined score

### 3. Individual Lineage Visualizations

**Files**: `{output}_positive_lineage_1.png`, `{output}_positive_lineage_2.png`, etc.

Horizontal bar chart showing:
- Each term in the lineage (y-axis, leaf at top, root at bottom)
- RLIPP score for each term (x-axis)
- Green bars = positive RLIPP
- Red bars = negative RLIPP
- Lineage metrics in title

## Interpretation Guide

### High Positive RLIPP Lineages

**What it means:**
- These pathways are **biologically important** for prediction
- Each term adds unique information beyond its children
- The pathway as a whole captures important biological processes

**Example interpretation:**
```
Lineage: DNA repair → DNA damage response → Cell cycle regulation → Root
Consistency: 0.90 | Mean: 0.65 | Combined: 2.35

Interpretation: This DNA repair pathway is highly predictive. Each level
of organization (specific repair, damage response, cell cycle) adds
unique predictive power. This suggests DNA damage response is a key
biological process for your phenotype.
```

### High Negative RLIPP Lineages

**What it means:**
- These pathways are **redundant** with their child processes
- Parent terms don't add much beyond what children capture
- May indicate overly granular ontology structure

**Example interpretation:**
```
Lineage: Specific kinase → Kinase family → Phosphorylation → Root
Consistency: 0.85 | Mean: -0.45 | Combined: 1.85

Interpretation: The parent terms (kinase family, phosphorylation) don't
add predictive power beyond the specific kinase. The model learns all
it needs from the specific kinase level; higher abstractions are redundant.
```

## Example Workflow

### Step 1: Calculate RLIPP Scores

```bash
# First, calculate RLIPP scores from your trained model
python model/calculate_rlipp.py \
    -model MODEL/trained_model.pt \
    -onto data/ontology.txt \
    -gene2id data/gene2ind.txt \
    -test data/test.pt \
    -output rlipp_scores.tsv
```

### Step 2: Identify Important Lineages

```bash
# Analyze lineages with moderate strictness
python model/analyze_rlipp_lineages.py \
    -rlipp rlipp_scores.tsv \
    -onto data/ontology.txt \
    -output important_pathways \
    -direction positive \
    -min_consistency 0.7 \
    -min_length 3 \
    -top_n 20
```

### Step 3: Examine Results

```bash
# View top lineages
head -20 important_pathways_positive.tsv | column -t

# Open summary plot
open important_pathways_positive_summary.png

# Open top lineage visualization
open important_pathways_positive_lineage_1.png
```

### Step 4: Biological Interpretation

Look at the top lineages and ask:
1. **Do the pathways make biological sense?**
   - Are they related to your phenotype?
   - Do they represent known biology?

2. **What level of organization is most important?**
   - Are specific genes most important? (leaf terms)
   - Are high-level processes most important? (near root)

3. **Are there consistent themes?**
   - Multiple lineages in same biological process?
   - Suggests that process is central to your phenotype

## Advanced Usage

### Find Very Strict Positive Lineages

```bash
# Only report lineages where >90% of terms are positive
python model/analyze_rlipp_lineages.py \
    -rlipp rlipp_scores.tsv \
    -onto data/ontology.txt \
    -output strict_positive \
    -direction positive \
    -min_consistency 0.9 \
    -min_length 4
```

### Compare Positive and Negative Lineages

```bash
# Analyze both directions
python model/analyze_rlipp_lineages.py \
    -rlipp rlipp_scores.tsv \
    -onto data/ontology.txt \
    -output comparison \
    -direction both \
    -min_consistency 0.6

# Compare the results
echo "Top positive lineage:"
head -2 comparison_positive.tsv | tail -1

echo "Top negative lineage:"
head -2 comparison_negative.tsv | tail -1
```

### Focus on Deep Pathways Only

```bash
# Only analyze long lineages (>5 terms)
python model/analyze_rlipp_lineages.py \
    -rlipp rlipp_scores.tsv \
    -onto data/ontology.txt \
    -output deep_pathways \
    -direction positive \
    -min_length 6
```

## Tips for Interpretation

### 1. Start with High Consistency

Begin with `min_consistency=0.7` or higher to find the clearest patterns. Lower the threshold if you find too few lineages.

### 2. Consider Path Length

- **Short paths** (2-3 terms): More specific, focused predictions
- **Long paths** (5+ terms): Broad biological processes

### 3. Look for Convergence

If multiple lineages pass through the same intermediate terms, those terms are likely important biological hubs.

### 4. Cross-reference with Literature

The most important lineages should align with known biology for your phenotype. If they don't, this could indicate:
- Novel biological insights
- Data quality issues
- Model artifacts

### 5. Use Combined Score

The combined score balances all factors. Lineages with high combined scores are:
- Consistent (most terms have expected sign)
- High magnitude (strong RLIPP scores)
- High mean (individual terms are strongly predictive)

## Troubleshooting

### "No lineages found"

**Problem**: Script reports 0 lineages meeting criteria

**Solutions**:
- Lower `min_consistency` (try 0.3-0.5)
- Lower `min_length` (try 2)
- Check that RLIPP scores loaded correctly
- Verify ontology file has same format as training

### "All lineages very short"

**Problem**: Most lineages are only 2-3 terms

**Possible reasons**:
- Your ontology is shallow (not deeply hierarchical)
- RLIPP scores are inconsistent across depths
- Consider using `min_length=2` and accepting shorter paths

### "Mixed results in same lineage"

**Problem**: Lineages have both high positive and negative terms

**Interpretation**:
- This is normal! Real biology is complex
- Focus on lineages with high consistency (>0.7)
- Or lower expectations and use `min_consistency=0.5`

## Integration with Other Analyses

### Combine with RLIPP Visualization

```bash
# 1. Calculate RLIPP scores
python model/calculate_rlipp.py -model model.pt -onto onto.txt -gene2id genes.txt -test test.pt -output rlipp.tsv

# 2. Visualize overall distribution
python model/visualize_rlipp.py -input rlipp.tsv -output rlipp_viz

# 3. Identify important lineages
python model/analyze_rlipp_lineages.py -rlipp rlipp.tsv -onto onto.txt -output lineages

# 4. Compare:
#    - Which high RLIPP terms appear in top lineages?
#    - Are high-RLIPP terms clustered in specific lineages?
```

## Citation

If you use this lineage analysis in your research, please cite:

1. The VNN/DCell method:
   > Ma et al. (2018). "Using deep learning to model the hierarchical structure and function of a cell."
   > Nature Methods, 15(4), 290-298. PMID: 29505029

2. Your own work implementing this analysis

## Questions?

For issues or questions:
1. Check the output TSV files for detailed metrics
2. Visualize individual lineages to understand patterns
3. Try different parameter combinations
4. See [RLIPP_README.md](RLIPP_README.md) for RLIPP calculation details
