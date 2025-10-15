# Lineage Driver Analysis Guide

## Overview

This guide describes how to identify:
1. **Which test samples** contribute most to high/low RLIPP lineages
2. **Which genes** in those samples activate the lineage pathways

This analysis answers critical biological questions:
- Which patients/samples exhibit the strongest pathway activity?
- What genetic variants drive these pathway phenotypes?
- Which genes are most important for high-RLIPP biological processes?

---

## Three-Step Workflow

### Step 1: Identify High-RLIPP Lineages

First, run the lineage analysis to find pathways with consistently high/low RLIPP scores:

```bash
python model/analyze_rlipp_lineages.py \
    -rlipp results/rlipp_scores.tsv \
    -onto data/ontology.txt \
    -output results/lineages \
    -direction positive \
    -min_consistency 0.6 \
    -top_n 20
```

**Outputs:**
- `results/lineages_positive.tsv` - Top lineages ranked by combined score
- `results/lineages_positive_summary.png` - Summary plots
- `results/lineages_positive_lineage_*.png` - Individual lineage plots

### Step 2: Visualize Lineages as Networks

Create network-style visualizations (similar to pathway diagrams):

```bash
python model/visualize_lineage_network.py \
    -lineages results/lineages_positive.tsv \
    -rlipp results/rlipp_scores.tsv \
    -onto data/ontology.txt \
    -output results/network \
    -direction positive \
    -top_n 3 \
    -combined 5 \
    -obo data/go-basic.obo
```

**Outputs:**
- `results/network_positive_network_1.png` - Top lineage network diagram
- `results/network_positive_network_2.png` - 2nd lineage network
- `results/network_positive_network_3.png` - 3rd lineage network
- `results/network_positive_network_combined.png` - Combined view of top 5 lineages

**Visualization Features:**
- Nodes sized by RLIPP magnitude (larger = more important)
- Nodes colored by RLIPP value (blue scale for positive)
- Hierarchical layout (root at top, leaves at bottom)
- Human-readable term names (if OBO file provided)
- Highlighted edges for primary lineage path

### Step 3: Identify Sample and Gene Drivers

Analyze which samples and genes drive the lineage results:

```bash
python model/analyze_lineage_drivers.py \
    -model results/model.pt \
    -onto data/ontology.txt \
    -gene2id data/gene2id.txt \
    -test data/test.pt \
    -lineages results/lineages_positive.tsv \
    -genotype_hiddens 6 \
    -output results/drivers \
    -top_lineages 3 \
    -top_samples 20
```

**Outputs (per lineage):**

Sample Analysis:
- `results/drivers_lineage1_samples.tsv` - Sample activation scores
- `results/drivers_lineage1_samples.png` - 4-panel sample analysis:
  - Top samples by total activation
  - Mean vs max activation scatter
  - Activation distribution histogram
  - Activation vs true phenotype values

Gene Analysis:
- `results/drivers_lineage1_genes.tsv` - Gene importance scores
- `results/drivers_lineage1_genes.png` - 4-panel gene analysis:
  - Genes directly annotated to lineage terms
  - Most frequently activated genes
  - Gene frequency vs lineage relevance
  - Gene value distribution

---

## Understanding the Results

### Sample Contribution Metrics

For each test sample, we compute:

1. **Mean Activation**: Average absolute embedding activation across lineage terms
   - Indicates overall engagement with the pathway

2. **Max Activation**: Maximum activation for any term in lineage
   - Identifies samples with strong activation at specific pathway points

3. **Total Activation**: Sum of activations across all lineage terms
   - Overall pathway involvement score

**Top samples** have the highest total activation - these are the patients/cases that most strongly exhibit this pathway's characteristics.

### Gene Contribution Metrics

For genes in the top samples, we compute:

1. **Lineage Relevance**: Number of lineage terms the gene is directly annotated to
   - 0 = gene not directly in pathway (but may affect it)
   - 1+ = gene is directly annotated to pathway terms

2. **Sample Frequency**: Number of top samples where this gene is activated
   - High frequency = gene commonly variant/activated in pathway-driven samples

3. **Mean/Total Value**: Average and sum of gene values across samples
   - For mutation data: mutation burden
   - For expression: expression level

**Most important genes** combine:
- High lineage relevance (directly annotated)
- High sample frequency (commonly activated)
- Present in samples with high activation

---

## Interpreting Results

### Example: DNA Repair Lineage

Suppose you find a high positive RLIPP lineage for DNA repair:

**Step 1 Output:**
```
Lineage 1:
  Leaf: GO:0006281 (DNA repair)
  Path: GO:0006281 -> GO:0006974 -> GO:0051716 -> ROOT
  Consistency: 0.85 (85% of terms positive)
  Mean RLIPP: 0.42
  Combined score: 2.34
```

**Interpretation**: This DNA repair pathway has consistently high RLIPP scores, meaning these terms add significant predictive power beyond their children.

**Step 2 Output:** Network diagram shows:
- Large nodes for GO:0006281 and GO:0006974 (high RLIPP)
- Hierarchical structure from specific (DNA repair) to general (ROOT)
- Blue gradient showing positive RLIPP values

**Step 3 Output:**

*Sample Analysis:*
```
Top 5 samples:
  Sample 42: activation=8.73
  Sample 127: activation=7.91
  Sample 8: activation=7.45
  ...
```

**Interpretation**: Samples 42, 127, and 8 have the strongest DNA repair pathway activation. These patients likely have characteristics related to DNA damage response.

*Gene Analysis:*
```
Top 5 lineage-relevant genes:
  BRCA1: in 15 samples, annotated to 3 terms
  BRCA2: in 12 samples, annotated to 3 terms
  ATM: in 18 samples, annotated to 2 terms
  TP53: in 20 samples, annotated to 1 term
  ...
```

**Interpretation**:
- BRCA1/BRCA2 are directly involved in the DNA repair pathway (annotated to 3 terms)
- These genes are mutated/activated in 12-15 of the top samples
- This explains why these samples have high DNA repair pathway activation
- TP53 appears in more samples (20) but fewer pathway terms - may be upstream regulator

---

## Biological Applications

### Use Case 1: Cancer Subtype Discovery

**Question**: Do different cancer samples activate different pathways?

**Approach**:
1. Run lineage analysis to find high-RLIPP pathways
2. Identify top samples for each pathway
3. Compare sample overlap between pathways
4. Check if samples cluster by pathway activation

**Insight**: Samples that activate "DNA repair" lineage may be different from those activating "cell cycle" lineage, suggesting distinct molecular subtypes.

### Use Case 2: Gene Prioritization

**Question**: Which genes should we validate for this pathway?

**Approach**:
1. Find genes with high lineage relevance (directly annotated)
2. Filter for genes in many top samples (high frequency)
3. Prioritize genes with both attributes

**Insight**: Genes like BRCA1 (annotated to 3 terms, in 15 samples) are better candidates than genes with low relevance or frequency.

### Use Case 3: Pathway Mechanism Discovery

**Question**: How do genetic variants lead to pathway activation?

**Approach**:
1. Identify top samples and their activated genes
2. Check which lineage terms these genes annotate to
3. Trace activation through the hierarchy

**Insight**: If BRCA1 mutations → high activation in "DNA double-strand break repair" (child term) → high activation in "DNA repair" (parent term), this confirms the mechanistic pathway.

### Use Case 4: Biomarker Development

**Question**: Which samples respond to DNA repair inhibitors?

**Approach**:
1. Find samples with high DNA repair lineage activation
2. Check if they have specific gene mutations (e.g., BRCA1/2)
3. Validate predictions experimentally

**Insight**: Samples with high DNA repair activation + BRCA mutations may be sensitive to PARP inhibitors.

---

## Advanced Analysis

### Comparing Positive vs Negative Lineages

Run analysis for both directions:

```bash
# Positive RLIPP lineages (terms add predictive power)
python model/analyze_lineage_drivers.py \
    -lineages results/lineages_positive.tsv \
    -output results/drivers_positive \
    ...

# Negative RLIPP lineages (terms are redundant)
python model/analyze_lineage_drivers.py \
    -lineages results/lineages_negative.tsv \
    -output results/drivers_negative \
    ...
```

**Interpretation**:
- **Positive lineages**: Important, non-redundant pathways
- **Negative lineages**: Redundant pathways (children already capture the information)

Samples driving positive lineages have unique pathway characteristics.
Samples driving negative lineages may have overlapping pathway activations.

### Cross-Lineage Sample Analysis

Identify samples that appear in multiple top lineages:

```bash
# Find common samples across lineages
comm -12 \
    <(tail -n +2 results/drivers_lineage1_samples.tsv | head -20 | cut -f1 | sort) \
    <(tail -n +2 results/drivers_lineage2_samples.tsv | head -20 | cut -f1 | sort)
```

**Interpretation**: Samples appearing in multiple high-RLIPP lineages may have:
- Complex phenotypes involving multiple pathways
- Pleiotropic genetic variants affecting multiple processes
- Coordinated pathway dysregulation

### Gene Network Analysis

Identify gene-gene relationships:

1. Extract top genes from multiple lineages
2. Check for:
   - Protein-protein interactions
   - Shared pathway annotations
   - Known regulatory relationships
3. Build gene network showing connections

**Interpretation**: Genes that co-occur in top samples across lineages may functionally interact.

---

## File Formats

### Input: lineages.tsv

```
rank  leaf      root  length  consistency  mean_score  combined_score  path                    rlipp_scores
1     GO:0006281  ROOT  5      0.85         0.42        2.34           GO:0006281 -> ... ROOT  0.521, 0.389, ...
```

### Output: drivers_lineage1_samples.tsv

```
sample_idx  mean_activation  max_activation  total_activation  true_value
42          2.91             4.52            8.73             0.823
127         2.64             4.01            7.91             0.745
```

### Output: drivers_lineage1_genes.tsv

```
gene    n_samples  mean_value  total_value  lineage_relevance  samples
BRCA1   15         1.23        18.45        3                  [42, 127, 8, ...]
BRCA2   12         1.08        12.96        3                  [42, 19, 103, ...]
```

---

## Troubleshooting

### Issue: "CUDA out of memory"

**Solution**: Process samples in batches:
```python
# Modify analyze_lineage_drivers.py to process in chunks
batch_size = 100
for i in range(0, len(x_test), batch_size):
    batch = x_test[i:i+batch_size]
    # Process batch...
```

### Issue: "No genes with lineage relevance"

**Cause**: Genes in test samples may not be annotated to the lineage terms.

**Solution**: This is expected if:
- Lineage is high-level (few direct gene annotations)
- Test genes are different from training genes
- Check "Most frequently activated genes" instead

### Issue: "All samples have similar activation"

**Cause**: Lineage may not be discriminative for test set.

**Solutions**:
- Try different lineages
- Check if test set is homogeneous
- Reduce `-top_samples` to focus on extreme cases

---

## Complete Example Workflow

```bash
# 1. Calculate RLIPP scores
python model/calculate_rlipp.py \
    -model results/model.pt \
    -onto data/ontology.txt \
    -gene2id data/gene2id.txt \
    -test data/test.pt \
    -genotype_hiddens 6 \
    -output results/rlipp_scores.tsv

# 2. Find top lineages
python model/analyze_rlipp_lineages.py \
    -rlipp results/rlipp_scores.tsv \
    -onto data/ontology.txt \
    -output results/lineages \
    -direction both \
    -min_consistency 0.5

# 3. Create network visualizations
python model/visualize_lineage_network.py \
    -lineages results/lineages_positive.tsv \
    -rlipp results/rlipp_scores.tsv \
    -onto data/ontology.txt \
    -output results/network \
    -obo data/go-basic.obo

# 4. Analyze sample and gene drivers
python model/analyze_lineage_drivers.py \
    -model results/model.pt \
    -onto data/ontology.txt \
    -gene2id data/gene2id.txt \
    -test data/test.pt \
    -lineages results/lineages_positive.tsv \
    -output results/drivers \
    -top_lineages 5

# 5. Review results
ls -lh results/drivers*
head -20 results/drivers_lineage1_samples.tsv
head -20 results/drivers_lineage1_genes.tsv
```

---

## Citation

If you use this lineage driver analysis in your research, please cite:

- The original VNN/DCell paper: Ma et al. (2018) Nature Methods
- RLIPP methodology: Ma et al. (2018) Nature Methods
- This analysis framework: [Your publication]

---

## Related Documentation

- [RLIPP_README.md](RLIPP_README.md) - RLIPP score calculation
- [LINEAGE_ANALYSIS_README.md](LINEAGE_ANALYSIS_README.md) - Finding high-RLIPP lineages
- [OBO_TERM_NAMES_README.md](OBO_TERM_NAMES_README.md) - Using GO term names
- [BUGFIX_LINEAGE_PATH_DIRECTION.md](BUGFIX_LINEAGE_PATH_DIRECTION.md) - Path finding bug fix
