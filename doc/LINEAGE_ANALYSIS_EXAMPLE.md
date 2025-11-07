# RLIPP Lineage Analysis - Example Walkthrough

## Scenario

You've trained a VNN model to predict drug response, and you want to identify which biological pathways are most important for the prediction.

## Step-by-Step Example

### Step 1: Calculate RLIPP Scores (if not done already)

```bash
python model/calculate_rlipp.py \
    -model MODEL/drug_response_model.pt \
    -onto data/GO_ontology.txt \
    -gene2id data/gene2ind.txt \
    -test data/test_patients.pt \
    -output drug_response_rlipp.tsv
```

**Output**: `drug_response_rlipp.tsv` with RLIPP scores for each GO term

---

### Step 2: Run Lineage Analysis

```bash
python model/analyze_rlipp_lineages.py \
    -rlipp drug_response_rlipp.tsv \
    -onto data/GO_ontology.txt \
    -output drug_pathways \
    -direction both \
    -top_n 20 \
    -min_consistency 0.6 \
    -min_length 3 \
    -visualize_top 5
```

**What this does**:
- Finds lineages where ≥60% of terms have consistent RLIPP sign
- Requires at least 3 terms per lineage
- Reports top 20 lineages for both positive and negative directions
- Creates visualizations for top 5 in each direction

---

### Step 3: Examine Results

#### 3.1 Check Console Output

```
================================================================================
RLIPP LINEAGE ANALYSIS
================================================================================

Loading RLIPP scores from: drug_response_rlipp.tsv
Loaded RLIPP scores for 1,234 terms

Loading ontology from: data/GO_ontology.txt
Loaded ontology with 1,234 terms and 2,456 edges
Root term: GO:0008150

Finding all paths from leaves to root...
Found 567 leaf nodes with 1,234 total paths to root

================================================================================
Analyzing POSITIVE lineages
================================================================================

Found 45 lineages meeting criteria:
  - Minimum consistency: 0.6
  - Minimum length: 3

Top 5 positive lineages:

1. Leaf: GO:0006281_DNA_repair
   Length: 5 terms
   Consistency: 0.80
   Mean score: 0.623
   Combined score: 2.494
   Path: GO:0006281 -> GO:0006974 -> GO:0051716 -> GO:0050789 -> GO:0008150

2. Leaf: GO:0045786_negative_regulation_of_cell_cycle
   Length: 6 terms
   Consistency: 0.83
   Mean score: 0.589
   Combined score: 2.458
   Path: GO:0045786 -> GO:0010564 -> GO:0051726 -> ...

...
```

#### 3.2 View Summary Plots

```bash
open drug_pathways_positive_summary.png
```

The summary shows:
- Most lineages are 4-6 terms long
- Consistency peaks around 0.7-0.8
- Clear separation between high and low combined scores

#### 3.3 View Top Lineage Details

```bash
open drug_pathways_positive_lineage_1.png
```

**Example visualization shows**:

```
GO:0006281_DNA_repair                    [████████████] 0.623
GO:0006974_DNA_damage_response          [██████████  ] 0.545
GO:0051716_cellular_response_to_stress  [███████     ] 0.401
GO:0050789_regulation_of_process        [████        ] 0.234
GO:0008150_biological_process           [██          ] 0.156

Consistency: 0.80 | Mean: 0.623 | Combined Score: 2.494
```

**Interpretation**: DNA repair pathway is highly important! Each level adds predictive power:
- Specific DNA repair mechanisms (0.623)
- Broader DNA damage response (0.545)
- General stress response (0.401)
- Even high-level regulation adds information (0.234)

---

### Step 4: Examine the Data File

```bash
head -10 drug_pathways_positive.tsv | column -t -s $'\t'
```

**Output**:
```
rank  leaf                   root          length  consistency  mean_score  combined_score  path
1     GO:0006281            GO:0008150    5       0.800        0.623       2.494          GO:0006281 -> GO:0006974 -> ...
2     GO:0045786            GO:0008150    6       0.833        0.589       2.458          GO:0045786 -> GO:0010564 -> ...
3     GO:0042493            GO:0008150    4       0.750        0.701       2.103          GO:0042493 -> GO:0001775 -> ...
```

---

### Step 5: Biological Interpretation

#### Top Positive Lineages (Important Pathways)

From your results, the top 5 positive lineages are:

1. **DNA repair pathway** (score: 2.494)
   - Makes biological sense for drug response
   - DNA repair defects affect chemotherapy sensitivity
   - Each hierarchical level adds unique information

2. **Cell cycle regulation** (score: 2.458)
   - Also biologically relevant
   - Cell cycle checkpoints affect drug efficacy
   - Hierarchical organization is informative

3. **Apoptosis response** (score: 2.103)
   - Expected for drug response
   - Cell death mechanisms vary by drug type

**Biological insight**: Your model has learned that drug response depends on hierarchical coordination of DNA damage response, from specific repair mechanisms up through stress response pathways.

#### Top Negative Lineages (Redundant Pathways)

If you also analyzed negative lineages:

1. **Metabolic enzyme lineage** (score: -1.876)
   - Specific enzyme → enzyme class → metabolism
   - Parent terms don't add beyond specific enzyme
   - Model only needs fine-grained enzyme information

**Biological insight**: For your drug response phenotype, knowing the specific metabolic enzyme is sufficient; broader metabolic categories are redundant.

---

### Step 6: Cross-Reference with Individual RLIPP Scores

```bash
# Look at RLIPP scores for terms in top lineage
grep -E 'GO:0006281|GO:0006974|GO:0051716' drug_response_rlipp.tsv
```

**Output**:
```
term                                    rlipp    term_power  children_power
GO:0006281_DNA_repair                  0.623    0.782       0.543
GO:0006974_DNA_damage_response         0.545    0.756       0.589
GO:0051716_cellular_response_to_stress 0.401    0.689       0.612
```

**Interpretation**:
- Each term has positive RLIPP (consistent with lineage analysis)
- DNA repair has highest score (most important level)
- Predictive power gradually decreases up the hierarchy (expected)

---

### Step 7: Report Findings

**Key Findings for Publication/Report**:

> "We identified 45 biological lineages with consistently high positive RLIPP scores (≥60% consistency, ≥3 terms). The top-ranked lineage was the DNA repair pathway (GO:0006281), with a combined score of 2.494 and 80% consistency. Each hierarchical level in this pathway—from specific DNA repair mechanisms to general stress response—contributed unique predictive information, suggesting that hierarchical coordination of DNA damage response is central to drug response phenotypes."

> "In contrast, metabolic enzyme pathways showed negative RLIPP lineages, indicating that specific enzyme-level information was sufficient and higher-level metabolic categories were redundant for prediction."

---

## Advanced Example: Finding Convergent Pathways

### Identify terms that appear in multiple top lineages

```bash
# Extract all paths from top 10 positive lineages
awk -F'\t' 'NR>1 && NR<=11 {print $12}' drug_pathways_positive.tsv | \
    tr ' -> ' '\n' | \
    sort | uniq -c | sort -rn | head -20
```

**Output**:
```
  8 GO:0050789_regulation_of_biological_process
  6 GO:0006974_DNA_damage_response
  5 GO:0051716_cellular_response_to_stress
  3 GO:0008150_biological_process
  ...
```

**Interpretation**:
- `GO:0050789` appears in 8 of top 10 lineages → key hub term
- Multiple lineages converge through DNA damage response
- Suggests this biological process is central to the phenotype

---

## Tips for Your Analysis

1. **Start with `min_consistency=0.6`**
   - Balanced between strict and lenient
   - Usually gives interpretable results

2. **Compare positive and negative lineages**
   - Positive = what's important
   - Negative = what's redundant
   - Both provide biological insights

3. **Cross-check with literature**
   - Do top lineages align with known biology?
   - Novel findings could be real discoveries or artifacts

4. **Iterate on parameters**
   - If too few results: lower `min_consistency` or `min_length`
   - If too many results: raise these parameters

5. **Focus on combined score**
   - Best overall metric
   - Balances consistency, magnitude, and mean score

---

## Complete Example Script

Save this as `run_complete_lineage_analysis.sh`:

```bash
#!/bin/bash

# Complete lineage analysis workflow

# 1. Calculate RLIPP scores
echo "Step 1: Calculating RLIPP scores..."
python model/calculate_rlipp.py \
    -model MODEL/my_model.pt \
    -onto data/ontology.txt \
    -gene2id data/gene2ind.txt \
    -test data/test.pt \
    -output my_rlipp.tsv

# 2. Analyze lineages
echo "Step 2: Analyzing lineages..."
python model/analyze_rlipp_lineages.py \
    -rlipp my_rlipp.tsv \
    -onto data/ontology.txt \
    -output my_lineages \
    -direction both \
    -min_consistency 0.6 \
    -min_length 3

# 3. Generate summary
echo "Step 3: Generating summary..."
echo "Top 5 Positive Lineages:" > lineage_summary.txt
head -6 my_lineages_positive.tsv >> lineage_summary.txt
echo "" >> lineage_summary.txt
echo "Top 5 Negative Lineages:" >> lineage_summary.txt
head -6 my_lineages_negative.tsv >> lineage_summary.txt

# 4. Open visualizations
open my_lineages_positive_summary.png
open my_lineages_positive_lineage_1.png

echo "Analysis complete! Check lineage_summary.txt for results."
```

Make it executable and run:
```bash
chmod +x run_complete_lineage_analysis.sh
./run_complete_lineage_analysis.sh
```
