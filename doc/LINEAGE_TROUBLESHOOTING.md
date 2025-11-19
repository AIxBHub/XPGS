# Lineage Analysis Troubleshooting Guide

## Problem: "Found 0 leaf nodes with 0 total paths to root"

This error occurs when the lineage analysis script cannot find any leaf nodes in the ontology graph. This document explains the cause and how to fix it.

---

## Root Cause

The lineage analysis depends on finding **leaf nodes** (terms with no children) in the ontology graph to trace paths from leaves to the root. The "0 leaf nodes" error can occur for several reasons:

### 1. **Ontology Format Issue**
The ontology file must follow this exact format:
```
parent_term  child_term  default
```

Where:
- Column 1: Parent term ID
- Column 2: Child term ID
- Column 3: Must be the word `default` (indicates term-to-term relationship)

**Common mistakes:**
- Missing the `default` keyword
- Using a different relationship type
- Incorrect column ordering
- Extra or missing columns

### 2. **Circular References**
If the ontology has circular relationships (A → B → C → A), then no true leaf nodes exist.

### 3. **Wrong File Format**
Using the wrong ontology file (e.g., a gene annotation file instead of the term hierarchy file).

### 4. **Term ID Mismatch**
The term IDs in your ontology file don't match the term IDs in your RLIPP scores file.

**Example mismatch:**
- Ontology uses: `GO:0006281`, `GO:0006974`
- RLIPP file uses: `DNA repair`, `DNA damage response`

This happens when:
- The OBO file was used to convert IDs to names in one file but not the other
- Different preprocessing was applied to different files

---

## Diagnosis

### Step 1: Run the Debug Script

```bash
python model/debug_lineage.py \
    -onto path/to/ontology.txt \
    -rlipp path/to/rlipp_scores.tsv
```

This script will check:
1. ✓ Number of nodes and edges in ontology
2. ✓ Number of root nodes (should be 1)
3. ✓ Number of leaf nodes (should be > 0)
4. ✓ Sample term IDs from ontology
5. ✓ Sample term IDs from RLIPP file
6. ✓ Whether term IDs match between files
7. ✓ Graph connectivity
8. ✓ Path finding capability

### Step 2: Interpret Results

**If you see "Found 0 leaf node(s)":**
- Check: Are all nodes intermediate nodes?
- Cause: Ontology file likely has wrong format or circular references

**If you see "NO OVERLAP between ontology and RLIPP terms!":**
- Check: Do the sample terms look different?
- Cause: Term ID mismatch between files

**If you see "No path found from leaf to root":**
- Check: Is the graph disconnected?
- Cause: Ontology may have multiple disconnected components

---

## Solutions

### Solution 1: Verify Ontology File Format

Check the first few lines of your ontology file:

```bash
head -20 path/to/ontology.txt
```

**Correct format:**
```
GO:0008150  GO:0009987  default
GO:0008150  GO:0008152  default
GO:0009987  GO:0071840  default
```

**Incorrect formats:**
```
# Missing 'default'
GO:0008150  GO:0009987

# Wrong relationship type
GO:0008150  GO:0009987  gene_annotation

# Reversed order
default  GO:0008150  GO:0009987
```

### Solution 2: Check Term ID Consistency

Compare term IDs between files:

```bash
# Show sample terms from ontology
head -5 path/to/ontology.txt | awk '{print $1, $2}'

# Show sample terms from RLIPP file
head -5 path/to/rlipp_scores.tsv
```

**The term IDs must match exactly!**

If they don't match:
- **Option A**: Re-run RLIPP calculation with correct ontology file
- **Option B**: Convert term IDs to be consistent (use GO IDs everywhere or names everywhere)

### Solution 3: Verify the Ontology Structure

The ontology must be a **Directed Acyclic Graph (DAG)** with:
- Exactly **1 root node** (no parents)
- Many **leaf nodes** (no children)
- **Connected** (all nodes reachable from root)

Check this with the debug script output:
```
Found 1 root node(s)     ✓ Good
Found 847 leaf node(s)   ✓ Good
Found 0 leaf node(s)     ✗ Problem!
```

### Solution 4: Use Correct Input Files

Make sure you're using the **same ontology file** that was used during VNN training:

```bash
# The ontology file used for training
grep "^onto" run_vnn.sh  # Check your training script

# Use this EXACT same file for lineage analysis
python model/analyze_rlipp_lineages.py \
    -rlipp rlipp_scores.tsv \
    -onto THE_SAME_ONTOLOGY_FILE.txt \
    -output lineages
```

---

## Common Scenarios

### Scenario 1: Using GO OBO File Incorrectly

**Problem:**
```bash
# You trained with:
python run_vnn.py -onto data/go_hierarchy.txt ...

# But RLIPP used names from OBO file internally
# Now ontology has GO IDs but RLIPP has names
```

**Solution:**
The `-obo` parameter in visualization scripts is **only for display purposes**. It does NOT change the term IDs used internally. All analysis must use the **same term IDs** as the training ontology.

### Scenario 2: Wrong Ontology File

**Problem:**
```bash
# Using gene annotation file instead of hierarchy file
# Gene annotations: GO:0006281 BRCA1 gene_annotation
# Hierarchy file:   GO:0006950 GO:0006281 default
```

**Solution:**
Use the hierarchy file (term-to-term relationships with "default"), not the gene annotation file.

### Scenario 3: Preprocessed Ontology

**Problem:**
Some preprocessing scripts convert GO IDs to names, causing mismatches.

**Solution:**
Keep term IDs consistent throughout the entire pipeline:
1. Training: Use ontology with IDs
2. RLIPP calculation: Same ontology (IDs preserved)
3. Lineage analysis: Same ontology (IDs preserved)
4. Visualization: Use `-obo` flag to show names in plots

---

## Detailed Example

### Before (Broken):

```bash
# Training used GO IDs
cat data/ontology.txt
# GO:0008150  GO:0009987  default
# GO:0008150  GO:0008152  default

# But RLIPP file has names (mismatch!)
cat rlipp_scores.tsv
# term                    rlipp
# biological_process      0.123
# cellular_process        0.456

# Result: 0 leaf nodes found (no term ID overlap)
```

### After (Fixed):

```bash
# Training used GO IDs
cat data/ontology.txt
# GO:0008150  GO:0009987  default

# RLIPP file also has GO IDs (match!)
cat rlipp_scores.tsv
# term          rlipp
# GO:0008150    0.123
# GO:0009987    0.456

# Lineage analysis uses same ontology
python model/analyze_rlipp_lineages.py \
    -rlipp rlipp_scores.tsv \
    -onto data/ontology.txt \
    -output lineages

# Visualization uses OBO for display only
python model/visualize_rlipp.py \
    -input rlipp_scores.tsv \
    -output viz \
    -obo go-basic.obo  # Names shown in plots only

# Result: Found 847 leaf nodes ✓
```

---

## Quick Checklist

Before running lineage analysis, verify:

- [ ] Ontology file has format: `parent child default`
- [ ] Ontology has exactly 1 root node
- [ ] Ontology has multiple leaf nodes (> 0)
- [ ] Term IDs in ontology match term IDs in RLIPP file
- [ ] Using the SAME ontology file that was used for training
- [ ] Graph is connected (all nodes reachable from root)

---

## Getting Help

If you're still seeing "0 leaf nodes" after following this guide:

1. Run the debug script and save output:
   ```bash
   python model/debug_lineage.py \
       -onto ontology.txt \
       -rlipp rlipp_scores.tsv \
       > debug_output.txt
   ```

2. Check the debug output for specific issues

3. Look at the actual file contents:
   ```bash
   head -20 ontology.txt
   head -20 rlipp_scores.tsv
   ```

4. Verify file formats match the expected patterns described above

---

## Related Files

- **[debug_lineage.py](debug_lineage.py)** - Diagnostic tool for identifying issues
- **[analyze_rlipp_lineages.py](analyze_rlipp_lineages.py)** - Main lineage analysis script
- **[LINEAGE_ANALYSIS_README.md](LINEAGE_ANALYSIS_README.md)** - General lineage analysis documentation
- **[OBO_TERM_NAMES_README.md](OBO_TERM_NAMES_README.md)** - How to use GO term names in visualizations
