# Lineage Analysis Quick Start

## TL;DR - Run Everything

```bash
# Edit the paths in this script first!
bash run_complete_lineage_analysis.sh
```

This runs the complete workflow:
1. ✓ Calculate RLIPP scores
2. ✓ Find high-RLIPP lineages
3. ✓ Create network visualizations
4. ✓ Analyze sample and gene drivers

---

## Individual Steps

### 1. Find High-RLIPP Lineages (5 minutes)

```bash
python analyze_rlipp_lineages.py \
    -rlipp rlipp_scores.tsv \
    -onto ontology.txt \
    -output lineages \
    -direction positive
```

**Output**: `lineages_positive.tsv` - Top pathways ranked by RLIPP

### 2. Visualize as Networks (2 minutes)

```bash
python visualize_lineage_network.py \
    -lineages lineages_positive.tsv \
    -rlipp rlipp_scores.tsv \
    -onto ontology.txt \
    -output network \
    -obo go-basic.obo
```

**Output**: Network diagrams similar to your reference image

### 3. Find Sample/Gene Drivers (10 minutes)

```bash
python analyze_lineage_drivers.py \
    -model model.pt \
    -onto ontology.txt \
    -gene2id gene2id.txt \
    -test test.pt \
    -lineages lineages_positive.tsv \
    -output drivers
```

**Output**:
- Which samples drive each lineage
- Which genes are activated in those samples

---

## What You Get

### Network Visualizations

Similar to your image.png:
- **Nodes** = GO terms (sized by RLIPP, colored by value)
- **Edges** = parent-child relationships
- **Layout** = hierarchical (root at top)
- **Labels** = readable term names (if OBO provided)

### Sample Analysis

Answers: **"Which patients exhibit this pathway?"**

Output files per lineage:
- `drivers_lineage1_samples.tsv` - Ranked by pathway activation
- `drivers_lineage1_samples.png` - 4 plots showing sample distribution

**Top samples** = Patients with strongest pathway activity

### Gene Analysis

Answers: **"Which genes drive this pathway in those patients?"**

Output files per lineage:
- `drivers_lineage1_genes.tsv` - Genes ranked by importance
- `drivers_lineage1_genes.png` - 4 plots showing gene contributions

**Key genes** = Those with:
- High lineage relevance (annotated to pathway terms)
- High frequency (in many top samples)

---

## Example Results

### Lineage: DNA Repair Pathway

**Top 3 Samples:**
```
Sample 42:  activation = 8.73
Sample 127: activation = 7.91
Sample 8:   activation = 7.45
```
→ These patients have strongest DNA repair pathway activity

**Top 3 Genes:**
```
BRCA1: 15 samples, annotated to 3 pathway terms
BRCA2: 12 samples, annotated to 3 pathway terms
ATM:   18 samples, annotated to 2 pathway terms
```
→ These genes are mutated in top samples and directly in pathway

**Biological Insight**:
Samples 42, 127, 8 have BRCA1/BRCA2 mutations → DNA repair pathway activation → May respond to PARP inhibitors

---

## Troubleshooting

**Problem**: "Found 0 leaf nodes"
**Solution**: See [LINEAGE_TROUBLESHOOTING.md](LINEAGE_TROUBLESHOOTING.md)

**Problem**: Path finding error
**Solution**: Bug was fixed in [BUGFIX_LINEAGE_PATH_DIRECTION.md](BUGFIX_LINEAGE_PATH_DIRECTION.md)

**Problem**: Want GO term names instead of IDs
**Solution**: Add `-obo go-basic.obo` to visualization commands

**Problem**: CUDA out of memory
**Solution**: Reduce batch size or process samples in chunks

---

## Files You Need

### Required:
- `model.pt` - Trained VNN model
- `ontology.txt` - GO hierarchy file
- `gene2id.txt` - Gene name to ID mapping
- `test.pt` - Test data
- `rlipp_scores.tsv` - From `calculate_rlipp.py`

### Optional:
- `go-basic.obo` - For readable term names in plots

### Where to get OBO file:
```bash
wget http://purl.obolibrary.org/obo/go/go-basic.obo
```

---

## Complete Documentation

- **[LINEAGE_DRIVER_ANALYSIS.md](LINEAGE_DRIVER_ANALYSIS.md)** - Detailed guide with biological interpretation
- **[LINEAGE_ANALYSIS_README.md](LINEAGE_ANALYSIS_README.md)** - Finding high-RLIPP lineages
- **[RLIPP_README.md](RLIPP_README.md)** - RLIPP score calculation
- **[OBO_TERM_NAMES_README.md](OBO_TERM_NAMES_README.md)** - Using GO term names
- **[LINEAGE_TROUBLESHOOTING.md](LINEAGE_TROUBLESHOOTING.md)** - Fixing common issues

---

## Citation

If you use this analysis:

```
Ma et al. (2018) Using deep learning to model the hierarchical structure
and function of a cell. Nature Methods, 15(4), 290-298.
```

---

## Questions?

Run the debug script to diagnose issues:
```bash
python debug_lineage.py -onto ontology.txt -rlipp rlipp_scores.tsv
```

Check connectivity:
```bash
python analyze_ontology_connectivity.py -onto ontology.txt
```
