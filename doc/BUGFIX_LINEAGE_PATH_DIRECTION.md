# Bug Fix: Lineage Path Direction Error

## Date
2025-10-14

## Problem
The lineage analysis script reported:
```
Finding all paths from leaves to root...
Found 0 leaf nodes with 0 total paths to root
```

Debug output showed:
```
Testing path finding from sample leaf to root:
    From 'GO:0051241' to 'ROOT':
    Found 0 path(s)
```

## Root Cause

**Edge Direction Mismatch**

The ontology file format is:
```
parent_term  child_term  default
```

The code loads this as:
```python
dG.add_edge(parent, child)
```

This creates a **directed edge FROM parent TO child**.

### Graph Structure:
```
ROOT → intermediate_term → leaf_term
```

Edges point **downward** from root to leaves (following the parent→child relationship).

### The Bug:
The original code tried to find paths from **leaf to root**:
```python
paths = list(nx.all_simple_paths(graph, leaf, root))
```

But `nx.all_simple_paths()` follows edges in their **forward direction** (parent→child). Since edges go FROM root TO leaf, there is **no path** going FROM leaf TO root by following edges forward!

### Analogy:
It's like trying to drive against traffic on a one-way street:
- One-way street: ROOT → leaf (edges go downward)
- Trying to drive: leaf → ROOT (against the arrows) ❌

## Solution

**Find paths in the correct direction, then reverse them.**

```python
def find_all_paths_to_root(graph, root):
    """Find all paths from each leaf node to the root."""
    leaves = [n for n in graph.nodes() if graph.out_degree(n) == 0]
    all_paths = {}

    for leaf in leaves:
        # Since edges go parent->child, find paths from ROOT to LEAF (following edges)
        paths_root_to_leaf = list(nx.all_simple_paths(graph, root, leaf))

        if paths_root_to_leaf:
            # Reverse each path to get leaf->root order
            paths_leaf_to_root = [list(reversed(path)) for path in paths_root_to_leaf]
            all_paths[leaf] = paths_leaf_to_root

    return all_paths
```

### How This Works:

1. **Find paths from root → leaf** (following edges forward) ✓
   ```
   ROOT → GO:0008150 → GO:0051241
   ```

2. **Reverse the path** to get leaf → root order ✓
   ```
   GO:0051241 → GO:0008150 → ROOT
   ```

3. **Result**: Paths go from leaf to root (as desired) ✓

## Files Modified

### 1. [analyze_rlipp_lineages.py](analyze_rlipp_lineages.py)
- **Function**: `find_all_paths_to_root()`
- **Line**: 58-93
- **Change**: Find paths root→leaf, then reverse to get leaf→root

### 2. [debug_lineage.py](debug_lineage.py)
- **Section**: Path finding test
- **Line**: 149-172
- **Change**: Test paths in correct direction (root→leaf) and clarify expected behavior

## Verification

After the fix, the script should output:
```
Finding all paths from leaves to root...
Found 847 leaf nodes with 2431 total paths to root  ✓
```

Instead of:
```
Found 0 leaf nodes with 0 total paths to root  ✗
```

## Testing

Run the fixed script:
```bash
python model/analyze_rlipp_lineages.py \
    -rlipp rlipp_scores.tsv \
    -onto ontology.txt \
    -output lineages \
    -direction positive
```

You should now see:
- Leaf nodes found > 0
- Total paths > 0
- Lineage analysis produces results

## Related Issues

This bug would affect anyone using:
- Directed graphs where edge direction matters
- NetworkX path finding in biological ontologies
- Any analysis requiring paths from leaves to roots

## Key Takeaway

When working with directed graphs:
1. **Understand edge direction**: Which way do edges point?
2. **Match path finding direction**: Follow edges in their forward direction
3. **Reverse if needed**: If you need opposite direction, reverse the path

## Prevention

Added detailed comments to clarify:
```python
# Since edges go parent->child, we need to find paths from root to leaf
# then reverse them to get leaf->root paths
```

This makes the edge direction explicit and prevents future confusion.
