#!/usr/bin/env python3
"""
Debug script to diagnose lineage analysis issues.

This script helps identify why no leaf nodes or paths are found by:
1. Checking the ontology structure
2. Showing sample term IDs
3. Analyzing the graph connectivity
4. Checking RLIPP score term format

Usage:
    python debug_lineage.py -onto ontology.txt -rlipp rlipp_scores.tsv
"""

import argparse
import pandas as pd
import networkx as nx


def load_ontology(ontology_file):
    """Load ontology file and create directed graph."""
    dG = nx.DiGraph()

    with open(ontology_file) as f:
        for line in f:
            parts = line.rstrip().split()
            if len(parts) >= 3 and parts[2] == 'default':
                # This is a term-term relationship
                parent, child = parts[0], parts[1]
                dG.add_edge(parent, child)

    return dG


def main():
    parser = argparse.ArgumentParser(description='Debug lineage analysis')
    parser.add_argument('-onto', help='Ontology file', type=str, required=True)
    parser.add_argument('-rlipp', help='RLIPP scores file', type=str, required=True)
    args = parser.parse_args()

    print("="*70)
    print("LINEAGE ANALYSIS DEBUG")
    print("="*70)

    # Load ontology
    print(f"\n1. Loading ontology from: {args.onto}")
    graph = load_ontology(args.onto)
    print(f"   Loaded {len(graph.nodes())} terms and {len(graph.edges())} edges")

    # Check for roots
    print("\n2. Checking for root nodes (in_degree = 0):")
    roots = [n for n in graph.nodes() if graph.in_degree(n) == 0]
    print(f"   Found {len(roots)} root node(s)")
    if roots:
        for root in roots[:5]:  # Show first 5
            print(f"      - {root}")

    # Check for leaves
    print("\n3. Checking for leaf nodes (out_degree = 0):")
    leaves = [n for n in graph.nodes() if graph.out_degree(n) == 0]
    print(f"   Found {len(leaves)} leaf node(s)")
    if leaves:
        print("   First 10 leaves:")
        for leaf in leaves[:10]:
            print(f"      - {leaf}")
    else:
        print("   WARNING: No leaf nodes found!")
        print("   This means all nodes have children, which is unusual.")
        print("   Check if your ontology file has the correct format:")
        print("      parent child default")

    # Check intermediate nodes
    print("\n4. Checking intermediate nodes (0 < out_degree):")
    intermediate = [n for n in graph.nodes() if 0 < graph.out_degree(n)]
    print(f"   Found {len(intermediate)} intermediate node(s)")
    if intermediate:
        print("   Sample intermediate nodes (with their out_degree):")
        for node in intermediate[:10]:
            print(f"      - {node} (out_degree={graph.out_degree(node)})")

    # Sample some terms to show format
    print("\n5. Sample term IDs from ontology:")
    sample_terms = list(graph.nodes())[:10]
    for term in sample_terms:
        print(f"      - '{term}' (type: {type(term).__name__})")

    # Check connectivity
    print("\n6. Checking graph connectivity:")
    if roots:
        root = roots[0]
        # Find nodes reachable from root
        reachable = nx.descendants(graph, root)
        reachable.add(root)
        print(f"   Nodes reachable from root '{root}': {len(reachable)}")
        print(f"   Total nodes in graph: {len(graph.nodes())}")
        if len(reachable) != len(graph.nodes()):
            print(f"   WARNING: {len(graph.nodes()) - len(reachable)} nodes are NOT reachable from root!")
            unreachable = set(graph.nodes()) - reachable
            print("   First 10 unreachable nodes:")
            for node in list(unreachable)[:10]:
                print(f"      - {node}")

    # Load RLIPP scores
    print(f"\n7. Loading RLIPP scores from: {args.rlipp}")
    rlipp_df = pd.read_csv(args.rlipp, sep='\t')
    print(f"   Loaded {len(rlipp_df)} RLIPP scores")

    # Show sample RLIPP terms
    print("\n8. Sample term IDs from RLIPP file:")
    for term in rlipp_df['term'].head(10):
        print(f"      - '{term}' (type: {type(term).__name__})")

    # Check overlap
    print("\n9. Checking term ID overlap:")
    ontology_terms = set(graph.nodes())
    rlipp_terms = set(rlipp_df['term'])
    overlap = ontology_terms & rlipp_terms
    print(f"   Terms in ontology: {len(ontology_terms)}")
    print(f"   Terms in RLIPP: {len(rlipp_terms)}")
    print(f"   Overlap: {len(overlap)} terms")

    if len(overlap) == 0:
        print("\n   ERROR: NO OVERLAP between ontology and RLIPP terms!")
        print("   This means the term IDs don't match.")
        print("\n   Ontology sample terms:")
        for term in list(ontology_terms)[:5]:
            print(f"      - '{term}'")
        print("\n   RLIPP sample terms:")
        for term in list(rlipp_terms)[:5]:
            print(f"      - '{term}'")
    elif len(overlap) < len(ontology_terms):
        print(f"\n   WARNING: Only {len(overlap)}/{len(ontology_terms)} ontology terms have RLIPP scores")
        missing = ontology_terms - rlipp_terms
        print(f"   {len(missing)} ontology terms missing from RLIPP:")
        for term in list(missing)[:10]:
            print(f"      - {term}")

    # Check if any leaves have RLIPP scores
    if leaves:
        print("\n10. Checking if leaf nodes have RLIPP scores:")
        leaves_with_rlipp = [leaf for leaf in leaves if leaf in rlipp_terms]
        print(f"    Leaves with RLIPP scores: {len(leaves_with_rlipp)}/{len(leaves)}")
        if leaves_with_rlipp:
            print("    Sample leaves with RLIPP:")
            for leaf in leaves_with_rlipp[:5]:
                rlipp_val = rlipp_df[rlipp_df['term'] == leaf]['rlipp'].values[0]
                print(f"       - {leaf}: RLIPP = {rlipp_val:.3f}")

    # Test path finding
    if roots and leaves:
        print("\n11. Testing path finding from sample leaf to root:")
        root = roots[0]
        sample_leaf = leaves[0]
        try:
            paths = list(nx.all_simple_paths(graph, sample_leaf, root))
            print(f"    From '{sample_leaf}' to '{root}':")
            print(f"    Found {len(paths)} path(s)")
            if paths:
                print(f"    First path ({len(paths[0])} terms):")
                for term in paths[0][:5]:
                    print(f"       - {term}")
                if len(paths[0]) > 5:
                    print(f"       ... ({len(paths[0]) - 5} more terms)")
        except nx.NetworkXNoPath:
            print(f"    ERROR: No path found from '{sample_leaf}' to '{root}'!")
            print("    This suggests the graph is not properly connected.")

    print("\n" + "="*70)
    print("DEBUG ANALYSIS COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()
