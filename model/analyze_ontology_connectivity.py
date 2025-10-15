#!/usr/bin/env python3
"""
Analyze ontology connectivity and structure.

This script provides detailed analysis of the ontology graph to identify
connectivity issues, disconnected components, and structural problems.

Usage:
    python analyze_ontology_connectivity.py -onto ontology.txt
"""

import argparse
import networkx as nx
from collections import defaultdict


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


def analyze_connectivity(graph):
    """Analyze graph connectivity in detail."""

    print("="*70)
    print("ONTOLOGY CONNECTIVITY ANALYSIS")
    print("="*70)

    # Basic stats
    print(f"\n1. BASIC STATISTICS:")
    print(f"   Total nodes: {len(graph.nodes())}")
    print(f"   Total edges: {len(graph.edges())}")

    # Find roots (nodes with no incoming edges)
    roots = [n for n in graph.nodes() if graph.in_degree(n) == 0]
    print(f"\n2. ROOT NODES (in_degree = 0):")
    print(f"   Found {len(roots)} root(s)")
    for i, root in enumerate(roots[:10], 1):
        descendants = len(nx.descendants(graph, root))
        print(f"   {i}. '{root}' - can reach {descendants} descendants")

    # Find leaves (nodes with no outgoing edges)
    leaves = [n for n in graph.nodes() if graph.out_degree(n) == 0]
    print(f"\n3. LEAF NODES (out_degree = 0):")
    print(f"   Found {len(leaves)} leaf node(s)")
    if leaves:
        print(f"   Sample leaves:")
        for leaf in leaves[:10]:
            print(f"      - '{leaf}'")

    # Check if it's a DAG
    print(f"\n4. GRAPH TYPE:")
    is_dag = nx.is_directed_acyclic_graph(graph)
    print(f"   Is Directed Acyclic Graph (DAG): {is_dag}")
    if not is_dag:
        print("   WARNING: Graph has cycles! This will cause problems.")
        try:
            cycles = list(nx.simple_cycles(graph))
            print(f"   Found {len(cycles)} cycle(s)")
            if cycles:
                print("   First cycle:")
                print(f"      {' -> '.join(cycles[0])}")
        except:
            pass

    # Analyze weakly connected components
    print(f"\n5. CONNECTED COMPONENTS:")
    undirected = graph.to_undirected()
    components = list(nx.connected_components(undirected))
    print(f"   Found {len(components)} weakly connected component(s)")

    if len(components) > 1:
        print("\n   WARNING: Graph has multiple disconnected components!")
        print("   Component sizes:")
        for i, comp in enumerate(sorted(components, key=len, reverse=True), 1):
            print(f"   {i}. {len(comp)} nodes")
            if len(comp) <= 5:
                print(f"      Nodes: {list(comp)}")
            else:
                print(f"      Sample nodes: {list(comp)[:5]}")

    # For each root, see what it can reach
    print(f"\n6. REACHABILITY FROM ROOTS:")
    for i, root in enumerate(roots[:5], 1):
        reachable = set(nx.descendants(graph, root))
        reachable.add(root)
        print(f"   Root {i}: '{root}'")
        print(f"      Can reach: {len(reachable)} / {len(graph.nodes())} nodes ({100*len(reachable)/len(graph.nodes()):.1f}%)")

        # Find which leaves this root can reach
        reachable_leaves = [leaf for leaf in leaves if leaf in reachable]
        print(f"      Can reach: {len(reachable_leaves)} / {len(leaves)} leaves")

    # Find orphaned nodes (not reachable from any root)
    print(f"\n7. ORPHANED NODES:")
    all_reachable = set()
    for root in roots:
        reachable = set(nx.descendants(graph, root))
        reachable.add(root)
        all_reachable.update(reachable)

    orphaned = set(graph.nodes()) - all_reachable
    print(f"   Found {len(orphaned)} orphaned node(s) (not reachable from any root)")
    if orphaned:
        print("   Sample orphaned nodes:")
        for node in list(orphaned)[:10]:
            in_deg = graph.in_degree(node)
            out_deg = graph.out_degree(node)
            print(f"      - '{node}' (in={in_deg}, out={out_deg})")

    # Analyze node degrees
    print(f"\n8. NODE DEGREE DISTRIBUTION:")
    in_degrees = [graph.in_degree(n) for n in graph.nodes()]
    out_degrees = [graph.out_degree(n) for n in graph.nodes()]

    print(f"   In-degree:  min={min(in_degrees)}, max={max(in_degrees)}, avg={sum(in_degrees)/len(in_degrees):.1f}")
    print(f"   Out-degree: min={min(out_degrees)}, max={max(out_degrees)}, avg={sum(out_degrees)/len(out_degrees):.1f}")

    # Find nodes with very high out-degree
    high_out = [(n, graph.out_degree(n)) for n in graph.nodes() if graph.out_degree(n) > 50]
    if high_out:
        high_out.sort(key=lambda x: x[1], reverse=True)
        print(f"\n   Nodes with >50 children:")
        for node, deg in high_out[:5]:
            print(f"      - '{node}': {deg} children")

    # Test path finding
    print(f"\n9. PATH FINDING TEST:")
    if roots and leaves:
        root = roots[0]
        # Try first 5 leaves
        success_count = 0
        fail_count = 0

        for leaf in leaves[:min(5, len(leaves))]:
            try:
                paths = list(nx.all_simple_paths(graph, leaf, root))
                if paths:
                    success_count += 1
                    print(f"   ✓ '{leaf}' → '{root}': {len(paths)} path(s), length {len(paths[0])}")
                else:
                    fail_count += 1
                    print(f"   ✗ '{leaf}' → '{root}': 0 paths")
            except nx.NetworkXNoPath:
                fail_count += 1
                print(f"   ✗ '{leaf}' → '{root}': No path (disconnected)")

        print(f"\n   Summary: {success_count} successful, {fail_count} failed")

        if fail_count > 0:
            print("\n   DIAGNOSIS: Some leaves cannot reach the root!")
            print("   Possible causes:")
            print("   - Multiple disconnected components")
            print("   - Wrong root identified")
            print("   - Leaf is in a different subgraph")

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)

    return {
        'roots': roots,
        'leaves': leaves,
        'components': components,
        'orphaned': orphaned,
        'is_dag': is_dag
    }


def suggest_fixes(graph, analysis):
    """Suggest fixes based on analysis."""

    print("\n" + "="*70)
    print("SUGGESTED FIXES")
    print("="*70)

    roots = analysis['roots']
    leaves = analysis['leaves']
    components = analysis['components']
    orphaned = analysis['orphaned']

    if len(roots) == 0:
        print("\n✗ PROBLEM: No root nodes found!")
        print("  FIX: Check if edges are in the correct direction")
        print("       Edges should go FROM parent TO child")

    elif len(roots) > 1:
        print(f"\n⚠ WARNING: Multiple roots found ({len(roots)})")
        print("  FIX OPTIONS:")
        print("  1. Add a super-root that connects all current roots")
        print("  2. Choose the root with most descendants as THE root")
        print("  3. Check if some roots are actually disconnected components")

        # Find which root has most descendants
        root_sizes = [(r, len(nx.descendants(graph, r))) for r in roots]
        root_sizes.sort(key=lambda x: x[1], reverse=True)
        print(f"\n  Largest root: '{root_sizes[0][0]}' with {root_sizes[0][1]} descendants")

    if len(leaves) == 0:
        print("\n✗ PROBLEM: No leaf nodes found!")
        print("  FIX: Check if all nodes have outgoing edges")
        print("       This is unusual for a biological ontology")

    if len(components) > 1:
        print(f"\n✗ PROBLEM: Graph has {len(components)} disconnected components!")
        print("  FIX: The ontology should be a single connected component")
        print("       Check for:")
        print("       - Missing parent-child relationships")
        print("       - Incorrect term IDs")
        print("       - Multiple separate ontologies in one file")

    if len(orphaned) > 0:
        print(f"\n✗ PROBLEM: {len(orphaned)} orphaned nodes (not reachable from root)")
        print("  FIX: These nodes need to be connected to the main hierarchy")

    if not analysis['is_dag']:
        print("\n✗ PROBLEM: Graph has cycles!")
        print("  FIX: Remove circular relationships")
        print("       Biological ontologies should be DAGs")

    # Check if the main issue is identified root vs actual structure
    if len(roots) > 0 and len(leaves) > 0:
        root = roots[0]
        reachable_leaves = []
        for leaf in leaves[:100]:  # Test first 100 leaves
            try:
                if nx.has_path(graph, leaf, root):
                    reachable_leaves.append(leaf)
            except:
                pass

        if len(reachable_leaves) == 0:
            print("\n✗ CRITICAL: NO leaves can reach the identified root!")
            print("  DIAGNOSIS: The root might be wrong, or edges are reversed")
            print("  FIX OPTIONS:")
            print("  1. Check if edges should be reversed (child -> parent instead of parent -> child)")
            print("  2. Try using a different node as root")
            print("  3. Verify the ontology file format")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze ontology connectivity',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-onto', help='Ontology file', type=str, required=True)
    args = parser.parse_args()

    print(f"Loading ontology from: {args.onto}\n")
    graph = load_ontology(args.onto)

    analysis = analyze_connectivity(graph)
    suggest_fixes(graph, analysis)


if __name__ == "__main__":
    main()
