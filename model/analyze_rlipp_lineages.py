#!/usr/bin/env python3
"""
Identify Lineages with Consistently High or Low RLIPP Scores

This script analyzes the ontology hierarchy to find lineages (paths from leaves to root)
that have consistently high positive or negative RLIPP scores. This helps identify:
- Important biological pathways (high positive RLIPP lineages)
- Redundant pathways (high negative RLIPP lineages)

A lineage is a path from a leaf term through intermediate terms to the root term.
We score each lineage based on the consistency and magnitude of RLIPP scores along the path.

Usage:
    python analyze_rlipp_lineages.py -rlipp rlipp_scores.tsv -onto ontology.txt -output lineages.tsv

Author: [Your name]
Date: 2024
"""

import argparse
import pandas as pd
import networkx as nx
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
try:
    from obo_parser import create_term_name_dict, format_term_label, get_term_name
    OBO_SUPPORT = True
except ImportError:
    OBO_SUPPORT = False
    print("Warning: obo_parser not found. Using GO IDs instead of names.")


def load_ontology(ontology_file):
    """
    Load ontology file and create directed graph.

    Args:
        ontology_file (str): Path to ontology file

    Returns:
        networkx.DiGraph: Directed graph of ontology
    """
    dG = nx.DiGraph()

    with open(ontology_file) as f:
        for line in f:
            parts = line.rstrip().split()
            if len(parts) >= 3 and parts[2] == 'default':
                # This is a term-term relationship
                parent, child = parts[0], parts[1]
                dG.add_edge(parent, child)

    return dG


def find_all_paths_to_root(graph, root):
    """
    Find all paths from each leaf node to the root.

    Args:
        graph (networkx.DiGraph): Ontology graph
        root (str): Root term

    Returns:
        dict: Maps leaf terms to list of paths to root
              Each path is a list of terms from leaf to root
    """
    # Find all leaf nodes (nodes with no children)
    leaves = [n for n in graph.nodes() if graph.out_degree(n) == 0]

    all_paths = {}

    for leaf in leaves:
        # Find all simple paths from leaf to root
        try:
            paths = list(nx.all_simple_paths(graph, leaf, root))
            if paths:
                all_paths[leaf] = paths
        except nx.NetworkXNoPath:
            # No path from this leaf to root (disconnected component)
            continue

    return all_paths


def score_lineage(path, rlipp_dict, direction='positive'):
    """
    Score a lineage based on RLIPP values along the path.

    Scoring considers:
    1. Mean absolute RLIPP score
    2. Consistency (what fraction of terms have the expected sign)
    3. Magnitude (sum of RLIPP scores in expected direction)

    Args:
        path (list): List of terms from leaf to root
        rlipp_dict (dict): Maps term to RLIPP score
        direction (str): 'positive' or 'negative'

    Returns:
        dict: Scoring metrics for this lineage
    """
    rlipp_scores = []

    for term in path:
        if term in rlipp_dict:
            rlipp_scores.append(rlipp_dict[term])

    if not rlipp_scores:
        return None

    rlipp_array = np.array(rlipp_scores)

    # Calculate metrics
    if direction == 'positive':
        # For positive lineages, we want high positive RLIPP scores
        consistency = (rlipp_array > 0).sum() / len(rlipp_array)
        magnitude = rlipp_array[rlipp_array > 0].sum() if (rlipp_array > 0).any() else 0
        mean_score = rlipp_array[rlipp_array > 0].mean() if (rlipp_array > 0).any() else 0
    else:
        # For negative lineages, we want high negative RLIPP scores
        consistency = (rlipp_array < 0).sum() / len(rlipp_array)
        magnitude = abs(rlipp_array[rlipp_array < 0].sum()) if (rlipp_array < 0).any() else 0
        mean_score = abs(rlipp_array[rlipp_array < 0].mean()) if (rlipp_array < 0).any() else 0

    # Combined score: consistency * magnitude * mean_score
    # This favors lineages with many terms in the expected direction with high scores
    combined_score = consistency * magnitude * mean_score

    return {
        'path': path,
        'rlipp_scores': rlipp_scores,
        'length': len(rlipp_scores),
        'consistency': consistency,
        'magnitude': magnitude,
        'mean_score': mean_score,
        'combined_score': combined_score,
        'median_score': np.median(np.abs(rlipp_array)),
        'min_score': rlipp_array.min(),
        'max_score': rlipp_array.max()
    }


def identify_top_lineages(all_paths, rlipp_dict, direction='positive', top_n=20,
                          min_consistency=0.5, min_length=2):
    """
    Identify top lineages with consistently high/low RLIPP scores.

    Args:
        all_paths (dict): Maps leaf to paths to root
        rlipp_dict (dict): Maps term to RLIPP score
        direction (str): 'positive' or 'negative'
        top_n (int): Number of top lineages to return
        min_consistency (float): Minimum fraction of terms with expected sign
        min_length (int): Minimum path length

    Returns:
        list: Top lineages sorted by combined score
    """
    lineage_scores = []

    for leaf, paths in all_paths.items():
        for path in paths:
            score_dict = score_lineage(path, rlipp_dict, direction)

            if score_dict is None:
                continue

            # Filter by consistency and length
            if (score_dict['consistency'] >= min_consistency and
                score_dict['length'] >= min_length):

                score_dict['leaf'] = leaf
                score_dict['root'] = path[-1] if path else None
                lineage_scores.append(score_dict)

    # Sort by combined score
    lineage_scores.sort(key=lambda x: x['combined_score'], reverse=True)

    return lineage_scores[:top_n]


def visualize_lineage(lineage, output_file, direction='positive', name_dict=None):
    """
    Create a visualization of a single lineage with RLIPP scores.

    Args:
        lineage (dict): Lineage score dictionary
        output_file (str): Output file path
        direction (str): 'positive' or 'negative'
        name_dict (dict, optional): Maps GO IDs to names
    """
    path = lineage['path']
    scores = lineage['rlipp_scores']

    fig, ax = plt.subplots(1, 1, figsize=(12, max(6, len(path) * 0.4)))

    # Create horizontal bar plot
    y_pos = np.arange(len(path))
    colors = ['green' if s > 0 else 'red' for s in scores]

    ax.barh(y_pos, scores, color=colors, alpha=0.7, edgecolor='black')
    ax.set_yticks(y_pos)

    # Use term names if available
    if name_dict:
        labels = [get_term_name(term, name_dict, max_length=50) for term in path]
    else:
        labels = path
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('RLIPP Score', fontsize=12)
    ax.set_title(
        f'{direction.capitalize()} Lineage\n'
        f'Consistency: {lineage["consistency"]:.2f} | '
        f'Mean: {lineage["mean_score"]:.3f} | '
        f'Combined Score: {lineage["combined_score"]:.3f}',
        fontsize=14, fontweight='bold'
    )
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    ax.grid(True, alpha=0.3, axis='x')
    ax.invert_yaxis()  # Leaf at top, root at bottom

    # Add annotations
    for i, (term, score) in enumerate(zip(path, scores)):
        ax.text(
            score + (0.01 if score > 0 else -0.01),
            i,
            f'{score:.3f}',
            va='center',
            ha='left' if score > 0 else 'right',
            fontsize=8
        )

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved lineage visualization: {output_file}")


def create_lineage_summary_plot(top_lineages, output_file, direction='positive', name_dict=None):
    """
    Create summary visualization of top lineages.

    Args:
        top_lineages (list): List of top lineages
        output_file (str): Output file path
        direction (str): 'positive' or 'negative'
        name_dict (dict, optional): Maps GO IDs to names
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Extract data
    lengths = [l['length'] for l in top_lineages]
    consistencies = [l['consistency'] for l in top_lineages]
    mean_scores = [l['mean_score'] for l in top_lineages]
    combined_scores = [l['combined_score'] for l in top_lineages]

    # Plot 1: Length distribution
    axes[0, 0].hist(lengths, bins=20, edgecolor='black', alpha=0.7)
    axes[0, 0].set_xlabel('Lineage Length (number of terms)')
    axes[0, 0].set_ylabel('Count')
    axes[0, 0].set_title('Distribution of Lineage Lengths')
    axes[0, 0].grid(True, alpha=0.3)

    # Plot 2: Consistency distribution
    axes[0, 1].hist(consistencies, bins=20, edgecolor='black', alpha=0.7, color='orange')
    axes[0, 1].set_xlabel('Consistency (fraction of terms with expected sign)')
    axes[0, 1].set_ylabel('Count')
    axes[0, 1].set_title('Distribution of Consistency Scores')
    axes[0, 1].grid(True, alpha=0.3)

    # Plot 3: Mean score vs consistency
    scatter = axes[1, 0].scatter(consistencies, mean_scores, c=combined_scores,
                                  cmap='viridis', s=100, alpha=0.6, edgecolors='black')
    axes[1, 0].set_xlabel('Consistency')
    axes[1, 0].set_ylabel('Mean RLIPP Score')
    axes[1, 0].set_title('Mean Score vs Consistency')
    axes[1, 0].grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=axes[1, 0], label='Combined Score')

    # Plot 4: Top lineages by combined score
    top_10 = sorted(top_lineages, key=lambda x: x['combined_score'], reverse=True)[:10]

    # Use term names if available
    if name_dict:
        labels = [get_term_name(l['leaf'], name_dict, max_length=30) for l in top_10]
    else:
        labels = [f"{l['leaf'][:20]}..." if len(l['leaf']) > 20 else l['leaf'] for l in top_10]
    scores = [l['combined_score'] for l in top_10]

    y_pos = np.arange(len(labels))
    axes[1, 1].barh(y_pos, scores, alpha=0.7, edgecolor='black')
    axes[1, 1].set_yticks(y_pos)
    axes[1, 1].set_yticklabels(labels, fontsize=8)
    axes[1, 1].set_xlabel('Combined Score')
    axes[1, 1].set_title('Top 10 Lineages by Combined Score')
    axes[1, 1].invert_yaxis()
    axes[1, 1].grid(True, alpha=0.3, axis='x')

    plt.suptitle(f'{direction.capitalize()} RLIPP Lineages - Summary',
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved summary plot: {output_file}")


def save_lineages_to_file(lineages, output_file):
    """
    Save lineage analysis results to TSV file.

    Args:
        lineages (list): List of lineage dictionaries
        output_file (str): Output file path
    """
    rows = []

    for i, lineage in enumerate(lineages, 1):
        path_str = ' -> '.join(lineage['path'])
        scores_str = ', '.join([f'{s:.3f}' for s in lineage['rlipp_scores']])

        rows.append({
            'rank': i,
            'leaf': lineage['leaf'],
            'root': lineage['root'],
            'length': lineage['length'],
            'consistency': lineage['consistency'],
            'mean_score': lineage['mean_score'],
            'median_score': lineage['median_score'],
            'magnitude': lineage['magnitude'],
            'combined_score': lineage['combined_score'],
            'min_score': lineage['min_score'],
            'max_score': lineage['max_score'],
            'path': path_str,
            'rlipp_scores': scores_str
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep='\t', index=False)

    print(f"Saved lineage analysis: {output_file}")
    return df


def main():
    parser = argparse.ArgumentParser(
        description='Identify lineages with consistently high or low RLIPP scores',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-rlipp', help='RLIPP scores TSV file from calculate_rlipp.py',
                       type=str, required=True)
    parser.add_argument('-onto', help='Ontology file used in VNN training',
                       type=str, required=True)
    parser.add_argument('-output', help='Output prefix for results files',
                       type=str, default='lineages')
    parser.add_argument('-direction', help='Direction to analyze (positive/negative/both)',
                       type=str, default='both', choices=['positive', 'negative', 'both'])
    parser.add_argument('-top_n', help='Number of top lineages to report',
                       type=int, default=20)
    parser.add_argument('-min_consistency', help='Minimum consistency score (0-1)',
                       type=float, default=0.5)
    parser.add_argument('-min_length', help='Minimum lineage length',
                       type=int, default=2)
    parser.add_argument('-visualize_top', help='Number of top lineages to visualize individually',
                       type=int, default=5)
    parser.add_argument('-obo', help='GO OBO file for term names (optional)', type=str, default=None)

    args = parser.parse_args()

    print("="*60)
    print("RLIPP LINEAGE ANALYSIS")
    print("="*60)

    # Load RLIPP scores
    print(f"\nLoading RLIPP scores from: {args.rlipp}")
    rlipp_df = pd.read_csv(args.rlipp, sep='\t')
    rlipp_dict = dict(zip(rlipp_df['term'], rlipp_df['rlipp']))
    print(f"Loaded RLIPP scores for {len(rlipp_dict)} terms")

    # Load GO term names if OBO file provided
    name_dict = None
    if args.obo and OBO_SUPPORT:
        print(f"\nLoading GO term names from: {args.obo}")
        try:
            name_dict = create_term_name_dict(args.obo)
            print(f"Loaded {len(name_dict)} term names")
        except Exception as e:
            print(f"Warning: Could not load OBO file: {e}")
            print("Continuing with GO IDs...")
            name_dict = None
    elif args.obo and not OBO_SUPPORT:
        print("\nWarning: OBO file provided but obo_parser module not available")
        print("Continuing with GO IDs...")

    # Load ontology
    print(f"\nLoading ontology from: {args.onto}")
    graph = load_ontology(args.onto)
    print(f"Loaded ontology with {len(graph.nodes())} terms and {len(graph.edges())} edges")

    # Find root
    roots = [n for n in graph.nodes() if graph.in_degree(n) == 0]
    if len(roots) != 1:
        print(f"ERROR: Expected 1 root, found {len(roots)}")
        return
    root = roots[0]
    print(f"Root term: {root}")

    # Find all paths
    print("\nFinding all paths from leaves to root...")
    all_paths = find_all_paths_to_root(graph, root)
    total_paths = sum(len(paths) for paths in all_paths.values())
    print(f"Found {len(all_paths)} leaf nodes with {total_paths} total paths to root")

    # Analyze lineages
    directions = ['positive', 'negative'] if args.direction == 'both' else [args.direction]

    for direction in directions:
        print(f"\n{'='*60}")
        print(f"Analyzing {direction.upper()} lineages")
        print(f"{'='*60}")

        top_lineages = identify_top_lineages(
            all_paths,
            rlipp_dict,
            direction=direction,
            top_n=args.top_n,
            min_consistency=args.min_consistency,
            min_length=args.min_length
        )

        print(f"\nFound {len(top_lineages)} lineages meeting criteria:")
        print(f"  - Minimum consistency: {args.min_consistency}")
        print(f"  - Minimum length: {args.min_length}")

        if len(top_lineages) == 0:
            print(f"No lineages found for {direction} direction. Try lowering min_consistency or min_length.")
            continue

        # Save to file
        output_file = f'{args.output}_{direction}.tsv'
        df = save_lineages_to_file(top_lineages, output_file)

        # Print top 5
        print(f"\nTop 5 {direction} lineages:")
        for i, lineage in enumerate(top_lineages[:5], 1):
            print(f"\n{i}. Leaf: {lineage['leaf']}")
            print(f"   Length: {lineage['length']} terms")
            print(f"   Consistency: {lineage['consistency']:.2f}")
            print(f"   Mean score: {lineage['mean_score']:.3f}")
            print(f"   Combined score: {lineage['combined_score']:.3f}")
            print(f"   Path: {' -> '.join(lineage['path'][:3])} ... {lineage['path'][-1]}")

        # Create summary plot
        summary_file = f'{args.output}_{direction}_summary.png'
        create_lineage_summary_plot(top_lineages, summary_file, direction, name_dict)

        # Visualize top individual lineages
        print(f"\nCreating individual visualizations for top {args.visualize_top} lineages...")
        for i, lineage in enumerate(top_lineages[:args.visualize_top], 1):
            viz_file = f'{args.output}_{direction}_lineage_{i}.png'
            visualize_lineage(lineage, viz_file, direction, name_dict)

    print("\n" + "="*60)
    print("LINEAGE ANALYSIS COMPLETE")
    print("="*60)
    print(f"\nOutput files created:")
    for direction in directions:
        print(f"  - {args.output}_{direction}.tsv (lineage data)")
        print(f"  - {args.output}_{direction}_summary.png (summary plots)")
        for i in range(1, min(args.visualize_top + 1, args.top_n + 1)):
            print(f"  - {args.output}_{direction}_lineage_{i}.png (individual lineage)")
    print("="*60)


if __name__ == "__main__":
    main()
