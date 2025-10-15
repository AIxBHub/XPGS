#!/usr/bin/env python3
"""
Create Network-Style Visualization of RLIPP Lineages

This script creates hierarchical network diagrams similar to biological pathway
visualizations, showing how terms are connected in high-RLIPP lineages.

The visualization shows:
- Nodes sized by RLIPP magnitude
- Nodes colored by RLIPP sign (positive/negative)
- Hierarchical layout from root to leaves
- Direct parent-child relationships
- Node labels with term names (if OBO provided)

Usage:
    python visualize_lineage_network.py -rlipp scores.tsv -onto ontology.txt -output network.png
"""

import argparse
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

try:
    from obo_parser import create_term_name_dict, get_term_name
    OBO_SUPPORT = True
except ImportError:
    OBO_SUPPORT = False


def load_ontology(ontology_file):
    """Load ontology and create directed graph."""
    dG = nx.DiGraph()
    with open(ontology_file) as f:
        for line in f:
            parts = line.rstrip().split()
            if len(parts) >= 3 and parts[2] == 'default':
                parent, child = parts[0], parts[1]
                dG.add_edge(parent, child)
    return dG


def extract_lineage_subgraph(graph, lineage_terms):
    """Extract subgraph containing only the lineage terms."""
    return graph.subgraph(lineage_terms).copy()


def create_hierarchical_layout(graph, root):
    """
    Create hierarchical layout with root at top.

    Returns:
        dict: Maps node to (x, y) position
    """
    # Use graphviz layout if available, otherwise hierarchical
    try:
        pos = nx.nx_agraph.graphviz_layout(graph, prog='dot')
    except:
        # Fallback: manual hierarchical layout
        pos = {}

        # Assign levels (distance from root)
        levels = {}
        for node in graph.nodes():
            try:
                levels[node] = nx.shortest_path_length(graph, root, node)
            except:
                levels[node] = 0

        # Group nodes by level
        level_nodes = {}
        for node, level in levels.items():
            if level not in level_nodes:
                level_nodes[level] = []
            level_nodes[level].append(node)

        # Assign positions
        max_level = max(levels.values()) if levels else 0
        for level, nodes in level_nodes.items():
            y = max_level - level  # Root at top
            for i, node in enumerate(nodes):
                x = (i - len(nodes)/2) * 2  # Spread horizontally
                pos[node] = (x, y)

    return pos


def plot_lineage_network(lineage_path, rlipp_dict, graph, output_file,
                         name_dict=None, direction='positive'):
    """
    Create network visualization of a lineage.

    Args:
        lineage_path (list): List of terms from leaf to root
        rlipp_dict (dict): Term to RLIPP score mapping
        graph (nx.DiGraph): Full ontology graph
        output_file (str): Output file path
        name_dict (dict): Optional GO ID to name mapping
        direction (str): 'positive' or 'negative'
    """
    # Extract subgraph for this lineage
    subgraph = extract_lineage_subgraph(graph, lineage_path)

    # Get root (first or last depending on path direction)
    root = lineage_path[-1]  # Assuming path goes leaf->root

    # Create layout
    pos = create_hierarchical_layout(subgraph, root)

    # Prepare node attributes
    node_colors = []
    node_sizes = []
    node_labels = {}

    for node in subgraph.nodes():
        rlipp = rlipp_dict.get(node, 0)

        # Color by RLIPP value (blue for positive, red for negative)
        node_colors.append(rlipp)

        # Size by absolute RLIPP (larger = more important)
        size = 300 + abs(rlipp) * 2000
        node_sizes.append(size)

        # Label with name or ID
        if name_dict:
            label = get_term_name(node, name_dict, max_length=20)
        else:
            label = node[:15] if len(node) > 15 else node
        node_labels[node] = label

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(14, 10))

    # Define colormap (red -> white -> blue)
    if direction == 'positive':
        cmap = plt.cm.Blues
        vmin, vmax = 0, max(node_colors) if node_colors else 1
    else:
        cmap = plt.cm.Reds
        vmin, vmax = min(node_colors) if node_colors else -1, 0

    # Draw edges
    nx.draw_networkx_edges(
        subgraph, pos,
        edge_color='gray',
        arrows=True,
        arrowsize=20,
        arrowstyle='->',
        width=2,
        alpha=0.6,
        ax=ax
    )

    # Draw nodes
    nodes = nx.draw_networkx_nodes(
        subgraph, pos,
        node_color=node_colors,
        node_size=node_sizes,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        alpha=0.9,
        edgecolors='black',
        linewidths=2,
        ax=ax
    )

    # Draw labels
    nx.draw_networkx_labels(
        subgraph, pos,
        labels=node_labels,
        font_size=9,
        font_weight='bold',
        ax=ax
    )

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('RLIPP Score', fontsize=12, weight='bold')

    # Title
    leaf_name = node_labels[lineage_path[0]] if lineage_path else "Unknown"
    ax.set_title(
        f'{direction.capitalize()} RLIPP Lineage Network\nLeaf: {leaf_name}',
        fontsize=14,
        fontweight='bold',
        pad=20
    )

    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved network visualization: {output_file}")


def plot_multiple_lineages_combined(top_lineages, rlipp_dict, graph, output_file,
                                    name_dict=None, direction='positive', max_lineages=5):
    """
    Create combined network showing multiple top lineages.

    Args:
        top_lineages (list): List of lineage dictionaries
        rlipp_dict (dict): Term to RLIPP mapping
        graph (nx.DiGraph): Full ontology graph
        output_file (str): Output file path
        name_dict (dict): Optional term name mapping
        direction (str): 'positive' or 'negative'
        max_lineages (int): Maximum number of lineages to include
    """
    # Collect all terms from top lineages
    all_terms = set()
    for lineage in top_lineages[:max_lineages]:
        all_terms.update(lineage['path'])

    # Extract subgraph
    subgraph = graph.subgraph(all_terms).copy()

    # Find root
    roots = [n for n in subgraph.nodes() if subgraph.in_degree(n) == 0]
    root = roots[0] if roots else list(subgraph.nodes())[0]

    # Create layout
    pos = create_hierarchical_layout(subgraph, root)

    # Prepare node attributes
    node_colors = []
    node_sizes = []
    node_labels = {}

    for node in subgraph.nodes():
        rlipp = rlipp_dict.get(node, 0)
        node_colors.append(rlipp)
        size = 300 + abs(rlipp) * 1500
        node_sizes.append(size)

        if name_dict:
            label = get_term_name(node, name_dict, max_length=15)
        else:
            label = node[:12] if len(node) > 12 else node
        node_labels[node] = label

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(16, 12))

    # Colormap
    if direction == 'positive':
        cmap = plt.cm.Blues
        vmin, vmax = 0, max(node_colors) if node_colors else 1
    else:
        cmap = plt.cm.Reds
        vmin, vmax = min(node_colors) if node_colors else -1, 0

    # Draw edges (color by whether it's in top lineage)
    edge_colors = []
    edge_widths = []
    for edge in subgraph.edges():
        # Check if this edge is in the top lineage
        in_top = any(
            edge[0] in lineage['path'] and edge[1] in lineage['path'] and
            abs(lineage['path'].index(edge[0]) - lineage['path'].index(edge[1])) == 1
            for lineage in top_lineages[:1]  # Highlight top lineage
        )
        edge_colors.append('darkblue' if in_top else 'lightgray')
        edge_widths.append(3 if in_top else 1.5)

    nx.draw_networkx_edges(
        subgraph, pos,
        edge_color=edge_colors,
        width=edge_widths,
        arrows=True,
        arrowsize=15,
        arrowstyle='->',
        alpha=0.7,
        ax=ax
    )

    # Draw nodes
    nodes = nx.draw_networkx_nodes(
        subgraph, pos,
        node_color=node_colors,
        node_size=node_sizes,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        alpha=0.9,
        edgecolors='black',
        linewidths=2,
        ax=ax
    )

    # Draw labels
    nx.draw_networkx_labels(
        subgraph, pos,
        labels=node_labels,
        font_size=8,
        font_weight='bold',
        ax=ax
    )

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('RLIPP Score', fontsize=12, weight='bold')

    # Title
    ax.set_title(
        f'Combined Network: Top {max_lineages} {direction.capitalize()} RLIPP Lineages\n'
        f'(Dark edges highlight top lineage)',
        fontsize=14,
        fontweight='bold',
        pad=20
    )

    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Saved combined network: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Create network-style lineage visualizations',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-lineages', help='Lineage TSV file from analyze_rlipp_lineages.py',
                       type=str, required=True)
    parser.add_argument('-rlipp', help='RLIPP scores TSV file',
                       type=str, required=True)
    parser.add_argument('-onto', help='Ontology file',
                       type=str, required=True)
    parser.add_argument('-output', help='Output prefix',
                       type=str, default='network')
    parser.add_argument('-direction', help='Direction (positive/negative)',
                       type=str, default='positive', choices=['positive', 'negative'])
    parser.add_argument('-top_n', help='Number of individual networks to create',
                       type=int, default=3)
    parser.add_argument('-combined', help='Number of lineages in combined network',
                       type=int, default=5)
    parser.add_argument('-obo', help='GO OBO file for term names',
                       type=str, default=None)

    args = parser.parse_args()

    print("="*70)
    print("NETWORK-STYLE LINEAGE VISUALIZATION")
    print("="*70)

    # Load data
    print(f"\nLoading RLIPP scores from: {args.rlipp}")
    rlipp_df = pd.read_csv(args.rlipp, sep='\t')
    rlipp_dict = dict(zip(rlipp_df['term'], rlipp_df['rlipp']))
    print(f"Loaded {len(rlipp_dict)} terms")

    print(f"\nLoading lineages from: {args.lineages}")
    lineages_df = pd.read_csv(args.lineages, sep='\t')
    print(f"Loaded {len(lineages_df)} lineages")

    print(f"\nLoading ontology from: {args.onto}")
    graph = load_ontology(args.onto)
    print(f"Loaded {len(graph.nodes())} terms, {len(graph.edges())} edges")

    # Load term names if provided
    name_dict = None
    if args.obo and OBO_SUPPORT:
        print(f"\nLoading term names from: {args.obo}")
        try:
            name_dict = create_term_name_dict(args.obo)
            print(f"Loaded {len(name_dict)} term names")
        except Exception as e:
            print(f"Warning: Could not load OBO file: {e}")

    # Convert lineage paths from string to list
    print("\nParsing lineage paths...")
    lineages = []
    for _, row in lineages_df.iterrows():
        path = row['path'].split(' -> ')
        rlipp_scores = [float(x) for x in row['rlipp_scores'].split(', ')]
        lineages.append({
            'path': path,
            'rlipp_scores': rlipp_scores,
            'leaf': row['leaf'],
            'combined_score': row['combined_score'],
            'consistency': row['consistency'],
            'mean_score': row['mean_score']
        })

    # Create individual network visualizations
    print(f"\nCreating {args.top_n} individual network visualizations...")
    for i, lineage in enumerate(lineages[:args.top_n], 1):
        output_file = f"{args.output}_{args.direction}_network_{i}.png"
        plot_lineage_network(
            lineage['path'],
            rlipp_dict,
            graph,
            output_file,
            name_dict,
            args.direction
        )

    # Create combined network
    print(f"\nCreating combined network with {args.combined} lineages...")
    combined_file = f"{args.output}_{args.direction}_network_combined.png"
    plot_multiple_lineages_combined(
        lineages,
        rlipp_dict,
        graph,
        combined_file,
        name_dict,
        args.direction,
        args.combined
    )

    print("\n" + "="*70)
    print("NETWORK VISUALIZATION COMPLETE")
    print("="*70)
    print(f"\nCreated {args.top_n + 1} network visualizations")


if __name__ == "__main__":
    main()
