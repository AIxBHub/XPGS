#!/usr/bin/env python3
"""
Analyze Which Samples and Genes Drive Lineage RLIPP Scores

This script identifies:
1. Which test samples contribute most to high/low RLIPP scores in a lineage
2. Which genes in those samples are most activated for the lineage terms

This helps answer:
- Which patients/samples drive the lineage results?
- What genetic variants activate these pathways?
- Which genes are most important for high-RLIPP lineages?

Usage:
    python analyze_lineage_drivers.py -model model.pt -data test.pt -lineages lineages.tsv
"""

import argparse
import torch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import sys
import os

# Import VNN modules
try:
    from data_wrapper import TrainingDataWrapper
    from dcell_nn import DCellNN
    import utils
except ImportError:
    print("Error: Could not import VNN modules. Make sure you're in the model directory.")
    sys.exit(1)


def load_model_and_data(model_path, data_wrapper):
    """Load trained model."""
    print(f"Loading model from: {model_path}")

    model = DCellNN(
        term_size_map=data_wrapper.term_size_map,
        term_direct_gene_map=data_wrapper.term_direct_gene_map,
        dG=data_wrapper.dG,
        gene_dim=len(data_wrapper.gene_id_mapping),
        num_hiddens_genotype=data_wrapper.num_hiddens_genotype,
        cuda_id=data_wrapper.cuda
    )

    checkpoint = torch.load(model_path, map_location=data_wrapper.cuda)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()

    print(f"Model loaded successfully")
    return model


def load_test_data(test_file):
    """Load test data from .pt file."""
    print(f"\nLoading test data from: {test_file}")

    data = torch.load(test_file)

    # Handle different possible formats
    if isinstance(data, dict):
        if 'x' in data and 'y' in data:
            x_test = data['x']
            y_test = data['y']
        else:
            raise ValueError("Test file dict must contain 'x' and 'y' keys")
    elif isinstance(data, (list, tuple)) and len(data) >= 2:
        x_test, y_test = data[0], data[1]
    else:
        raise ValueError("Unsupported test data format")

    print(f"Test data: {x_test.shape[0]} samples, {x_test.shape[1]} features")
    return x_test, y_test


def extract_term_embeddings_per_sample(model, x_data, device):
    """
    Extract embeddings for each term for each sample.

    Returns:
        dict: {term: tensor of shape (n_samples, hidden_dim)}
    """
    print("\nExtracting term embeddings for each sample...")

    model.eval()
    with torch.no_grad():
        x_data = x_data.to(device)

        # Forward pass through the model
        term_embeddings = {}

        # Get gene embeddings
        gene_input = x_data

        # Process through each layer
        for i, layer_terms in enumerate(model.term_layer_list):
            for term in layer_terms:
                # Get the term's computation module
                term_nn_dict = model.term_nn_dict[term]

                # Get inputs for this term
                if i == 0:
                    # First layer: directly from genes
                    term_mask = model.term_direct_gene_map[term]
                    term_input = utils.build_input_vector(gene_input, term_mask)
                else:
                    # Higher layers: from children terms
                    child_terms = model.term_neighbor_map[term]
                    child_inputs = [term_embeddings[child] for child in child_terms]
                    if child_inputs:
                        term_input = torch.cat(child_inputs, dim=1)
                    else:
                        continue

                # Forward through this term's network
                out = term_nn_dict['dropout'](term_input)
                out = term_nn_dict['linear'](out)
                out = torch.tanh(out)
                if 'batchnorm' in term_nn_dict:
                    out = term_nn_dict['batchnorm'](out)

                # Store embedding (shape: n_samples x hidden_dim)
                term_embeddings[term] = out.cpu()

    print(f"Extracted embeddings for {len(term_embeddings)} terms")
    return term_embeddings


def analyze_sample_contributions(term_embeddings, lineage_path, y_test):
    """
    Identify which samples have highest activation for lineage terms.

    Returns:
        DataFrame: Sample contributions with metrics
    """
    print("\nAnalyzing sample contributions to lineage...")

    n_samples = y_test.shape[0]

    # For each sample, compute activation metrics across lineage
    sample_scores = []

    for sample_idx in range(n_samples):
        # Compute metrics for this sample
        activations = []

        for term in lineage_path:
            if term in term_embeddings:
                embedding = term_embeddings[term][sample_idx]  # (hidden_dim,)
                # Use mean absolute activation as importance
                activation = torch.mean(torch.abs(embedding)).item()
                activations.append(activation)

        if activations:
            mean_activation = np.mean(activations)
            max_activation = np.max(activations)
            total_activation = np.sum(activations)
        else:
            mean_activation = 0
            max_activation = 0
            total_activation = 0

        sample_scores.append({
            'sample_idx': sample_idx,
            'mean_activation': mean_activation,
            'max_activation': max_activation,
            'total_activation': total_activation,
            'true_value': y_test[sample_idx].item() if y_test is not None else None
        })

    df = pd.DataFrame(sample_scores)
    df = df.sort_values('total_activation', ascending=False)

    return df


def analyze_gene_contributions(model, x_test, lineage_path, top_samples, gene_id_mapping):
    """
    For top samples, identify which genes are most important for lineage activation.

    Returns:
        DataFrame: Gene importance scores
    """
    print("\nAnalyzing gene contributions for top samples...")

    # Reverse gene mapping (id -> gene name)
    id_to_gene = {v: k for k, v in gene_id_mapping.items()}

    gene_scores = defaultdict(lambda: {'activations': [], 'samples': []})

    # For each top sample
    for sample_idx in top_samples:
        sample_genotype = x_test[sample_idx]  # (n_genes,)

        # Find which genes are mutated/activated
        nonzero_genes = torch.nonzero(sample_genotype).squeeze()

        if nonzero_genes.dim() == 0:
            nonzero_genes = [nonzero_genes.item()]
        else:
            nonzero_genes = nonzero_genes.tolist()

        # For each activated gene
        for gene_idx in nonzero_genes:
            gene_name = id_to_gene.get(gene_idx, f"Gene_{gene_idx}")
            gene_value = sample_genotype[gene_idx].item()

            # Check if this gene is annotated to any lineage term
            gene_relevance = 0
            for term in lineage_path:
                if term in model.term_direct_gene_map:
                    term_genes = model.term_direct_gene_map[term]
                    if gene_idx in term_genes or gene_name in term_genes:
                        gene_relevance += 1

            gene_scores[gene_name]['activations'].append(gene_value)
            gene_scores[gene_name]['samples'].append(sample_idx)
            gene_scores[gene_name]['relevance'] = gene_relevance

    # Aggregate gene scores
    gene_summary = []
    for gene, data in gene_scores.items():
        gene_summary.append({
            'gene': gene,
            'n_samples': len(data['samples']),
            'mean_value': np.mean(data['activations']),
            'total_value': np.sum(data['activations']),
            'lineage_relevance': data['relevance'],
            'samples': data['samples'][:10]  # First 10 samples
        })

    df = pd.DataFrame(gene_summary)
    df = df.sort_values(['lineage_relevance', 'n_samples'], ascending=False)

    return df


def plot_sample_contributions(sample_df, output_file, top_n=20):
    """Plot sample contribution analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    top_samples = sample_df.head(top_n)

    # Plot 1: Total activation by sample
    axes[0, 0].barh(range(len(top_samples)), top_samples['total_activation'].values)
    axes[0, 0].set_yticks(range(len(top_samples)))
    axes[0, 0].set_yticklabels([f"Sample {idx}" for idx in top_samples['sample_idx']])
    axes[0, 0].set_xlabel('Total Lineage Activation')
    axes[0, 0].set_title(f'Top {top_n} Samples by Total Activation')
    axes[0, 0].invert_yaxis()

    # Plot 2: Mean vs max activation
    axes[0, 1].scatter(sample_df['mean_activation'], sample_df['max_activation'],
                      alpha=0.5, s=50)
    axes[0, 1].set_xlabel('Mean Activation')
    axes[0, 1].set_ylabel('Max Activation')
    axes[0, 1].set_title('Sample Activation Distribution')
    axes[0, 1].grid(True, alpha=0.3)

    # Plot 3: Activation distribution
    axes[1, 0].hist(sample_df['total_activation'], bins=50, edgecolor='black', alpha=0.7)
    axes[1, 0].set_xlabel('Total Activation')
    axes[1, 0].set_ylabel('Number of Samples')
    axes[1, 0].set_title('Distribution of Sample Activations')
    axes[1, 0].axvline(top_samples.iloc[0]['total_activation'],
                       color='red', linestyle='--', label='Top sample')
    axes[1, 0].legend()

    # Plot 4: True values vs activation (if available)
    if sample_df['true_value'].notna().any():
        axes[1, 1].scatter(sample_df['total_activation'], sample_df['true_value'],
                          alpha=0.5, s=50)
        axes[1, 1].set_xlabel('Total Lineage Activation')
        axes[1, 1].set_ylabel('True Target Value')
        axes[1, 1].set_title('Activation vs True Value')
        axes[1, 1].grid(True, alpha=0.3)
    else:
        axes[1, 1].text(0.5, 0.5, 'True values not available',
                       ha='center', va='center', transform=axes[1, 1].transAxes)
        axes[1, 1].axis('off')

    plt.suptitle('Sample Contribution Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved sample analysis: {output_file}")


def plot_gene_contributions(gene_df, output_file, top_n=30):
    """Plot gene contribution analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    top_genes = gene_df.head(top_n)

    # Plot 1: Top genes by lineage relevance
    relevant_genes = top_genes[top_genes['lineage_relevance'] > 0].head(20)
    if len(relevant_genes) > 0:
        axes[0, 0].barh(range(len(relevant_genes)), relevant_genes['lineage_relevance'].values,
                       color='steelblue')
        axes[0, 0].set_yticks(range(len(relevant_genes)))
        axes[0, 0].set_yticklabels(relevant_genes['gene'].values, fontsize=8)
        axes[0, 0].set_xlabel('Number of Lineage Terms Annotated')
        axes[0, 0].set_title('Genes Directly Annotated to Lineage Terms')
        axes[0, 0].invert_yaxis()
    else:
        axes[0, 0].text(0.5, 0.5, 'No genes directly annotated to lineage',
                       ha='center', va='center', transform=axes[0, 0].transAxes)

    # Plot 2: Top genes by sample frequency
    top_freq = gene_df.nlargest(20, 'n_samples')
    axes[0, 1].barh(range(len(top_freq)), top_freq['n_samples'].values, color='coral')
    axes[0, 1].set_yticks(range(len(top_freq)))
    axes[0, 1].set_yticklabels(top_freq['gene'].values, fontsize=8)
    axes[0, 1].set_xlabel('Number of Samples')
    axes[0, 1].set_title('Most Frequently Activated Genes')
    axes[0, 1].invert_yaxis()

    # Plot 3: Sample frequency vs lineage relevance
    axes[1, 0].scatter(gene_df['n_samples'], gene_df['lineage_relevance'],
                      alpha=0.6, s=80, c=gene_df['mean_value'], cmap='viridis')
    axes[1, 0].set_xlabel('Number of Samples')
    axes[1, 0].set_ylabel('Lineage Relevance')
    axes[1, 0].set_title('Gene Frequency vs Lineage Annotation')
    axes[1, 0].grid(True, alpha=0.3)

    # Plot 4: Gene value distribution
    axes[1, 1].hist(gene_df['mean_value'], bins=30, edgecolor='black', alpha=0.7)
    axes[1, 1].set_xlabel('Mean Gene Value')
    axes[1, 1].set_ylabel('Number of Genes')
    axes[1, 1].set_title('Distribution of Gene Values')

    plt.suptitle('Gene Contribution Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved gene analysis: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze which samples and genes drive lineage RLIPP scores',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-model', help='Trained model file (model.pt)',
                       type=str, required=True)
    parser.add_argument('-onto', help='Ontology file',
                       type=str, required=True)
    parser.add_argument('-gene2id', help='Gene to ID mapping file',
                       type=str, required=True)
    parser.add_argument('-test', help='Test data file (.pt)',
                       type=str, required=True)
    parser.add_argument('-lineages', help='Lineages TSV file',
                       type=str, required=True)
    parser.add_argument('-genotype_hiddens', help='Number of genotype hidden units',
                       type=int, default=6)
    parser.add_argument('-output', help='Output prefix',
                       type=str, default='drivers')
    parser.add_argument('-top_lineages', help='Number of top lineages to analyze',
                       type=int, default=3)
    parser.add_argument('-top_samples', help='Number of top samples to analyze per lineage',
                       type=int, default=20)

    args = parser.parse_args()

    print("="*70)
    print("LINEAGE DRIVER ANALYSIS")
    print("="*70)

    # Create minimal args for data wrapper
    class MinimalArgs:
        def __init__(self):
            # Required files
            self.onto = args.onto
            self.gene2id = args.gene2id

            # Model architecture
            self.genotype_hiddens = args.genotype_hiddens
            self.min_dropout_layer = 2
            self.dropout_fraction = 0.3

            # Training hyperparameters (dummy values, not used for inference)
            self.lr = 0.001
            self.wd = 0.001
            self.alpha = 0.3
            self.epoch = 100
            self.batchsize = 64

            # Output configuration (dummy values, not used)
            self.modeldir = '.'
            self.delta = 0.001
            self.metric_output = 'metrics.tsv'

            # Dataset configuration (dummy values, not used)
            self.train = None
            self.test = None
            self.testsetratio = 0.2

            # Optimization (not used for inference)
            self.optimize = 0

    minimal_args = MinimalArgs()

    # Load data wrapper (for ontology and gene mapping)
    print("\nInitializing data wrapper...")
    data_wrapper = TrainingDataWrapper(minimal_args, logger=None)

    # Load model
    model = load_model_and_data(args.model, data_wrapper)

    # Load test data
    x_test, y_test = load_test_data(args.test)

    # Load lineages
    print(f"\nLoading lineages from: {args.lineages}")
    lineages_df = pd.read_csv(args.lineages, sep='\t')
    print(f"Loaded {len(lineages_df)} lineages")

    # Extract term embeddings for all samples
    term_embeddings = extract_term_embeddings_per_sample(
        model, x_test, data_wrapper.cuda
    )

    # Analyze each top lineage
    for i, row in lineages_df.head(args.top_lineages).iterrows():
        print(f"\n{'='*70}")
        print(f"ANALYZING LINEAGE {i+1}")
        print(f"{'='*70}")

        # Parse lineage path
        lineage_path = row['path'].split(' -> ')
        print(f"Lineage: {row['leaf']}")
        print(f"Length: {len(lineage_path)} terms")
        print(f"Combined score: {row['combined_score']:.3f}")

        # Analyze sample contributions
        sample_df = analyze_sample_contributions(term_embeddings, lineage_path, y_test)

        # Save sample results
        sample_file = f"{args.output}_lineage{i+1}_samples.tsv"
        sample_df.to_csv(sample_file, sep='\t', index=False)
        print(f"\nSaved sample analysis: {sample_file}")

        # Plot sample contributions
        plot_file = f"{args.output}_lineage{i+1}_samples.png"
        plot_sample_contributions(sample_df, plot_file)

        # Analyze gene contributions for top samples
        top_sample_indices = sample_df.head(args.top_samples)['sample_idx'].tolist()
        gene_df = analyze_gene_contributions(
            model, x_test, lineage_path, top_sample_indices,
            data_wrapper.gene_id_mapping
        )

        # Save gene results
        gene_file = f"{args.output}_lineage{i+1}_genes.tsv"
        gene_df.to_csv(gene_file, sep='\t', index=False)
        print(f"Saved gene analysis: {gene_file}")

        # Plot gene contributions
        gene_plot_file = f"{args.output}_lineage{i+1}_genes.png"
        plot_gene_contributions(gene_df, gene_plot_file)

        # Print summary
        print(f"\nTop 5 samples:")
        for idx, sample_row in sample_df.head(5).iterrows():
            print(f"  Sample {sample_row['sample_idx']}: "
                  f"activation={sample_row['total_activation']:.3f}")

        print(f"\nTop 5 lineage-relevant genes:")
        relevant = gene_df[gene_df['lineage_relevance'] > 0].head(5)
        for idx, gene_row in relevant.iterrows():
            print(f"  {gene_row['gene']}: "
                  f"in {gene_row['n_samples']} samples, "
                  f"annotated to {gene_row['lineage_relevance']} terms")

    print("\n" + "="*70)
    print("DRIVER ANALYSIS COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()
