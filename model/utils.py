"""
Utility Functions for VNN Training

This module provides utility functions for:
- Loading gene/term mappings from files
- Creating masks for ontology-based sparse connections
- Building input vectors from genotype data
- Calculating Pearson correlation for evaluation

Functions:
    load_mapping: Load gene or term ID mappings from file
    create_term_mask: Create binary masks for enforcing ontology structure
    build_input_vector: Convert genotype indices to one-hot vectors
    pearson_corr: Calculate Pearson correlation coefficient
"""

import torch


def load_mapping(mapping_file, mapping_type):
    """
    Load ID mapping from a tab-separated file.

    The file should have two columns:
        Column 1: Name (gene name, term name, etc.)
        Column 2: Integer ID

    Args:
        mapping_file (str): Path to mapping file
        mapping_type (str): Type of mapping (for logging, e.g., 'genes', 'terms')

    Returns:
        dict: Maps names (str) to integer IDs (int)

    Example:
        Gene mapping file:
            BRCA1   0
            TP53    1
            MYC     2
    """
    mapping = {}
    file_handle = open(mapping_file)

    for line in file_handle:
        line = line.rstrip().split()
        # Format: name  id
        mapping[line[0]] = int(line[1])

    file_handle.close()

    print(f'Total number of {mapping_type} = {len(mapping)}')
    return mapping


def create_term_mask(term_direct_gene_map, gene_dim, cuda_id):
    """
    Create binary masks for gene-to-term connections.

    For each term, creates a mask matrix that enforces sparse connections:
    only genes directly annotated with the term can connect to it.

    Args:
        term_direct_gene_map (dict): Maps term names to sets of gene IDs
        gene_dim (int): Total number of genes (input dimension)
        cuda_id (torch.device): Device for tensor placement

    Returns:
        dict: Maps term names to mask tensors
            Each mask has shape (num_annotated_genes, gene_dim)
            mask[i, j] = 1 if gene j is the i-th annotated gene for this term
            mask[i, j] = 0 otherwise

    Example:
        If term 'GO:0001234' has genes [5, 12, 23] annotated:
        mask shape: (3, gene_dim)
        mask[0, 5] = 1   # First annotated gene
        mask[1, 12] = 1  # Second annotated gene
        mask[2, 23] = 1  # Third annotated gene
        All other entries = 0

    This mask is element-wise multiplied with the weight matrix to enforce
    that only relevant genes contribute to each term.
    """
    term_mask_map = {}

    for term, gene_set in term_direct_gene_map.items():
        # Create zero mask: (num_annotated_genes, total_genes)
        mask = torch.zeros(len(gene_set), gene_dim, requires_grad=False).to(cuda_id)

        # Set mask[i, gene_id] = 1 for each annotated gene
        for i, gene_id in enumerate(gene_set):
            mask[i, gene_id] = 1

        term_mask_map[term] = mask

    return term_mask_map


def build_input_vector(inputdata, gene_id_mapping):
    """
    Convert gene indices to one-hot encoded vectors.

    This function is used when input data contains gene indices
    rather than binary gene presence/absence vectors.

    Args:
        inputdata (torch.Tensor): Gene indices, shape (batch_size, num_genes_per_sample)
        gene_id_mapping (dict): Maps gene names to IDs (used for dimension)

    Returns:
        torch.Tensor: One-hot encoded features, shape (batch_size, gene_dim)
            features[i, j] = 1 if gene j is present in sample i
            features[i, j] = 0 otherwise

    Example:
        If sample has genes [2, 5, 7] and gene_dim = 10:
        output[0, 2] = 1
        output[0, 5] = 1
        output[0, 7] = 1
        All other values = 0
    """
    # Convert to int64 for scatter operation
    inputdata = inputdata.to(torch.int64)

    genedim = len(gene_id_mapping)
    sample_num = inputdata.size()[0]

    # Initialize zero matrix
    features = torch.zeros((sample_num, genedim), requires_grad=False)

    # Scatter 1s at gene indices
    # scatter_(dim, index, value)
    features.scatter_(1, inputdata, 1)

    return features


def pearson_corr(x, y):
    """
    Calculate Pearson correlation coefficient between two tensors.

    The Pearson correlation measures the linear relationship between
    predictions and true values, ranging from -1 to 1:
        1: Perfect positive correlation
        0: No linear correlation
       -1: Perfect negative correlation

    Args:
        x (torch.Tensor): First tensor (e.g., predictions)
        y (torch.Tensor): Second tensor (e.g., true labels)

    Returns:
        torch.Tensor: Scalar Pearson correlation coefficient

    Formula:
        corr = Σ((x - x̄)(y - ȳ)) / (||x - x̄|| * ||y - ȳ||)

    where x̄ and ȳ are means of x and y respectively.
    """
    # Center the data (subtract means)
    xx = x - torch.mean(x)
    yy = y - torch.mean(y)

    # Calculate correlation
    # numerator: covariance
    # denominator: product of standard deviations (L2 norms)
    return torch.sum(xx * yy) / (torch.norm(xx, 2) * torch.norm(yy, 2))
