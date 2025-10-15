"""
DCell Neural Network Architecture Module

This module implements the hierarchical neural network architecture based on biological ontologies.
The network structure mirrors the directed acyclic graph (DAG) of the ontology, where each
biological term/subsystem becomes a layer in the network.

Key Features:
    - Hierarchical architecture guided by ontology structure
    - Masked connections that respect gene-term relationships
    - Bottom-up information flow from genes through terms to root
    - Auxiliary outputs at each term for multi-task learning
    - Interpretable intermediate representations

Based on DCell architecture from:
Ma et al. (2018). "Using deep learning to model the hierarchical structure and function of a cell."
Nature Methods, 15(4), 290-298. PMID: 29505029

Classes:
    DCellNN: Main neural network class implementing hierarchical architecture
"""

import sys
import os
import copy
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class DCellNN(nn.Module):
    """
    Hierarchical Neural Network based on Biological Ontology (DCell architecture).

    This network builds a multi-layer architecture where:
        1. Input layer: Gene-level features
        2. Hidden layers: Biological terms/subsystems organized hierarchically
        3. Output layer: Final phenotype prediction

    Each term in the ontology becomes a "neuron subsystem" that:
        - Receives input from its directly annotated genes
        - Receives input from all its child terms
        - Produces a hidden representation (embedding)
        - Generates an auxiliary prediction for multi-task learning

    Architecture Details:
        - Gene-to-term connections are masked to respect annotations
        - Term-to-term connections follow parent-child relationships
        - Each term has: linear layer → tanh → batch norm → auxiliary output
        - Dropout applied from specified layer onwards
        - Final layer aggregates root term representation

    Attributes:
        root (str): Root term of the ontology
        num_hiddens_genotype (int): Number of hidden neurons per term
        input_size_map (dict): Maps each term to its input dimension
        term_direct_gene_map (dict): Maps terms to directly annotated gene IDs
        min_dropout_layer (int): Layer index to start applying dropout
        dropout_fraction (float): Dropout probability
        term_dim_map (dict): Maps terms to their hidden dimension
        gene_id_mapping (dict): Maps gene names to IDs
        gene_dim (int): Total number of genes (input dimension)
        term_layer_list (list): List of term lists, organized bottom-up
        term_neighbor_map (dict): Maps each term to its children
    """

    def __init__(self, data_wrapper):
        """
        Initialize the hierarchical neural network based on ontology.

        Args:
            data_wrapper (TrainingDataWrapper): Container with ontology and configuration

        Side Effects:
            - Constructs all network layers based on ontology structure
            - Adds modules for each term's processing
            - Creates final output layer
        """
        super().__init__()

        # Store configuration
        self.root = data_wrapper.root
        self.num_hiddens_genotype = data_wrapper.num_hiddens_genotype

        # Initialize input size map
        self.input_size_map = {}

        # Store term-gene relationships
        self.term_direct_gene_map = data_wrapper.term_direct_gene_map

        # Dropout configuration
        self.min_dropout_layer = data_wrapper.min_dropout_layer
        self.dropout_fraction = data_wrapper.dropout_fraction

        # Calculate dimensions for each term
        self.cal_term_dim(data_wrapper.term_size_map)

        # Store gene mapping
        self.gene_id_mapping = data_wrapper.gene_id_mapping
        # Total number of genes (input dimension)
        self.gene_dim = len(self.gene_id_mapping)

        # Build network architecture
        # 1. Add gene-to-term layers
        self.contruct_direct_gene_layer()

        # 2. Build term hierarchy layers
        self.construct_NN_graph(copy.deepcopy(data_wrapper.dG))

        # 3. Add final output layers
        self.add_module('final_aux_linear_layer', nn.Linear(data_wrapper.num_hiddens_genotype, 1))
        self.add_module('final_linear_layer_output', nn.Linear(1, 1))


    def cal_term_dim(self, term_size_map):
        """
        Calculate the hidden dimension for each term.

        Currently, all terms use the same hidden dimension (num_hiddens_genotype).
        This could be modified to use term-specific dimensions based on term_size_map.

        Args:
            term_size_map (dict): Maps term IDs to number of annotated genes

        Side Effects:
            Sets self.term_dim_map with hidden dimensions for each term
        """
        self.term_dim_map = {}

        for term, term_size in term_size_map.items():
            # All terms use the same number of hidden neurons
            num_output = self.num_hiddens_genotype

            # Store the hidden dimension for this term
            num_output = int(num_output)
            self.term_dim_map[term] = num_output


    def contruct_direct_gene_layer(self):
        """
        Build gene-to-term connection layers.

        For each term, creates a linear layer that:
            - Takes all genes as input (gene_dim)
            - Outputs only the genes directly annotated with the term
            - Will be masked during training to enforce sparse connections

        These layers extract gene-level features relevant to each term.

        Side Effects:
            Adds modules named '{term}_direct_gene_layer' for each term

        Raises:
            SystemExit: If any term has no directly annotated genes
        """
        for term, gene_set in self.term_direct_gene_map.items():
            if len(gene_set) == 0:
                print(f'ERROR: No directly annotated genes for term: {term}')
                sys.exit(1)

            # Create layer: all genes → directly annotated genes for this term
            # Note: This layer will be masked during training to enforce sparse connections
            self.add_module(
                term + '_direct_gene_layer',
                nn.Linear(self.gene_dim, len(gene_set))
            )


    def construct_NN_graph(self, dG):
        """
        Build the hierarchical neural network layers based on ontology DAG.

        Constructs layers bottom-up (leaves to root) by:
            1. Finding leaf nodes (no children)
            2. Creating layers for each term
            3. Removing processed nodes
            4. Repeating until all nodes processed

        Each term gets:
            - Optional dropout layer
            - Linear transformation layer
            - Batch normalization layer
            - Two auxiliary output layers (for multi-task learning)

        Args:
            dG (networkx.DiGraph): Directed graph of ontology (will be modified)

        Side Effects:
            - Sets self.term_layer_list: List of term layers (bottom to top)
            - Sets self.term_neighbor_map: Maps terms to their children
            - Sets self.input_size_map: Maps terms to input dimensions
            - Adds modules for each term's processing layers
        """
        self.term_layer_list = []  # Stores layers from bottom (leaves) to top (root)
        self.term_neighbor_map = {}  # Maps each term to its children

        # Build neighbor map (term -> list of children)
        for term in dG.nodes():
            self.term_neighbor_map[term] = []
            for child in dG.neighbors(term):
                self.term_neighbor_map[term].append(child)

        # Build layers bottom-up
        layer_idx = 0
        while True:
            # Find leaf nodes (nodes with no outgoing edges)
            leaves = [n for n in dG.nodes() if dG.out_degree(n) == 0]

            if len(leaves) == 0:
                # All nodes processed
                break

            # Add this layer to the list
            self.term_layer_list.append(leaves)

            # Create neural network modules for each term in this layer
            for term in leaves:
                # Calculate input size for this term
                # Input = concatenation of:
                #   1. Embeddings from all child terms
                #   2. Direct gene features for this term
                input_size = 0

                # Add dimensions from child terms
                for child in self.term_neighbor_map[term]:
                    input_size += self.term_dim_map[child]

                # Add dimensions from directly annotated genes
                if term in self.term_direct_gene_map:
                    input_size += len(self.term_direct_gene_map[term])

                self.input_size_map[term] = input_size

                # Get hidden dimension for this term
                term_hidden = self.term_dim_map[term]

                # Add dropout layer if we're past the minimum dropout layer
                if layer_idx >= self.min_dropout_layer:
                    self.add_module(
                        term + '_dropout_layer',
                        nn.Dropout(p=self.dropout_fraction)
                    )

                # Main transformation layer: input → hidden representation
                self.add_module(
                    term + '_linear_layer',
                    nn.Linear(input_size, term_hidden)
                )

                # Batch normalization for hidden representation
                self.add_module(
                    term + '_batchnorm_layer',
                    nn.BatchNorm1d(term_hidden)
                )

                # Auxiliary output layers (for multi-task learning)
                # These allow the network to make predictions at each subsystem level
                self.add_module(
                    term + '_aux_linear_layer1',
                    nn.Linear(term_hidden, 1)
                )
                self.add_module(
                    term + '_aux_linear_layer2',
                    nn.Linear(1, 1)
                )

            # Move to next layer
            layer_idx += 1

            # Remove processed nodes from graph
            dG.remove_nodes_from(leaves)


    def forward(self, gene_input):
        """
        Forward pass through the hierarchical network.

        Information flows bottom-up through the ontology hierarchy:
            1. Gene features are extracted for each term
            2. Leaf terms process their gene features
            3. Parent terms process concatenated [children embeddings, gene features]
            4. Process continues up to root term
            5. Final layer produces phenotype prediction

        At each term, the network:
            - Concatenates inputs from children and genes
            - Applies dropout (if configured)
            - Linear transformation → tanh activation → batch norm
            - Generates auxiliary prediction

        Args:
            gene_input (torch.Tensor): Gene-level features, shape (batch_size, gene_dim)

        Returns:
            tuple: (aux_out_map, hidden_embeddings_map)
                - aux_out_map (dict): Maps term names to auxiliary predictions
                    Includes 'final' key for final output
                - hidden_embeddings_map (dict): Maps term names to hidden representations
        """
        # Initialize output dictionaries
        hidden_embeddings_map = {}  # Store hidden representations for each term
        aux_out_map = {}  # Store auxiliary predictions for each term
        term_gene_out_map = {}  # Store gene features for each term

        # Step 1: Extract gene features for each term
        # Each term gets its directly annotated gene features
        for term, _ in self.term_direct_gene_map.items():
            term_gene_out_map[term] = self._modules[term + '_direct_gene_layer'](gene_input)

        # Step 2: Process terms layer by layer, bottom-up
        for layer_idx, layer in enumerate(self.term_layer_list):

            for term in layer:
                # Collect inputs for this term
                child_input_list = []

                # Add embeddings from all child terms
                for child in self.term_neighbor_map[term]:
                    child_input_list.append(hidden_embeddings_map[child])

                # Add gene features for this term
                if term in self.term_direct_gene_map:
                    child_input_list.append(term_gene_out_map[term])

                # Concatenate all inputs
                child_input = torch.cat(child_input_list, dim=1)

                # Apply dropout if configured for this layer
                if layer_idx >= self.min_dropout_layer:
                    dropout_out = self._modules[term + '_dropout_layer'](child_input)
                    term_NN_out = self._modules[term + '_linear_layer'](dropout_out)
                else:
                    term_NN_out = self._modules[term + '_linear_layer'](child_input)

                # Apply activation and batch normalization
                Tanh_out = torch.tanh(term_NN_out)
                hidden_embeddings_map[term] = self._modules[term + '_batchnorm_layer'](Tanh_out)

                # Generate auxiliary prediction for this term
                aux_layer1_out = torch.tanh(
                    self._modules[term + '_aux_linear_layer1'](hidden_embeddings_map[term])
                )
                aux_out_map[term] = self._modules[term + '_aux_linear_layer2'](aux_layer1_out)

        # Step 3: Generate final output from root term embedding
        final_input = hidden_embeddings_map[self.root]
        aux_layer_out = torch.tanh(self._modules['final_aux_linear_layer'](final_input))
        aux_out_map['final'] = self._modules['final_linear_layer_output'](aux_layer_out)

        return aux_out_map, hidden_embeddings_map
