"""
Data Wrapper Module for VNN Training

This module handles all data loading and preprocessing for VNN training:
- Loading gene-to-ID mappings
- Loading and parsing biological ontologies (e.g., Gene Ontology)
- Constructing directed acyclic graphs (DAGs) of term hierarchies
- Validating ontology structure

The TrainingDataWrapper class serves as a central container for all
data and configuration needed by the VNN model.

Classes:
    TrainingDataWrapper: Main data container and loader
"""

import sys
import numpy as np
import pandas as pd
import networkx as nx
import networkx.algorithms.components.connected as nxacc
import networkx.algorithms.dag as nxadag
import time
import os
import torch
import utils


class TrainingDataWrapper():
    """
    Container class that loads and stores all data and configuration for VNN training.

    This class handles:
        1. Loading gene-to-ID mappings
        2. Loading and parsing ontology files (e.g., GO terms and their relationships)
        3. Creating a directed acyclic graph (DAG) representing the ontology hierarchy
        4. Mapping genes to ontology terms
        5. Validating ontology structure (single root, connected components)
        6. Storing all training hyperparameters

    Attributes:
        gene_id_mapping (dict): Maps gene names to integer IDs
        num_hiddens_genotype (int): Number of hidden neurons per term
        min_dropout_layer (int): Layer index to start dropout
        dropout_fraction (float): Dropout probability
        lr (float): Learning rate
        wd (float): Weight decay
        alpha (float): Weight for auxiliary losses
        epochs (int): Number of training epochs
        batchsize (int): Batch size
        cuda (torch.device): Device for computation (CPU or GPU)
        modeldir (str): Directory for saving models
        delta (float): Minimum improvement threshold for model saving
        dG (networkx.DiGraph): Directed graph of ontology
        root (str): Root term of ontology
        term_size_map (dict): Maps terms to number of annotated genes
        term_direct_gene_map (dict): Maps terms to directly annotated genes
        train (str): Path to training data
        test (str): Path to test data
        testsetratio (float): Test set ratio if splitting from training
        strfime (str): Timestamp string for file naming
        outfile (str): Path to metrics output file
    """

    def __init__(self, args, logger=None):
        """
        Initialize the data wrapper with configuration and load all necessary data.

        Args:
            args (argparse.Namespace): Parsed command-line arguments containing all configuration
            logger (logging.Logger, optional): Logger instance for tracking progress

        Raises:
            FileNotFoundError: If any required data files are not found
            ValueError: If ontology structure is invalid (multiple roots, disconnected components)
        """
        self.logger = logger

        if self.logger:
            self.logger.debug("Loading gene-to-ID mapping...")

        # Load gene to ID mapping
        self.gene_id_mapping = utils.load_mapping(args.gene2id, 'genes')

        if self.logger:
            self.logger.debug(f"Loaded {len(self.gene_id_mapping)} genes")

        # Model architecture parameters
        self.num_hiddens_genotype = args.genotype_hiddens
        self.min_dropout_layer = args.min_dropout_layer
        self.dropout_fraction = args.dropout_fraction

        # Training hyperparameters
        self.lr = args.lr
        self.wd = args.wd
        self.alpha = args.alpha
        self.epochs = args.epoch
        self.batchsize = args.batchsize

        # Device configuration
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.cuda = device

        if self.logger:
            self.logger.info(f"Using device: {device}")

        # Output configuration
        self.modeldir = args.modeldir
        self.delta = args.delta

        # Load ontology structure
        if self.logger:
            self.logger.info("Loading ontology...")

        self.load_ontology(args.onto)

        if self.logger:
            self.logger.info(f"Ontology loaded: {len(self.dG.nodes())} terms, root={self.root}")

        # Dataset configuration
        self.train = args.train
        self.test = args.test
        self.testsetratio = args.testsetratio

        # Optimization mode
        self.optimize = args.optimize

        # Output file naming
        self.strfime = time.strftime('%Y%m%d%H%M')
        self.outfile = os.path.join(self.modeldir, f"{self.strfime}_{args.metric_output}")

        if self.logger:
            self.logger.debug(f"Metrics will be saved to: {self.outfile}")


    def load_ontology(self, file_name):
        """
        Load and parse ontology file to create hierarchical structure.

        The ontology file should have three columns per line:
            1. Parent term
            2. Child term or gene name
            3. Relationship type ('default' for term-term, otherwise term-gene)

        This method:
            1. Builds a directed graph of term relationships
            2. Maps terms to their directly annotated genes
            3. Propagates gene annotations up the hierarchy
            4. Validates the ontology structure (single root, connected)

        Args:
            file_name (str): Path to ontology file

        Raises:
            FileNotFoundError: If ontology file doesn't exist
            ValueError: If ontology has multiple roots or disconnected components
            SystemExit: If any terms have no annotated genes

        Side Effects:
            Sets self.dG, self.root, self.term_size_map, self.term_direct_gene_map
        """
        if self.logger:
            self.logger.debug(f"Parsing ontology file: {file_name}")

        # Initialize data structures
        dG = nx.DiGraph()  # Directed graph for term hierarchy
        term_direct_gene_map = {}  # Maps term -> set of directly annotated gene IDs
        term_size_map = {}  # Maps term -> total number of genes (including children)
        gene_set = set()  # All genes mentioned in ontology

        # Parse ontology file
        file_handle = open(file_name)
        line_count = 0
        edge_count = 0
        gene_annotation_count = 0

        for line in file_handle:
            line_count += 1
            line = line.rstrip().split()

            if len(line) < 3:
                if self.logger:
                    self.logger.warning(f"Skipping malformed line {line_count}: {line}")
                continue

            # Line format: parent child relationship_type
            parent_term = line[0]
            child_or_gene = line[1]
            relationship_type = line[2]

            if relationship_type == 'default':
                # This is a term-to-term relationship
                dG.add_edge(parent_term, child_or_gene)
                edge_count += 1
            else:
                # This is a term-to-gene annotation
                # Only include genes that are in our gene mapping
                if child_or_gene not in self.gene_id_mapping:
                    continue

                if parent_term not in term_direct_gene_map:
                    term_direct_gene_map[parent_term] = set()

                term_direct_gene_map[parent_term].add(self.gene_id_mapping[child_or_gene])
                gene_set.add(child_or_gene)
                gene_annotation_count += 1

        file_handle.close()

        if self.logger:
            self.logger.debug(f"Parsed {line_count} lines")
            self.logger.debug(f"Created {edge_count} term-term edges")
            self.logger.debug(f"Found {gene_annotation_count} gene annotations")
            self.logger.debug(f"Total unique genes in ontology: {len(gene_set)}")

        # Propagate gene annotations up the hierarchy
        # Each term includes genes from all its descendants
        empty_terms = []

        for term in dG.nodes():
            term_gene_set = set()

            # Add directly annotated genes
            if term in term_direct_gene_map:
                term_gene_set = term_direct_gene_map[term]

            # Add genes from all descendant terms
            deslist = nxadag.descendants(dG, term)
            for child in deslist:
                if child in term_direct_gene_map:
                    term_gene_set = term_gene_set | term_direct_gene_map[child]

            # Check for empty terms
            if len(term_gene_set) == 0:
                if self.logger:
                    self.logger.warning(f'Empty term found: {term}')
                empty_terms.append(term)
            else:
                term_size_map[term] = len(term_gene_set)

        if empty_terms and self.logger:
            self.logger.warning(f"Found {len(empty_terms)} empty terms (will be skipped)")

        # Validate ontology structure
        # 1. Check for single root
        roots = [n for n in dG.nodes() if dG.in_degree(n) == 0]

        # 2. Check for connected components
        uG = dG.to_undirected()
        connected_subG_list = list(nxacc.connected_components(uG))

        # Log structure information
        if self.logger:
            self.logger.info(f"Ontology structure:")
            self.logger.info(f"  Number of terms: {len(dG.nodes())}")
            self.logger.info(f"  Number of roots: {len(roots)}")
            if len(roots) > 0:
                self.logger.info(f"  Root term: {roots[0]}")
            self.logger.info(f"  Connected components: {len(connected_subG_list)}")
            self.logger.info(f"  Empty terms: {len(empty_terms)}")

        # Validate structure
        if len(roots) > 1:
            error_msg = f'Ontology has {len(roots)} roots. Please use an ontology with a single root.'
            if self.logger:
                self.logger.error(error_msg)
            print('ERROR:', error_msg)
            print('Roots:', roots)
            sys.exit(1)

        if len(roots) == 0:
            error_msg = 'Ontology has no root terms. Please check the ontology file.'
            if self.logger:
                self.logger.error(error_msg)
            print('ERROR:', error_msg)
            sys.exit(1)

        if len(connected_subG_list) > 1:
            error_msg = f'Ontology has {len(connected_subG_list)} disconnected components. Please use a connected ontology.'
            if self.logger:
                self.logger.error(error_msg)
            print('ERROR:', error_msg)
            sys.exit(1)

        # Store ontology structure
        self.dG = dG
        self.root = roots[0]
        self.term_size_map = term_size_map
        self.term_direct_gene_map = term_direct_gene_map

        if self.logger:
            self.logger.debug("Ontology loaded and validated successfully")
