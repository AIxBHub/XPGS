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
    def __init__(self, args):
        self.gene_id_mapping = utils.load_mapping(args.gene2id, 'genes')
        self.num_hiddens_genotype = args.genotype_hiddens
        self.min_dropout_layer = args.min_dropout_layer
        self.dropout_fraction = args.dropout_fraction
        self.lr = args.lr
        self.wd = args.wd
        self.alpha = args.alpha
        self.epochs = args.epoch
        self.batchsize = args.batchsize
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.cuda = device
        self.modeldir = args.modeldir
        self.delta = args.delta
        self.load_ontology(args.onto)
        # Train test dataset sort out. If no test set given, use training dataset to make hold outs for test
        self.train = args.train
        self.test = args.test
        self.testsetratio = args.testsetratio
        #self.trainset, self.testset = prepare_dataloader.get_data(self.train, self.test, self.testsetratio, self.gene_id_mapping, self.batchsize)
            
        self.optimize = args.optimize
        self.strfime = time.strftime('%Y%m%d%H%M')
        self.outfile = os.path.join(self.modeldir,  f"{self.strfime}_{args.metric_output}")
    
    
    def load_ontology(self, file_name):
        dG = nx.DiGraph()
        term_direct_gene_map = {}
        term_size_map = {}
        gene_set = set()

        file_handle = open(file_name)
        for line in file_handle:
            line = line.rstrip().split()
            if line[2] == 'default':
                dG.add_edge(line[0], line[1])
            else:
                if line[1] not in self.gene_id_mapping:
                    continue
                if line[0] not in term_direct_gene_map:
                    term_direct_gene_map[line[0]] = set()
                term_direct_gene_map[line[0]].add(self.gene_id_mapping[line[1]])
                gene_set.add(line[1])
        file_handle.close()
        empty_terms = []
        for term in dG.nodes():
            term_gene_set = set()
            if term in term_direct_gene_map:
                term_gene_set = term_direct_gene_map[term]
            deslist = nxadag.descendants(dG, term)
            for child in deslist:
                if child in term_direct_gene_map:
                    term_gene_set = term_gene_set | term_direct_gene_map[child]
            # jisoo
            if len(term_gene_set) == 0:
                print('There is empty terms, please delete term:', term)
                empty_terms.append(term)
                sys.exit(1)
            else:
                term_size_map[term] = len(term_gene_set)
        print(empty_terms)        
        roots = [n for n in dG.nodes if dG.in_degree(n) == 0]

        uG = dG.to_undirected()
        connected_subG_list = list(nxacc.connected_components(uG))

        print('There are', len(roots), 'roots:', roots[0])
        print('There are', len(dG.nodes()), 'terms')
        print('There are', len(connected_subG_list), 'connected componenets')

        if len(roots) > 1:
            print('There are more than 1 root of ontology. Please use only one root.')
            sys.exit(1)
        if len(connected_subG_list) > 1:
            print('There are more than connected components. Please connect them.')
            sys.exit(1)

        self.dG = dG
        self.root = roots[0]
        self.term_size_map = term_size_map
        self.term_direct_gene_map = term_direct_gene_map
        
        
    