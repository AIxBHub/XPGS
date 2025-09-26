import sys
import os
import copy
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

class DCellNN(nn.Module):

    def __init__(self, 
                 data_wrapper, 
                 hidden_layers = 10, 
                 dropout_fraction = 0.1,
                 activation_fn = 'tanh'):

        super().__init__()
        self.activation_fn = activation_fn
        self.root = data_wrapper.root
        # input size map initiation
        self.input_size_map = {}
        # Dropout Params
        self.min_dropout_layer = data_wrapper.min_dropout_layer
        self.dropout_fraction = data_wrapper.dropout_fraction

        # This map is still missing. Need a function to get it from input file. This map needs to include genes from train AND test sets.
        self.gene_id_mapping = data_wrapper.gene_id_mapping  
        # ngenes, gene_dim are the number of all genes
        self.gene_dim = len(self.gene_id_mapping)

        self.layers = nn.ModuleList()

        ## input layer
        first_hidden_size = self.gene_dim * 0.5 ## first hidden 50% of input
        self.layers.append(
            nn.Squential(
                nn.Linear(self.gene_dim, first_hidden_size),
                nn.BatchNorm1d(first_hidden_size),
                nn.Dropout(p = dropout_fraction)
        ))

        ## hidden layers
        hidden_size = first_hidden_size
        while hidden_size >= 4:
            for layer in hidden_layers:
                self.layers.append(
                    nn.Squential(
                        nn.Linear(hidden_size, hidden_size * 0.5),
                        nn.BatchNorm1d(hidden_size * 0.5),
                        nn.Dropout(p = dropout_fraction) 
                ))
            hidden_size = hidden_size * 0.5

        # add module for final layer
#        self.add_module('final_aux_linear_layer', nn.Linear(data_wrapper.num_hiddens_genotype, 1))
        self.output_layer(, 1))

        
    # definition of forward function
    def forward(self, gene_input):

        hidden_embeddings_map = {}
        aux_out_map = {}
        term_gene_out_map = {}
        for term, _ in self.term_direct_gene_map.items():
            term_gene_out_map[term] = self._modules[term + '_direct_gene_layer'](gene_input)

        for i, layer in enumerate(self.term_layer_list):

            for term in layer:

                child_input_list = []
                for child in self.term_neighbor_map[term]:
                    child_input_list.append(hidden_embeddings_map[child])

                if term in self.term_direct_gene_map:
                    child_input_list.append(term_gene_out_map[term])

                child_input = torch.cat(child_input_list, 1)
                if i >= self.min_dropout_layer:
                    dropout_out = self._modules[term + '_dropout_layer'](child_input)
                    term_NN_out = self._modules[term + '_linear_layer'](dropout_out)
                else:
                    term_NN_out = self._modules[term + '_linear_layer'](child_input)
                Tanh_out = torch.tanh(term_NN_out)
                hidden_embeddings_map[term] = self._modules[term + '_batchnorm_layer'](Tanh_out)
                aux_layer1_out = torch.tanh(self._modules[term + '_aux_linear_layer1'](hidden_embeddings_map[term]))
                aux_out_map[term] = self._modules[term + '_aux_linear_layer2'](aux_layer1_out)

        final_input = hidden_embeddings_map[self.root]
        aux_layer_out = torch.tanh(self._modules['final_aux_linear_layer'](final_input))
        aux_out_map['final'] = self._modules['final_linear_layer_output'](aux_layer_out)

        return aux_out_map, hidden_embeddings_map
