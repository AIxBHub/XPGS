import sys
import os
import copy
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

class blackbox(nn.Module):

    def __init__(self, data_wrapper):

        super().__init__()
        self.activation_fn = nn.Tanh()
        #self.root = data_wrapper.root
        self.num_hiddens_genotype = data_wrapper.num_hiddens_genotype

        self.dropout_fraction = data_wrapper.dropout_fraction
        self.gene_id_mapping = data_wrapper.gene_id_mapping  
        # ngenes, gene_dim are the number of all genes
        self.gene_dim = len(self.gene_id_mapping)

        self.layers = nn.ModuleList()

        ## input layer
        first_input_layer = self.gene_dim
        print(first_input_layer)
        self.layers.append(nn.Sequential(
            nn.Linear(first_input_layer, int(first_input_layer * 0.5)),
            nn.BatchNorm1d(int(first_input_layer * 0.5)),
            nn.Dropout(self.dropout_fraction)
        ))

        ## add hidden layers
        hidden = int(first_input_layer * 0.5)
        for l in range(data_wrapper.hidden_layers):
            next_hidden = int(hidden * 0.5)
            while hidden >= 4:
                self.layers.append(nn.Sequential(
                    nn.Linear(hidden, next_hidden),
                    nn.BatchNorm1d(next_hidden),
                    nn.Dropout(self.dropout_fraction)
                ))
                hidden = next_hidden 
                next_hidden = int(hidden * 0.5)

        self.output_layer = nn.Linear(hidden, 1)
        # add module for final layer
#        self.add_module('final_aux_linear_layer', nn.Linear(data_wrapper.num_hiddens_genotype, 1))
#        self.add_module('final_linear_layer_output', nn.Linear(1, 1))

    def forward(self, x):
        for layer in self.layers:
            x = layer(x)
            x = self.activation_fn(x)
        x = self.output_layer(x)
        return(x)