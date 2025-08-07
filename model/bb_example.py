import torch
from torch import nn
import pandas as pd
import numpy as np

#https://pytorch.org/tutorials/beginner/basics/buildmodel_tutorial.html
class blackbox(nn.Module):
    def __init__(self, input_size,
                  hidden_layers, 
                  activation_fn,
                  weight_initialization,
                  embedding_dim=64,
                  dropout_rate=0.3):

        super(blackbox, self).__init__() # my understanding is that inherits all the methods from the nn.Module class
        self.activation_fn = activation_fn 
        self.input_size = input_size
        self.embedding_dim = embedding_dim
        self.dropout_rate = dropout_rate

        self.gene_embedding = nn.Embedding(input_size, embedding_dim)

        self.layers = nn.ModuleList()  # Create a ModuleList for the layers
        #neurons = input_size

        self.layers.append(nn.Sequential(
            nn.Linear(embedding_dim * 2, embedding_dim * 2),
            nn.BatchNorm1d(embedding_dim * 2),
            self.activation_fn(),
            nn.Dropout(dropout_rate)
        ))

#        neurons = int(neurons * 0.5)

        # Add hidden layers with batch norm
        current_size = embedding_dim * 2
        for _ in range(hidden_layers):
            if current_size < 4:
                break
            next_size = int(current_size * 0.7)
            self.layers.append(nn.Sequential(
                nn.Linear(current_size, next_size),
                nn.BatchNorm1d(next_size),
                self.activation_fn(),
                nn.Dropout(dropout_rate)
            ))
            current_size = next_size

        # define the output layer
        self.output_layer = nn.Linear(current_size, 1)

        if weight_initialization:
            self.layers.apply(self.init_weights) #TODO add arg to pass act fn for correct weight initilization

    def forward(self, query_gene_idx, array_gene_idx):
        query_emb = self.gene_embedding(query_gene_idx)
        array_emb = self.gene_embedding(array_gene_idx)

        combined = torch.cat((query_emb, array_emb), dim = 1)
#        out = self.layers(input)
        x = combined
        for layer in self.layers:
            x = layer(x)

        out = self.output_layer(x)
        return out

    def init_weights(self, m):
        if isinstance(m, nn.Linear):
            # Kaiming initialization with larger variance to compensate for dropout
            # This is particularly important when using dropout
            nn.init.kaiming_normal_(m.weight, mode='fan_in', nonlinearity='leaky_relu')
            if m.bias is not None:
                nn.init.constant_(m.bias, 0.01)
                
