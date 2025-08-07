import torch
from torch import nn
import torch.nn.functional as F
import pandas as pd
import numpy as np

### draft -- with gpt4.1 help, mostly reviewed but need to test

class genes_input(nn.Module):
    def __init__(self, ins, outs, mask, activation_fn):
        """
        masks is a list of ontology based masks
        """
        super().__init__()
        self.activation_fn = activation_fn
        self.input_size = ins
        self.out_size = outs
        self.register_buffer('mask',mask) ## model params that should be saved and restored in state_dict but not used for training
        self.weight = nn.Parameters(torch.zeros(outs, ins))
        nn.init.xavier_uniform_(self.weight.data)
        self.weight.data *= self.mask

    def forward(self, x):
        masked_weights = self.weight * self.mask
        return F.linear(x, masked_weights)

class sparse_net(nn.Module):
    def __init__(self, in_size, out_size, masks):
        super().__init__()
        self.num_layers = len(masks)
        layers = []
        for i, (h_dim, mask) in enumerate(zip(len(masks), masks)):
            layers.append(genes_input(in_dim, h_dim, mask))
            in_dim = h_dim
        self.sparse_layers = nn.ModuleList(layers)
        self.final = nn.Linear(in_dim, out_size)
    
    def forward(self, x):
        for layer in self.sparse_layers:
            x = torch.relu(layer(x))
        x = self.final(x)
        return torch.sigmoid(x)