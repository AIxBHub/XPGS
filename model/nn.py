import torch
from torch import nn
import torch.nn.functional as F
import pandas as pd
import numpy as np
import synthetic

### draft -- with gpt4.1 help, mostly reviewed but need to test

#class genes_input(nn.Module):
#    def __init__(self, ins, outs, mask, activation_fn):
#        """
#        masks is a list of ontology based masks
#        """
#        super().__init__()
#        self.activation_fn = activation_fn
#        self.input_size = ins
#        self.out_size = outs
#        self.register_buffer('mask',mask) ## model params that should be saved and restored in state_dict but not used for training
#        self.weight = nn.Parameters(torch.zeros(outs, ins))
#        nn.init.xavier_uniform_(self.weight.data)
#        self.weight.data *= self.mask
#
#    def forward(self, x):
#        masked_weights = self.weight * self.mask
#        return F.linear(x, masked_weights)

class SparseNet(nn.Module):
    def __init__(self,
                 activation_fn, 
                 masks,
                 in_size, ## encoded genes
                 out_size = 1 ## default to 1 for prediction 
                 ):
        """
        NN class for sparse network.
        Masks is a list of binary masks derived from graph ontology.
        The total number of hidden layers is defined by the total number of masks. 
        """
        super().__init__()
        self.activation_fn = activation_fn
        self.num_layers = len(masks) ## defines the total number of layers
        layers = []
        ## loop through all masks and create a new sparse layer for each
        ## uses helper function _create_sparese_layer
        ## each hidden layer is the size of the binary mask
        for i, mask in enumerate(masks): 
            hidden_size = len(mask)
            layers.append(self._create_sparse_layer(in_size, hidden_size, mask))
            in_size = hidden_size

        self.sparse_layers = nn.ModuleList(layers)
        self.final = nn.Linear(in_size, out_size)

    def _create_sparse_layer(self, in_dim, out_dim, mask):
        """
        helper function to make the sparse layer
        for each layer, takes the mask and sets initial weight and mask params, without applying the mask
        """
        weight = nn.Parameter(torch.zeros(out_dim, in_dim))
        nn.init.xavier_uniform_(weight.data)
        weight.data *= mask 

        layer = nn.Linear(in_dim, out_dim, bias=False)
        layer.weight = weight
        layer.register_buffer('mask', mask)

        return layer

    def forward(self, x):
        """
        on each forward pass the layer weights are masked with the params set by _create_sparse_layer
        """
        for layer in self.sparse_layers:
            x = self.activation_fn(F.linear(x, layer.weight * layer.mask))
        x = self.final(x)
        return torch.sigmoid(x)


if __name__ == "__main__":

    X, y = synthetic.generate_synthetic_data(num_samples=100, num_features=10)

    # Initialize the network
    masks = [(torch.rand(5,10)>0.2).float(), (torch.rand(3,5)>0.2).float()]
    print(masks)
#    masks = [torch.ones(5, 10), torch.ones(3, 5)]  # Example masks, replace with actual if needed
    net = SparseNet(activation_fn=torch.relu, in_size=10, out_size=1, masks=masks)

    # Perform a forward pass
    output = net(X)
    print("Network output:", output[:5])