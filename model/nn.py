import torch
from torch import nn
import pandas as pd
import numpy as np

class sparse(nn.Module):
    def __init__(self, input_size, out_size, mask, activation_fn):
        super(sparse, self).__init__()
        self.activation_fn = activation_fn
        self.input_size = input_size
        self.out_size = out_size
        self.register_buffer('mask',mask) ## model params that should be saved and restored in state_dict but not used for training
        self.layers = nn.ModuleList() ## create a module list to hold layers

        self.layers.append(nn.Sequential(
            nn.Linear
        ))