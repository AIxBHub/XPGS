import numpy as np
import pandas as pd
import time
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.data as du
from torch.autograd import Variable

import utils
import prepare_dataloader
from dcell_nn import DCellNN
from bb_nn import blackbox
from ccc_loss import *
import os

class ANNTrainer():

    def __init__(self, data_wrapper):
        self.data_wrapper = data_wrapper

    def train_epoch(self, train_loader, optimizer, loss_fn):
        self.model.train()
        running_loss = 0.0
        all_predictions = []
        all_labels = []
        for inputdata, labels in train_loader:
            features = inputdata.to(torch.float32).cuda(self.data_wrapper.cuda)
            labels = labels.to(torch.float32).unsqueeze(1).cuda(self.data_wrapper.cuda)
            optimizer.zero_grad()
            outputs = self.model(features)
            loss = loss_fn(outputs, labels)
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
            all_predictions.append(outputs.detach())
            all_labels.append(labels.detach())
        all_predictions = torch.cat(all_predictions, dim=0)
        all_labels = torch.cat(all_labels, dim=0)
        avg_loss = running_loss / len(train_loader) if len(train_loader) > 0 else 0.0
        
        return avg_loss, all_predictions, all_labels

    def validate_epoch(self, val_loader, loss_fn):
        self.model.eval()
        running_loss = 0.0
        all_predictions = []
        all_labels = []
        with torch.no_grad():
            for inputdata, labels in val_loader:
                features = inputdata.to(torch.float32).cuda(self.data_wrapper.cuda)
                labels = labels.to(torch.float32).unsqueeze(1).cuda(self.data_wrapper.cuda)
                outputs = self.model(features)
                loss = loss_fn(outputs, labels)
                running_loss += loss.item()
                all_predictions.append(outputs)
                all_labels.append(labels)
        all_predictions = torch.cat(all_predictions, dim=0)
        all_labels = torch.cat(all_labels, dim=0)
        avg_loss = running_loss / len(val_loader) if len(val_loader) > 0 else 0.0
        
        return avg_loss, all_predictions, all_labels 

    def train_model(self):
        loss_fn = nn.MSELoss()
        self.model = blackbox(self.data_wrapper)
        
        # Move model to GPU if available
        if torch.cuda.is_available() and hasattr(self.data_wrapper, 'cuda'):
            self.model.cuda(self.data_wrapper.cuda)

        train_loader, val_loader = prepare_dataloader.get_from_pt(
            self.data_wrapper.train,
            self.data_wrapper.test,
            self.data_wrapper.batchsize,
            self.data_wrapper.subset
        )

        optimizer = torch.optim.Adam(
            self.model.parameters(), 
            lr=self.data_wrapper.lr,
            betas=(0.9, 0.99),
            eps=1e-05,
            weight_decay=self.data_wrapper.wd
        )
        
        # Open log file for results
        logout = open(self.data_wrapper.outfile, "w")
        logout.write("epoch\ttrain_corr\ttrain_loss\tval_corr\tval_loss\telapsed_time\n")
        
        for epoch in range(self.data_wrapper.epochs):
            epoch_start_time = time.time()
            # Training phase
            train_loss, train_preds, train_labels = self.train_epoch(
                train_loader, optimizer, loss_fn
            )
            
            val_loss, val_preds, val_labels = self.validate_epoch(
                val_loader, loss_fn
            )

            ### get corrs
            val_corr = utils.pearson_corr(val_preds, val_labels)
            train_corr = utils.pearson_corr(train_preds, train_labels)

            elapsed_time = time.time() - epoch_start_time
            # Log results
            logout.write(f"{epoch}\t{train_corr:.6f}\t{train_loss:.6f}\t{val_corr:.6f}\t{val_loss:.6f}\t{elapsed_time:.2f}\n")
            logout.flush()
            # Print epoch summary
            print(f"Epoch {epoch}: Train Loss: {train_loss:.6f}, Train Corr: {train_corr:.6f}, "
                  f"Val Loss: {val_loss:.6f}, Val Corr: {val_corr:.6f}, Time: {elapsed_time:.2f}s")
  
        logout.close()
        

    