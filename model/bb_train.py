import numpy as np
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


    def train_epoch(self, train_loader, optimizer, loss_fn, epoch):
        """
        Execute a single training epoch
        
        Args:
            train_loader: DataLoader for training data
            optimizer: Optimizer for updating model parameters
            loss_fn: Loss function
            epoch: Current epoch number
            
        Returns:
            tuple: (average_loss, predictions, true_labels)
        """
        self.model.train()
        running_loss = 0.0
        all_predictions = []
        all_labels = []
        
        for i, (inputdata, labels) in enumerate(train_loader):
            # Prepare input features using the utility function
            features = utils.build_input_vector(inputdata, self.data_wrapper.gene_id_mapping)
            features = features.to(torch.float32).cuda(self.data_wrapper.cuda)
            labels = labels.to(torch.float32).unsqueeze(1).cuda(self.data_wrapper.cuda)
            
            # Zero gradients
            optimizer.zero_grad()
            
            # Forward pass
            outputs = self.model(features)
            # Compute loss
            loss = loss_fn(outputs, labels)
            
            # Backward pass
            loss.backward()
            
            # Update parameters
            optimizer.step()
            
            # Track metrics
            running_loss += loss.item()
            all_predictions.append(outputs.detach())
            all_labels.append(labels.detach())
            
            # Print progress every 100 batches
            if i % 100 == 99:
                avg_loss = running_loss / 100
                print(f'  Epoch {epoch}, Batch {i + 1}, Loss: {avg_loss:.6f}')
                running_loss = 0.0
        
        # Concatenate all predictions and labels for correlation calculation
        all_predictions = torch.cat(all_predictions, dim=0)
        all_labels = torch.cat(all_labels, dim=0)
        
        # Calculate average loss for the epoch
        avg_epoch_loss = running_loss / len(train_loader) if len(train_loader) > 0 else 0.0
        
        return avg_epoch_loss, all_predictions, all_labels

    def train_model(self):
        """
        Main training method for the ANN trainer
        """
        # Initialize the model
        self.model = blackbox(self.data_wrapper)
        
        # Move model to GPU if available
        if torch.cuda.is_available() and hasattr(self.data_wrapper, 'cuda'):
            self.model.cuda(self.data_wrapper.cuda)

        # Get data loaders
        train_loader, val_loader = prepare_dataloader.get_data(
            self.data_wrapper.train,
            self.data_wrapper.test,
            self.data_wrapper.testsetratio,
            self.data_wrapper.gene_id_mapping,
            self.data_wrapper.batchsize
        )

        # Initialize optimizer
        optimizer = torch.optim.Adam(
            self.model.parameters(), 
            lr=self.data_wrapper.lr,
            betas=(0.9, 0.99),
            eps=1e-05,
            weight_decay=self.data_wrapper.wd
        )
        
        loss_fn = nn.MSELoss()
        
        # Open log file for results
        logout = open(self.data_wrapper.outfile, "w")
        logout.write("epoch\ttrain_corr\ttrain_loss\tval_corr\tval_loss\telapsed_time\n")
        
        best_val_corr = -float('inf')
        
        print(f"Starting training for {self.data_wrapper.epochs} epochs...")
        
        for epoch in range(self.data_wrapper.epochs):
            epoch_start_time = time.time()
            
            # Training phase
            train_loss, train_preds, train_labels = self.train_epoch(
                train_loader, optimizer, loss_fn, epoch
            )
            
            # Calculate training correlation
            train_corr = self.calculate_correlation(train_preds, train_labels)
            
            # Validation phase
            val_loss, val_corr = self.validate_epoch(val_loader, loss_fn)
            
            # Calculate elapsed time
            elapsed_time = time.time() - epoch_start_time
            
            # Log results
            logout.write(f"{epoch}\t{train_corr:.6f}\t{train_loss:.6f}\t{val_corr:.6f}\t{val_loss:.6f}\t{elapsed_time:.2f}\n")
            logout.flush()
            
            # Print epoch summary
            print(f"Epoch {epoch}: Train Loss: {train_loss:.6f}, Train Corr: {train_corr:.6f}, "
                  f"Val Loss: {val_loss:.6f}, Val Corr: {val_corr:.6f}, Time: {elapsed_time:.2f}s")
            
            ## Save best model
            #if val_corr > best_val_corr:
            #    best_val_corr = val_corr
            #    torch.save(self.model.state_dict(), f"{self.data_wrapper.outfile}_best_model.pt")
            #    print(f"New best validation correlation: {best_val_corr:.6f}")
        
        logout.close()
        
    def validate_epoch(self, val_loader, loss_fn):
        """
        Execute validation for one epoch
        
        Returns:
            tuple: (average_loss, correlation)
        """
        self.model.eval()
        running_loss = 0.0
        all_predictions = []
        all_labels = []
        
        with torch.no_grad():
            for inputdata, labels in val_loader:
                # Prepare input features
                features = utils.build_input_vector(inputdata, self.data_wrapper.gene_id_mapping)
                features = features.to(torch.float32).cuda(self.data_wrapper.cuda)
                labels = labels.to(torch.float32).unsqueeze(1).cuda(self.data_wrapper.cuda)
                
                # Forward pass
                outputs = self.model(features)
                
                # Compute loss
                loss = loss_fn(outputs, labels)
                running_loss += loss.item()
                
                # Store predictions and labels
                all_predictions.append(outputs)
                all_labels.append(labels)
        
        # Concatenate all predictions and labels
        all_predictions = torch.cat(all_predictions, dim=0)
        all_labels = torch.cat(all_labels, dim=0)
        
        # Calculate average loss and correlation
        avg_loss = running_loss / len(val_loader) if len(val_loader) > 0 else 0.0
        correlation = self.calculate_correlation(all_predictions, all_labels)
        
        return avg_loss, correlation
    
    def calculate_correlation(self, predictions, labels):
        """
        Calculate Pearson correlation coefficient between predictions and labels
        
        Returns:
            float: correlation coefficient
        """
        # Convert to numpy for correlation calculation
        pred_np = predictions.cpu().numpy().flatten()
        label_np = labels.cpu().numpy().flatten()
        
        # Calculate correlation
        correlation_matrix = np.corrcoef(pred_np, label_np)
        correlation = correlation_matrix[0, 1] if not np.isnan(correlation_matrix[0, 1]) else 0.0
        
        return correlation

