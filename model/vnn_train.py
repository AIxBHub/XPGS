"""
VNN Training Module

This module contains the core training loop for the Visible Neural Network (VNN).
It handles:
- Model initialization and setup
- Training loop with auxiliary losses
- Validation and metric tracking
- Model checkpointing based on validation performance
- Comprehensive logging of training progress

Classes:
    VNNTrainer: Main trainer class for VNN models
"""

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
from ccc_loss import *
import os


class VNNTrainer():
    """
    Trainer class for Visible Neural Network models.

    This class implements the complete training pipeline including:
        - Model initialization with masked connections based on ontology
        - Multi-task learning with auxiliary outputs at each ontology term
        - Training with backpropagation and masked gradient updates
        - Validation with correlation metrics
        - Model checkpointing based on validation performance
        - Comprehensive logging and metric tracking

    Attributes:
        data_wrapper (TrainingDataWrapper): Container with all data and configuration
        model (DCellNN): The hierarchical neural network model
        logger (logging.Logger): Logger for tracking training progress
    """

    def __init__(self, data_wrapper, logger=None):
        """
        Initialize the VNN trainer.

        Args:
            data_wrapper (TrainingDataWrapper): Data and configuration container
            logger (logging.Logger, optional): Logger instance for tracking progress
        """
        self.data_wrapper = data_wrapper
        self.logger = logger

        if self.logger:
            self.logger.info("Initializing VNN Trainer")


    def train_model(self):
        """
        Main training loop for the VNN model.

        This method implements the complete training pipeline:
            1. Initialize the model architecture based on ontology
            2. Set up masked connections for gene-term relationships
            3. Create data loaders for training and validation
            4. Configure optimizer (AdamW)
            5. Training loop:
                - Forward pass with auxiliary outputs
                - Multi-task loss calculation
                - Backward pass with masked gradients
                - Validation
                - Metric tracking and logging
                - Model checkpointing
            6. Save final model and metrics

        The training uses:
            - MSE loss at final output and all intermediate subsystem outputs
            - Auxiliary loss weight (alpha) to balance final vs intermediate predictions
            - Masked gradients to respect ontology structure
            - Pearson correlation for evaluation
            - Model saving based on validation correlation improvement

        Returns:
            float: Final minimum validation loss (or None if not tracked)
        """
        if self.logger:
            self.logger.info("="*60)
            self.logger.info("BUILDING VNN MODEL")
            self.logger.info("="*60)

        # Initialize the model
        self.model = DCellNN(self.data_wrapper)
        self.model.to(self.data_wrapper.cuda)

        if self.logger:
            self.logger.info(f"Model architecture created")
            self.logger.info(f"  Total parameters: {sum(p.numel() for p in self.model.parameters()):,}")
            self.logger.info(f"  Trainable parameters: {sum(p.numel() for p in self.model.parameters() if p.requires_grad):,}")
            self.logger.info(f"  Number of hierarchy layers: {len(self.model.term_layer_list)}")
            self.logger.info(f"  Root term: {self.model.root}")
            self.logger.info(f"  Input dimension (genes): {self.model.gene_dim}")

        # Initialize tracking variables
        min_loss = None
        max_corr = None

        # Create term masks for enforcing ontology structure
        # Masks ensure that each term only receives input from its annotated genes
        if self.logger:
            self.logger.debug("Creating term masks for ontology structure...")

        term_mask_map = utils.create_term_mask(
            self.model.term_direct_gene_map,
            self.model.gene_dim,
            self.data_wrapper.cuda
        )

        # Initialize gene-to-term connection weights with masks
        # This ensures that each term only connects to its relevant genes
        for name, param in self.model.named_parameters():
            term_name = name.split('_')[0]
            if '_direct_gene_layer.weight' in name:
                # Apply mask to gene-to-term connections
                param.data = torch.mul(param.data, term_mask_map[term_name])
            else:
                # Other weights are initialized normally
                param.data = param.data

        if self.logger:
            self.logger.info("Term masks applied to enforce ontology structure")

        # Load data
        if self.logger:
            self.logger.info("\nLoading training and validation data...")

        train_loader, val_loader = prepare_dataloader.get_from_pt(
            self.data_wrapper.train,
            self.data_wrapper.test,
            self.data_wrapper.batchsize
        )

        if self.logger:
            self.logger.info(f"  Training batches: {len(train_loader)}")
            self.logger.info(f"  Validation batches: {len(val_loader)}")
            self.logger.info(f"  Batch size: {self.data_wrapper.batchsize}")

        # Set up optimizer
        optimizer = torch.optim.AdamW(
            self.model.parameters(),
            lr=self.data_wrapper.lr,
            betas=(0.9, 0.99),
            eps=1e-05,
            weight_decay=self.data_wrapper.wd
        )
        optimizer.zero_grad()

        if self.logger:
            self.logger.info(f"\nOptimizer: AdamW")
            self.logger.info(f"  Learning rate: {self.data_wrapper.lr}")
            self.logger.info(f"  Weight decay: {self.data_wrapper.wd}")
            self.logger.info(f"  Betas: (0.9, 0.99)")

        # Open log file for metrics
        logout = open(self.data_wrapper.outfile, "w")
        logout.write("epoch\ttrain_corr\ttrain_loss\ttrue_auc\tpred_auc\tval_corr\tval_loss\telapsed_time\n")

        if self.logger:
            self.logger.info(f"\nMetrics will be logged to: {self.data_wrapper.outfile}")

        if self.logger:
            self.logger.info("\n" + "="*60)
            self.logger.info("STARTING TRAINING")
            self.logger.info("="*60)
            self.logger.info(f"Total epochs: {self.data_wrapper.epochs}")
            self.logger.info(f"Alpha (auxiliary loss weight): {self.data_wrapper.alpha}")
            self.logger.info(f"Delta (improvement threshold): {self.data_wrapper.delta}")
            self.logger.info("="*60 + "\n")

        epoch_start_time = time.time()

        # ====================================================================
        # MAIN TRAINING LOOP
        # ====================================================================
        for epoch in range(self.data_wrapper.epochs):
            if self.logger:
                self.logger.info(f"Epoch {epoch+1}/{self.data_wrapper.epochs}")

            # ================================================================
            # TRAINING PHASE
            # ================================================================
            self.model.train()
            train_predict = torch.zeros(0, 0).to(self.data_wrapper.cuda)
            train_label_gpu = None

            for i, (inputdata, labels) in enumerate(train_loader):
                # Prepare input data
                features = inputdata.to(torch.float32).to(self.data_wrapper.cuda)
                labels = labels.to(torch.float32).unsqueeze(1).to(self.data_wrapper.cuda)

                # Forward pass
                optimizer.zero_grad()
                cur_out = self.model(features)
                aux_out_map, _ = cur_out

                # Accumulate predictions for correlation calculation
                if train_predict.size()[0] == 0:
                    train_predict = aux_out_map['final'].data
                    train_label_gpu = labels
                else:
                    train_predict = torch.cat([train_predict, aux_out_map['final'].data], dim=0)
                    train_label_gpu = torch.cat([train_label_gpu, labels], dim=0)

                # Calculate multi-task loss
                # Final output gets full loss, auxiliary outputs get weighted loss
                total_loss = 0
                loss_fn = nn.MSELoss()

                for name, output in aux_out_map.items():
                    if name == 'final':
                        # Final output loss (full weight)
                        total_loss += loss_fn(output, labels)
                    else:
                        # Auxiliary output loss (weighted by alpha)
                        total_loss += self.data_wrapper.alpha * loss_fn(output, labels)

                # Backward pass
                total_loss.backward()

                # Apply masks to gradients to maintain ontology structure
                for name, param in self.model.named_parameters():
                    if '_direct_gene_layer.weight' not in name:
                        continue
                    term_name = name.split('_')[0]
                    if param.grad is not None:
                        param.grad.data = torch.mul(param.grad.data, term_mask_map[term_name])

                # Update weights
                optimizer.step()

            # Calculate training metrics
            train_corr = utils.pearson_corr(train_predict, train_label_gpu)

            # ================================================================
            # VALIDATION PHASE
            # ================================================================
            self.model.eval()
            val_predict = torch.zeros(0, 0).to(self.data_wrapper.cuda)
            val_label_gpu = None
            val_loss = 0

            with torch.no_grad():
                for i, (inputdata, labels) in enumerate(val_loader):
                    # Prepare validation data
                    vfeatures = inputdata.to(torch.float32).to(self.data_wrapper.cuda)
                    labels = labels.to(torch.float32).unsqueeze(1).to(self.data_wrapper.cuda)

                    # Forward pass
                    aux_out_map, hidden_embeddings_map = self.model(vfeatures)

                    # Accumulate predictions
                    if val_predict.size()[0] == 0:
                        val_predict = aux_out_map['final'].detach().cpu()
                        val_label_gpu = labels.cpu()
                    else:
                        val_predict = torch.cat([val_predict, aux_out_map['final'].detach().cpu()], dim=0)
                        val_label_gpu = torch.cat([val_label_gpu, labels.cpu()], dim=0)

                    # Calculate validation loss (final output only)
                    loss_fn = nn.MSELoss()
                    val_loss += loss_fn(aux_out_map['final'], labels)

            # Calculate validation correlation
            val_corr = utils.pearson_corr(val_predict, val_label_gpu)

            # ================================================================
            # EPOCH LOGGING AND CHECKPOINTING
            # ================================================================
            epoch_end_time = time.time()
            epoch_duration = epoch_end_time - epoch_start_time

            # Calculate additional metrics
            true_auc = torch.mean(train_label_gpu)
            pred_auc = torch.mean(train_predict)

            # Log to console
            if self.logger:
                self.logger.info(
                    f"  Train Loss: {total_loss:.4f} | "
                    f"Train Corr: {train_corr:.4f} | "
                    f"Val Loss: {val_loss:.4f} | "
                    f"Val Corr: {val_corr:.4f} | "
                    f"Time: {epoch_duration:.1f}s"
                )

            # Log to file
            logout.write(
                f"{epoch}\t{train_corr:.4f}\t{total_loss:.4f}\t"
                f"{true_auc:.4f}\t{pred_auc:.4f}\t{val_corr:.4f}\t"
                f"{val_loss:.4f}\t{epoch_duration:.4f}\n"
            )
            logout.flush()

            # Model checkpointing based on validation correlation
            pt_path = os.path.join(
                self.data_wrapper.modeldir,
                f'{self.data_wrapper.strfime}_model.pt'
            )

            # Save model if validation correlation improves
            if max_corr is None and val_corr > 0.1:
                max_corr = val_corr
                torch.save(self.model, pt_path)
                if self.logger:
                    self.logger.info(f"  *** Model saved (initial, val_corr={val_corr:.4f})")
            elif max_corr and (val_corr - max_corr) > self.data_wrapper.delta:
                improvement = val_corr - max_corr
                max_corr = val_corr
                torch.save(self.model, pt_path)
                if self.logger:
                    self.logger.info(
                        f"  *** Model saved (improved by {improvement:.4f}, "
                        f"new best val_corr={val_corr:.4f})"
                    )

            # Update epoch start time for next iteration
            epoch_start_time = epoch_end_time

        # ====================================================================
        # TRAINING COMPLETE
        # ====================================================================
        logout.close()

        if self.logger:
            self.logger.info("\n" + "="*60)
            self.logger.info("TRAINING COMPLETE")
            self.logger.info("="*60)
            if max_corr:
                self.logger.info(f"Best validation correlation: {max_corr:.4f}")
            self.logger.info(f"Final model saved to: {pt_path}")
            self.logger.info(f"Metrics saved to: {self.data_wrapper.outfile}")
            self.logger.info("="*60 + "\n")

        # Future work: Calculate RLIPP scores
        if self.logger:
            self.logger.debug("RLIPP calculation not yet implemented")
            self.logger.debug("Use calculate_rlipp.py script to compute RLIPP scores from saved model")

        return min_loss
