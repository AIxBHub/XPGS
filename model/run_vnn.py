"""
VNN (Visible Neural Network) Training Script

This script trains a hierarchical neural network based on biological ontologies (e.g., Gene Ontology).
The architecture respects the parent-child relationships between biological terms/subsystems,
creating an interpretable deep learning model for genotype-phenotype prediction.

Based on the DCell model from:
Ma et al. (2018). "Using deep learning to model the hierarchical structure and function of a cell."
Nature Methods, 15(4), 290-298. PMID: 29505029

Usage:
    python run_vnn.py -onto ontology.txt -gene2id gene2ind.txt -train train.pt -test test.pt

Key Features:
    - Hierarchical neural network guided by biological ontology
    - Interpretable subsystem-level predictions
    - Multi-task learning with auxiliary losses at each subsystem
    - Comprehensive logging of training progress and configurations

Author: [Your name]
Date: 2024
"""

import argparse
import copy
import os
import logging
import json
import sys
from datetime import datetime
from data_wrapper import TrainingDataWrapper
from vnn_train import VNNTrainer


def setup_logging(modeldir, timestamp):
    """
    Set up comprehensive logging system for the training process.

    Creates both file and console handlers to log:
    - Configuration parameters
    - Training progress
    - Model architecture details
    - Validation metrics

    Args:
        modeldir (str): Directory where logs will be saved
        timestamp (str): Timestamp string for unique log file naming

    Returns:
        logging.Logger: Configured logger instance
    """
    # Create logs directory if it doesn't exist
    log_dir = os.path.join(modeldir, 'logs')
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    # Create logger
    logger = logging.getLogger('VNN_Training')
    logger.setLevel(logging.DEBUG)

    # Remove any existing handlers to avoid duplicates
    logger.handlers = []

    # File handler - detailed logging
    log_file = os.path.join(log_dir, f'{timestamp}_training.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler.setFormatter(file_formatter)

    # Console handler - important messages only
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(levelname)s: %(message)s')
    console_handler.setFormatter(console_formatter)

    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


def log_configuration(logger, args):
    """
    Log all configuration parameters to both file and console.

    This creates a permanent record of training configuration for reproducibility.

    Args:
        logger (logging.Logger): Logger instance
        args (argparse.Namespace): Parsed command-line arguments
    """
    logger.info("="*80)
    logger.info("VNN TRAINING CONFIGURATION")
    logger.info("="*80)

    # Data files
    logger.info("\nData Configuration:")
    logger.info(f"  Ontology file:        {args.onto}")
    logger.info(f"  Gene2ID mapping:      {args.gene2id}")
    logger.info(f"  Training data:        {args.train}")
    logger.info(f"  Test data:            {args.test if args.test else 'None (will split from training)'}")
    logger.info(f"  Test set ratio:       {args.testsetratio}")

    # Model architecture
    logger.info("\nModel Architecture:")
    logger.info(f"  Hidden neurons/term:  {args.genotype_hiddens}")
    logger.info(f"  Min dropout layer:    {args.min_dropout_layer}")
    logger.info(f"  Dropout fraction:     {args.dropout_fraction}")

    # Training hyperparameters
    logger.info("\nTraining Hyperparameters:")
    logger.info(f"  Learning rate:        {args.lr}")
    logger.info(f"  Weight decay:         {args.wd}")
    logger.info(f"  Alpha (aux loss):     {args.alpha}")
    logger.info(f"  Batch size:           {args.batchsize}")
    logger.info(f"  Epochs:               {args.epoch}")
    logger.info(f"  Delta (improvement):  {args.delta}")

    # Output configuration
    logger.info("\nOutput Configuration:")
    logger.info(f"  Model directory:      {args.modeldir}")
    logger.info(f"  Metrics output:       {args.metric_output}")
    logger.info(f"  GPU device:           {args.cuda}")
    logger.info(f"  Optimization mode:    {args.optimize}")

    logger.info("="*80 + "\n")

    # Also save configuration as JSON for easy parsing
    config_dict = vars(args)
    config_file = os.path.join(args.modeldir, 'logs', f'{datetime.now().strftime("%Y%m%d%H%M")}_config.json')
    with open(config_file, 'w') as f:
        json.dump(config_dict, f, indent=4)
    logger.info(f"Configuration saved to: {config_file}\n")

from bb_train import ANNTrainer

def main():
    """
    Main training function.

    Workflow:
        1. Parse command-line arguments
        2. Set up logging system
        3. Log all configurations
        4. Load and prepare data using TrainingDataWrapper
        5. Initialize and train VNN model
        6. Save trained model and metrics

    The training process includes:
        - Hierarchical neural network construction based on ontology
        - Multi-task learning with auxiliary outputs at each subsystem level
        - Validation at each epoch
        - Model checkpointing based on validation performance
    """
    # Create argument parser with detailed descriptions
    parser = argparse.ArgumentParser(
        description='Train Visible Neural Network (VNN) for genotype-phenotype prediction',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # ========== Required Arguments ==========
    parser.add_argument('-onto',
                       help='Ontology file used to guide neural network structure (e.g., GO terms)',
                       type=str,
                       required=True)
    parser.add_argument('-gene2id',
                       help='Gene to ID mapping file (tab-separated: gene_name gene_id)',
                       type=str,
                       required=True)
    parser.add_argument('-train',
                       help='Training dataset (.pt file with X and y tensors)',
                       type=str,
                       required=True)

    # ========== Data Arguments ==========
    parser.add_argument('-test',
                       help='Test dataset (.pt file). If not provided, will split from training data',
                       type=str,
                       default=None)
    parser.add_argument('-testsetratio',
                       help='Test set ratio (0-1) if no separate test set provided',
                       type=float,
                       default=0.2)

    # ========== Model Architecture Arguments ==========
    parser.add_argument('-genotype_hiddens',
                       help='Number of hidden neurons for each term/subsystem in the hierarchy',
                       type=int,
                       default=50)
    parser.add_argument('-min_dropout_layer',
                       help='Layer index to start applying dropout (0=all layers, 2=skip first 2 layers)',
                       type=int,
                       default=2)
    parser.add_argument('-dropout_fraction',
                       help='Dropout probability (0-1)',
                       type=float,
                       default=0.1)

    # ========== Training Hyperparameters ==========
    parser.add_argument('-lr',
                       help='Learning rate for AdamW optimizer',
                       type=float,
                       default=0.0002)
    parser.add_argument('-wd',
                       help='Weight decay for L2 regularization',
                       type=float,
                       default=0.00005)
    parser.add_argument('-alpha',
                       help='Weight for auxiliary losses (0-1). Final loss = main_loss + alpha * aux_losses',
                       type=float,
                       default=0.1)
    parser.add_argument('-batchsize',
                       help='Batch size for training and validation',
                       type=int,
                       default=64)
    parser.add_argument('-epoch',
                       help='Number of training epochs',
                       type=int,
                       default=50)
    parser.add_argument('-delta',
                       help='Minimum improvement in validation correlation to save model',
                       type=float,
                       default=0.01)

    # ========== Output Arguments ==========
    parser.add_argument('-modeldir',
                       help='Directory to save trained models and logs',
                       type=str,
                       default='MODEL/')
    parser.add_argument('-metric_output',
                       help='Filename for epoch-by-epoch metrics (saved in modeldir)',
                       type=str,
                       default="metrics_output.tsv")

    # ========== System Arguments ==========
    parser.add_argument('-cuda',
                       help='GPU device ID (0, 1, 2, ...). Will use CPU if CUDA unavailable',
                       type=int,
                       default=0)
    parser.add_argument('-optimize',
                       help='Optimization mode: 1=Standard training, 2=Hyperparameter optimization (not implemented)',
                       type=int,
                       default=1)

    # Parse arguments
    opt = parser.parse_args()

    # Create model directory if it doesn't exist
    if not os.path.exists(opt.modeldir):
        os.makedirs(opt.modeldir)
        print(f"Created model directory: {opt.modeldir}")

    # Set up logging
    timestamp = datetime.now().strftime('%Y%m%d%H%M')
    logger = setup_logging(opt.modeldir, timestamp)

    logger.info("Starting VNN Training Pipeline")
    logger.info(f"Timestamp: {timestamp}\n")

    # Log all configuration parameters
    log_configuration(logger, opt)

    if not opt.black_box: 
        VNNTrainer(data_wrapper).train_model()
    else:
        ANNTrainer(data_wrapper).train_model()
    #if opt.optimize == 1:
    #    VNNTrainer(data_wrapper).train_model()

    # Initialize data wrapper (loads ontology, genes, and datasets)
    logger.info("Initializing data wrapper...")
    try:
        data_wrapper = TrainingDataWrapper(opt, logger=logger)
        logger.info("Data wrapper initialized successfully\n")
    except Exception as e:
        logger.error(f"Failed to initialize data wrapper: {e}")
        raise

    # Train the model
    logger.info("Starting model training...")
    try:
        trainer = VNNTrainer(data_wrapper, logger=logger)
        trainer.train_model()
        logger.info("Training completed successfully!")
    except Exception as e:
        logger.error(f"Training failed: {e}")
        raise

    logger.info("\n" + "="*80)
    logger.info("VNN TRAINING COMPLETE")
    logger.info("="*80)
    logger.info(f"Model saved in: {opt.modeldir}")
    logger.info(f"Logs saved in: {os.path.join(opt.modeldir, 'logs')}")
    logger.info(f"Metrics saved in: {os.path.join(opt.modeldir, opt.metric_output)}")
    logger.info("="*80)


if __name__ == "__main__":
    main()
