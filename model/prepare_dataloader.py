"""
Data Loading Module for VNN Training

This module provides functions to create PyTorch DataLoaders from different data sources:
- PySpark DataFrames (get_data) - for large-scale distributed data
- PyTorch tensor files (.pt files) (get_from_pt) - for preprocessed data

The module supports:
- Loading training and test datasets
- Automatic train/test splitting
- Batch creation for efficient training

Functions:
    get_data: Load data from CSV files using PySpark (commented out, requires Spark)
    get_from_pt: Load data from PyTorch .pt tensor files (currently used)
"""

from torch.utils.data import Dataset, DataLoader, TensorDataset
import torch
import os


def get_data(train_file, test_file, testratio, g2imap, batchsize):
    """
    Load data from CSV files using PySpark for distributed processing.

    NOTE: This function is currently not used. The get_from_pt() function
    is used instead for loading preprocessed .pt files.

    This function uses PySpark to:
        1. Load genotype-phenotype data from CSV files
        2. Split into train/test if no test file provided
        3. Convert to PyTorch tensors
        4. Create DataLoaders for batching

    Args:
        train_file (str): Path to training CSV file
        test_file (str): Path to test CSV file (or None to auto-split)
        testratio (float): Proportion of data to use for testing if no test_file
        g2imap (dict): Gene name to ID mapping
        batchsize (int): Batch size for DataLoader

    Returns:
        tuple: (train_loader, test_loader)
            Both are torch.utils.data.DataLoader objects

    CSV Format Expected:
        Index1  Index2  Value
        gene1   gene2   0.5
        gene3   gene4   0.8

    Note: Requires PySpark configuration with sufficient memory.
    """
    # NOTE: PySpark code is commented out but kept for reference
    # The current implementation uses get_from_pt() instead

    # Initialize Spark Session
    # spark = SparkSession.builder \
    #     .appName("PySpark PyTorch Batching") \
    #     .config("spark.driver.memory", "32g").config("spark.executor.memory", "32g") \
    #     .getOrCreate()

    # Define schema for CSV data
    # schema = StructType() \
    #     .add("Index1", StringType(), True) \
    #     .add("Index2", StringType(), True) \
    #     .add("Value", DoubleType(), True)

    # Convert Spark DataFrame to RDD for PyTorch Processing
    # def convert_to_tensors(row):
    #     features = torch.tensor([g2i(row["Index1"]), g2i(row["Index2"])], dtype=torch.int)
    #     label = torch.tensor(row['Value'], dtype=torch.double)
    #     return features, label

    # def g2i(genename):
    #     return g2imap[genename]

    # Create a PyTorch Dataset from Spark RDD
    # class PySparkDataset(Dataset):
    #     def __init__(self, rdd):
    #         self.data = rdd.repartition(10).collect()
    #
    #     def __len__(self):
    #         return len(self.data)
    #
    #     def __getitem__(self, idx):
    #         return self.data[idx]

    # Load and optionally split data
    # if not test_file:
    #     df = spark.read.csv(train_file, header=True, sep='\t', schema=schema)
    #     train_data, test_data = df.randomSplit([(1-testratio), testratio], seed=42)
    # else:
    #     train_data = spark.read.csv(train_file, header=True, sep='\t', schema=schema)
    #     test_data = spark.read.csv(test_file, header=True, sep='\t', schema=schema)

    # Convert to PyTorch datasets and create DataLoaders
    # rddtrain = train_data.rdd.map(convert_to_tensors)
    # trainset = PySparkDataset(rddtrain)
    #
    # rddtest = test_data.rdd.map(convert_to_tensors)
    # testset = PySparkDataset(rddtest)
    #
    # return DataLoader(trainset, batch_size=batchsize, shuffle=True, drop_last=False), \
    #        DataLoader(testset, batch_size=batchsize, shuffle=False)

    pass  # Not currently implemented


def get_from_pt(train_file, test_file, batchsize):
    """
    Load preprocessed data from PyTorch .pt tensor files.

    This is the currently used data loading function. It loads pre-processed
    genotype-phenotype data that has been saved as PyTorch tensors.

    Expected .pt file format:
        A dictionary with keys:
            "X": Tensor of shape (n_samples, n_features) - genotype features
            "y": Tensor of shape (n_samples,) - phenotype labels

    Args:
        train_file (str): Path to training .pt file
        test_file (str): Path to test .pt file
        batchsize (int): Batch size for DataLoader

    Returns:
        tuple: (train_dataloader, test_dataloader)
            - train_dataloader: DataLoader for training set (shuffled)
            - test_dataloader: DataLoader for test set (not shuffled)

    Example:
        >>> train_loader, test_loader = get_from_pt(
        ...     "data/train.pt",
        ...     "data/test.pt",
        ...     batchsize=64
        ... )
        >>> for batch_X, batch_y in train_loader:
        ...     # batch_X shape: (64, n_features)
        ...     # batch_y shape: (64,)
        ...     pass
    """
    # Load .pt files
    # Each file contains a dictionary with "X" and "y" keys
    traindata = torch.load(train_file)
    testdata = torch.load(test_file)

    # Create TensorDatasets
    # TensorDataset pairs features (X) with labels (y)
    train_dataset = TensorDataset(traindata["X"], traindata["y"])
    test_dataset = TensorDataset(testdata["X"], testdata["y"])

    # Create DataLoaders
    # Training data is shuffled for better generalization
    # Test data is not shuffled to maintain consistent evaluation
    train_dataloader = DataLoader(
        train_dataset,
        batch_size=batchsize,
        shuffle=True  # Shuffle training data each epoch
    )
    test_dataloader = DataLoader(
        test_dataset,
        batch_size=batchsize,
        shuffle=False  # Don't shuffle test data
    )

    return train_dataloader, test_dataloader
