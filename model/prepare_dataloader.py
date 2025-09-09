from pyspark.sql import SparkSession
from pyspark.sql import functions as F
from pyspark.sql.types import StructType,StructField, StringType, DoubleType
from torch.utils.data import Dataset, DataLoader
import torch
import os

def get_data(train_file, test_file, testratio, g2imap, batchsize):
    #Initialize Spark Session
    spark = SparkSession.builder \
        .appName("PySpark PyTorch Batching") \
        .config("spark.driver.memory", "32g").config("spark.executor.memory", "32g") \
        .getOrCreate()

    # Define schema
    schema = StructType() \
        .add("Index1", StringType(), True) \
        .add("Index2", StringType(), True) \
        .add("Value", DoubleType(), True) \
            
    #Convert Spark DataFrame to RDD for PyTorch Processing
    def convert_to_tensors(row):
        features = torch.tensor([g2i(row["Index1"]), g2i(row["Index2"])], dtype=torch.int)
        label = torch.tensor(row['Value'], dtype=torch.double)  # Example label
        return features, label

    # Gene 2 id map
    def g2i(genename):
        return g2imap[genename]
    #Create a PyTorch Dataset
    class PySparkDataset(Dataset):
        def __init__(self, rdd):
            self.data = rdd.repartition(10).collect()  # Collect only the partitioned data

        def __len__(self):
            return len(self.data)

        def __getitem__(self, idx):
            return self.data[idx]
    if not test_file:
        
        df = spark.read.csv(train_file, 
                            header=True,
                            sep = '\t', 
                            schema = schema)
        train_data, test_data = df.randomSplit([(1-testratio), testratio], seed=42)
        #split_name = "split_"+ os.path.basename(train_file).split(".")[0]
    
    else:
        train_data = spark.read.csv(train_file, 
                            header=True,
                            sep = '\t', 
                            schema = schema)
        
        test_data = spark.read.csv(test_file, 
                            header=True,
                            sep = '\t', 
                            schema = schema)
        
        
    #split_name = "split_"+ os.path.basename(train_file).split(".")[0]
        
    rddtrain = train_data.rdd.map(convert_to_tensors)
    #PyTorch DataLoader for Efficient Batching
    trainset = PySparkDataset(rddtrain)
    
    rddtest = test_data.rdd.map(convert_to_tensors)
    testset = PySparkDataset(rddtest)
    
    #if not os.path.exists(f"sample/test_{split_name}"):
    #    test_data.coalesce(1).write.option("delimiter", "\t").csv(f"sample/test_{split_name}")
    #if not os.path.exists(f"sample/train_{split_name}"):
    #    train_data.coalesce(1).write.option("delimiter", "\t").csv(f"sample/train_{split_name}")  
    
            
    return DataLoader(trainset, batch_size=batchsize, shuffle=True, drop_last=False), DataLoader(testset, batch_size=batchsize, shuffle=False)
        

    
