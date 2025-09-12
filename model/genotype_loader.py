#import glow
#from pyspark.sql import SparkSession
import polars as pl
import re
import torch
import numpy as np

def loader1_spark():
    spark = (
        SparkSession.builder
        .appName("VCF Loader")
        # Tell Spark to pull in the Glow JAR
        .config("spark.jars.packages", "io.projectglow:glow-spark3_2.12:2.0.0")
        # Required so Glow can read/write BGZF-compressed files (VCF bgz/gz)
        .config("spark.hadoop.io.compression.codecs", "io.projectglow.sql.util.BGZFCodec")
        .getOrCreate()
    )

    # Register Glow functions with Spark
    glow.register(spark)

    vcf_file = "PGS001990/ukb_imp_bct_vars_merged_effects.vcf"

    # Example: load VCF
    df = (
        spark.read.format("vcf")
        .option("includeSampleIds", True)
        .load(vcf_file)
    )


    #vcf_reader = vcfpy.Reader.from_path(vcf)

    #samples = vcf_reader.header.samples.names
    ## for testing
    #samples = samples[0:5]
    #print(samples)

    #sample_to_genotypes = { samp: {} for samp in samples }

    #print(sample_to_genotypes)

    #i = 0
    #for record in vcf_reader:
            ## Determine an identifier for the variant
            ## If record.ID is set and not “.”, use that; else make one
        #varid = record.ID
        #print(varid)
        #for call in record.calls[0:5]:
            #samp = call.sample
            ## call.data is an OrderedDict containing FORMAT fields
            ## The genotype is under key “GT”
            #gt = call.data.get('GT')
            #print(samp, gt)

            #if gt == '0/0':
                #gt_score = 0
            #elif gt == '0/1':
                #gt_score = 1
            #elif gt == '1/1':
                #gt_score = 2
            #else:
                #print(f'{gt} ??')

            #sample_to_genotypes[samp][varid[0]] = gt_score

        #if i == 4:
            #print(sample_to_genotypes)
            #exit()
        #i += 1
        
def loader1_polars():
    with open("/proj/jchunglab/projects/genesight/ukb_imp_bct_vars_merged_effects.vcf") as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.lstrip("#").strip().split("\t")
                break
            
    # Load the data with the fixed header
    df = pl.scan_csv(
        "/proj/jchunglab/projects/genesight/ukb_imp_bct_vars_merged_effects.vcf",
        has_header=False,
        new_columns=header,
        separator="\t",
        skip_rows=38
    )
    df = df.rename({col: re.sub(r"_.*", "", col) for col in df.columns})
    columns = df.collect_schema().names()[9:]
    df_100 = df.select(df.columns[9:100])
    
    # Select only data columns
    dfdata = df.select(df.columns[9:])
    # Convert to a numpy array of strings first
    arr = dfdata.collect().to_numpy().astype(str)

    # Replace missing genotypes "./." with "0/0"
    arr[arr == "./."] = "0/0"

    # fixed width string array, 3 bytes each
    arr_bytes = np.ascontiguousarray(arr.astype("S3"))

    # reinterpret memory as raw uint8
    arr_uint8 = arr_bytes.view("uint8").reshape(-1, 3)

    # extract first and third character
    left = arr_uint8[:, 0]
    right = arr_uint8[:, 2]

    # ASCII '0' == 48
    encoded = (left == 48).astype(np.int8) + (right == 48).astype(np.int8)

    encoded = encoded.reshape(arr.shape)
    
    # Store in parquet format for further preprocessing
    dfgeno = pl.from_numpy(encoded)  # or pl.from_numpy(encoded)
    dfgeno = dfgeno.rename({i: col for i, col in enumerate(columns)})
    dfgeno = dfgeno.with_columns(df["ID"].alias("rsID"))
    
    dfgeno.write_parquet("genotypes.parquet", compression="zstd")

    # Store in tensors for training
    #t = torch.from_numpy(encoded).to(torch.int8)  # values -1,0,1,2 fit in int8
    #torch.save(t, "genotypes.pt")


    
if __name__ == "__main__":
    loader1_polars()