import polars as pl
import re
import torch
import numpy as np

#TODO add test for duplicate family ids before droping subfamily id -- look for duplicates in family ids        
def loader1_polars(vcf_file, outdir):
    #with open("../../../ukb_for_pgs/ukb_imp_bct_vars_merged_effects.vcf") as f:
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.lstrip("#").strip().split("\t")
                break
            
    # Load the data with the fixed header
    df = pl.scan_csv(
        #"../../../ukb_for_pgs/ukb_imp_bct_vars_merged_effects.vcf",
        vcf_file,
        has_header=False,
        new_columns=header,
        separator="\t",
        skip_rows=38
    )
    df = df.rename({col: re.sub(r"_.*", "", col) for col in df.columns})
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

    ## transpose so variants are columns
    encoded_transposed = encoded.T
    dfgeno = pl.DataFrame(encoded_transposed)
    ## set column names to variant ids for genotype mapping
    dfgeno.columns = df.select('ID').collect().to_series().to_list()
    # Store in parquet format for further preprocessing
    dfgeno.write_parquet(f"{outdir}/genotypes.parquet", compression="zstd")

    # Store in tensors for training
    t = torch.from_numpy(encoded).to(torch.int8)  # values -1,0,1,2 fit in int8
    torch.save(t, f"{outdir}/genotypes.pt")


if __name__ == "__main__":

    vcf_file = "PGS001990/ukb_imp_bct_vars_merged_effects.vcf"
    loader1_polars(vcf_file, outdir = "PGS001990")
#    map_variant(annotation_file)