#import glow
#from pyspark.sql import SparkSession
import polars as pl
import re
import torch
import numpy as np

def map_variants(genotypes, annotations):
    geno = pl.scan_parquet(genotypes)
    anno = pl.read_csv(annotations) 

    print(geno.head().collect())

    variant_map = {}
    for row in anno.iter_rows(named = True):
        #if row['id_x'] not in geno.columns:
        #    continue

        if row['biotype'] != 'protein_coding':
            continue

        gene = row['name'] if row['name'] is not None else row['gene']

        if 'Hsap' in gene:
            gene = gene.split(' ')[0]

        if gene not in variant_map:
            variant_map[gene] = []

        if row['id_x'] not in variant_map[gene]:
            variant_map[gene].append(row['id_x'])

    #df = pl.DataFrame(schema = [(gene, pl.Int8) for gene in genes])
    nsamples = geno.collect().height

    arrs = []
    mapped_genes = []
    for gene, variants in variant_map.items():
        arr = np.zeros((nsamples,1), np.int8) ## empty array for n samples
        for var in variants: ## for each variant aggregate scores
            if var in geno.columns:
                arr += geno.select(var).collect().to_numpy().reshape(nsamples, 1)
        if np.sum(arr) > 0: 
            arrs.append(arr)
            mapped_genes.append(gene)

    arrs = np.concatenate(arrs, axis = 1)
    arrs = arrs.T
    print(arrs.shape)
    t = torch.from_numpy(arrs).to(torch.int8)  # values -1,0,1,2 fit in int8
    print(t)
    print(len(mapped_genes))

if __name__ == "__main__":
    genotypes = "PGS001990/genotypes.parquet"
    annotations = "PGS001990/ukb_imp_bct_vars_merged_clean_annotations-NodeNorm.csv"
    map_variants(genotypes, annotations)

