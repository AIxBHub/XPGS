#import glow
#from pyspark.sql import SparkSession
import polars as pl
import re
import torch
import numpy as np

def map_variants(genotypes, annotations, outdir):
    geno = pl.scan_parquet(genotypes).collect()

    ## get variant ids col
    rs_ids = geno.select('rsID').to_series().to_list()
    ## get genotypes 
    data = geno.drop('rsID')
    ## transpose to (samples, variants)
    geno = data.transpose()
    ## Set column names to rs_ids
    geno.columns = rs_ids

    print("Transposed Data Shape:", geno.shape)

    anno = pl.read_csv(annotations) 

    variant_map = {}
    ## build a variant map that has {gene_name : [variants]}
    for row in anno.iter_rows(named = True):

        if row['biotype'] != 'protein_coding':
            continue

        gene = row['name'] if row['name'] is not None else row['gene']

        if 'Hsap' in gene:
            gene = gene.split(' ')[0]

        if gene not in variant_map:
            variant_map[gene] = []

        if row['id_x'] not in variant_map[gene]:
            variant_map[gene].append(row['id_x'])

    ## use the variant map to aggregate genotypes to genes for each sample
    nsamples = geno.height ## number of individual genotypes
    print(f"{nsamples} samples")
    arrs = []
    mapped_genes = []
    for gene, variants in variant_map.items():
        arr = np.zeros((nsamples,1), np.int8) ## empty array for n samples
        for var in variants: ## for each variant aggregate scores
            if var in geno.columns:
                arr += geno.select(var).to_numpy().reshape(nsamples, 1)
        ## only keep arrays where variants in the annotation file (which has unfiltered variants) where matched to variants in the processed genotype file
        ## append these aggregated gene-level arrays to list of arrays
        if np.sum(arr) > 0: 
            arrs.append(arr)
            mapped_genes.append(gene)
    ## concat the array list by x axis -- (samples, genes)
    arrs = np.concatenate(arrs, axis = 1)
    print(f"{arrs.shape} before transpose.")
    ## transpose to (genes, samples)
    arrs = arrs.T
    print(f"{arrs.shape} after transpose.")

    t = torch.from_numpy(arrs).to(torch.int8)
    torch.save(t, f"{outdir}/genotypes_mapped.pt")
    print(t)
    
    print(f"{mapped_genes} mapped genes")
    with open(f'{outdir}/mapped_genes.txt', 'w') as f:
        for line in mapped_genes:
            f.write(f"{line}\n")

if __name__ == "__main__":
    genotypes = "PGS001990/genotypes.parquet"
    annotations = "PGS001990/ukb_imp_bct_vars_merged_clean_annotations-NodeNorm.csv"
    outdir = "PGS001990"
    map_variants(genotypes, annotations, outdir)

