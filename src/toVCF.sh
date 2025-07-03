#!/bin/bash
#SBATCH -p barc
#SBATCH -t 4:00:00

dir=$1
mkdir -p vcf

for chr in $(seq 1 22) X; do 
    prefix="ukb_imp_chr${chr}_bct_vars"
    path="${dir}/${prefix}"
    if [ ! -f "${path}.vcf" ]
    then 
        echo "converting ${prefix}..."
        sbatch -p barc --mem=50G -t 2:00:00 --wrap="module load plink; plink2 --bgen ${path}.bgen ref-first --sample ${path}.sample --export vcf --out vcf/${prefix}" 
    else
        echo "vcf file already exists! skipping..."
    fi
done

## the vcf output has many additional columns  past INFO -- FORMAT and many columns that look like "2409086_2409086 3514581_3514581" -- unclear what these are but not needed for snpEff
## see src/cleanVCF.sh to awk these extra columns out; that way the snpEff annotation will be more manageable 