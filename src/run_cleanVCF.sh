#!/bin/bash

## loop chromosomes and remove those unwanted columns in the VCF files
## runs cleanVCF.sh
for chr in $(seq 1 22) X; do 
  prefix="ukb_imp_chr${chr}_bct_vars"
  input_path="vcf/${prefix}.vcf"
  output_path="vcf/${prefix}-clean.vcf"
  echo $input_path
  echo $output_path
  sbatch -p barc --mem=8G -t 1:00:00 --wrap="src/cleanVCF.sh ${input_path} ${output_path}"
done