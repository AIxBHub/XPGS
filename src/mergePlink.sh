#! /bin/bash
#SBATCH --mem=100G
#SBATCH -p barc
#SBATCH --constraint=rhel9
#SBATCH -n 12
#SBATCH -t 12:00:00

#dir='vcf-GRCh37'
dir='/proj/jchunglab/data/ukb_geno_bct_pheno/bct_vars/ukb_imp_bct_vars_bgen/'

module load plink/2.00 bcftools

cd $dir

## loop the chromsome vcf files and sort before concat
for chr in {1..22} X; do
  bcftools sort ukb_imp_chr${chr}_bct_vars_clean.ann.vcf -Oz -o ukb_imp_chr${chr}_bct_vars_clean.ann.vcf.gz
done

## concat all the sorted vcf files together 
bcftools concat *ann.vcf.gz -o ukb_imp_bct_vars_merge.vcf.gz

## convert vcf to binary
plink2 --vcf ukb_imp_bct_vars_merge.vcf.gz --make-bed --threads 12 --out ukb_imp_bct_vars_merge
