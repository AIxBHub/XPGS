library(bigsnpr)
library(magrittr)
library(s3)
options(bigstatsr.check.parallel.blas = FALSE)


# Read external summary statistics
sumstats <- bigreadr::fread2("data/mcv.assoc")
str(sumstats)

sumstats <- sumstats %>%
  #  dplyr::select(CHR_num, ID, BP, REF, ALT, EFFECT, SE) %>%
  dplyr::rename('chr' = 'CHR_num', 
                'rsid' = 'ID',
                'pos' = 'BP',
                'a0' = 'REF',
                'a1' = 'ALT',
                'beta' = 'EFFECT',
                'beta_se' = 'SE') %>%
  dplyr::mutate(n_eff = 408112)

# Read from bed/bim/fam, it generates .bk and .rds files.
snp_readBed("vcf-GRCh37/")

obj.bigSNP <- snp_attach("../data/ukb_geno_bct_pheno/ukb_cal_merged.rds")

map <- obj.bigSNP$map %>%
  #  dplyr::select(chromosome, marker.ID, physical.pos, allele1, allele2) %>%
  dplyr::mutate(chromosome = as.integer(chromosome)) %>%
  dplyr::rename('chr' = 'chromosome',
                'rsid' = 'marker.ID',
                'pos' = 'physical.pos',
                'a0' = 'allele1',
                'a1' = 'allele2')