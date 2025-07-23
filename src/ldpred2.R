renv::load()
library(bigsnpr)
library(magrittr)
library(s3)
options(bigstatsr.check.parallel.blas = FALSE)

#map_ldref <- readRDS('data/map_hms_plus.rds')
#str(map_ldref)

# Read external summary statistics
sumstats <- bigreadr::fread2("data/mcv.assoc", nThread = 12)
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

#TODO move path to params file
if(!file.exists('PGS001990/ukb_imp_bct_vars_merged.rds')) {
  snp_readBed("PGS001990/ukb_imp_bct_vars_merged.bed")
}

obj.bigSNP <- snp_attach("PGS001990/ukb_imp_bct_vars_merged.rds")
#str(obj.bigSNP)

map <- obj.bigSNP$map %>%
  #  dplyr::select(chromosome, marker.ID, physical.pos, allele1, allele2) %>%
  dplyr::mutate(chromosome = as.integer(chromosome)) %>%
  dplyr::rename('chr' = 'chromosome',
                'rsid' = 'marker.ID',
                'pos' = 'physical.pos',
                'a0' = 'allele1',
                'a1' = 'allele2')

G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
RSID <- obj.bigSNP$map$marker.ID
y   <- obj.bigSNP$fam$affection

set.seed(123)
ind.train <- sample(nrow(G), 350000) # 350k individuals for tuning hyper param
ind.val <- setdiff(rows_along(G), ind.train) # leaves ~150k for testing

#sumstats$n_eff <- sumstats$n_eff
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map)

POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = 12) 

tmp <- tempfile(tmpdir = "tmp-data")
do.call(file.remove, list(list.files('tmp-data','*.sbk', full.names = T)))
for (chr in 1:22) {
  print(chr)
  
  ind.chr <- which(df_beta$chr == chr)
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  
  # Add this sorting step
  ord <- order(POS2[ind.chr2])
  
  # Use the sorted indices
  corr0 <- snp_cor(G, 
                  ind.col = ind.chr2[ord], 
                  infos.pos = POS2[ind.chr2[ord]],
                  size = 3 / 1000, 
                  ncores = 12)
  
  # Also need to keep df_beta in the same sorted order
  if (chr == 1) {
    df_sort <- df_beta[ind.chr[ord], ]
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    df_sort <- rbind(df_sort, df_beta[ind.chr[ord], ])
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

saveRDS(df_sort, 'data/df_beta_sort.RDS')

## is this the correct sample size?? using value i found for sum stats in paper
(ldsc <- with(df_sort, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = 408112, blocks = NULL)))


#ldsc_h2_est <- ldsc[['h2']]
#df_sort$n_eff <- 408112 ## right??
#beta_inf <- snp_ldpred2_inf(corr, df_sort, h2 = ldsc_h2_est)
#
#G_imputed <- snp_fastImputeSimple(G, method = "mean2")
#pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.val, ind.col = df_sort[["_NUM_ID_"]])
#pcor(pred_inf, y[ind.val], NULL)
#
##ind.row <- rows_along(G)
##maf <- snp_MAF(G, ind.row = ind.row, ind.col = df_beta$`_NUM_ID_`, ncores = 12)
##maf_thr <- 1 / sqrt(length(ind.row))  # threshold I like to use
##df_beta <- df_beta[maf > maf_thr, ]
#
##in_test <- vctrs::vec_in(info_snp[, c("chr", "pos")], map[, c("chr", "pos")])
##info_snp <- info_snp[in_test, ]
#
#tmp <- tempfile(tmpdir = "tmp-data")
#
#for (chr in 1:22) {
#  
#  cat(chr, ".. ", sep = "")
#  
#  ## indices in 'df_beta'
#  ind.chr <- which(info_snp$chr == chr)
#  ## indices in 'map_ldref'
#  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
#  ## indices in 'corr_chr'
#  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
#  
#  corr_chr <- readRDS(paste0("data/corr_hm3_plus/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]
#  
#  if (chr == 1) {
#    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
#  } else {
#    corr$add_columns(corr_chr, nrow(corr))
#  }
#}
#
#(ldsc <- with(info_snp, snp_ldsc(ld, ld_size = nrow(map_ldref),
#                                chi2 = (beta / beta_se)^2,
#                                sample_size = nrow(info_snp),
#                                ncores = 12)))
#
#
#