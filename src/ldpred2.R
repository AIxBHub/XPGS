renv::load()
library(ggplot2)
library(bigsnpr)
library(magrittr)
library(s3)
options(bigstatsr.check.parallel.blas = FALSE)
neff = 563085
set.seed(1234)
NCORES = 24
#map_ldref <- readRDS('data/map_hms_plus.rds')
#str(map_ldref)

# Read external summary statistics
sumstats <- bigreadr::fread2("data/mcv.assoc", nThread = NCORES)
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
  dplyr::mutate(n_eff = neff)



# Read from bed/bim/fam, it generates .bk and .rds files.

#TODO move path to params file
if(!file.exists('PGS001990/ukb_imp_bct_vars_merged.rds')) {
  snp_readBed2("PGS001990/ukb_imp_bct_vars_merged.bed", ncores = NCORES)
}

obj.bigSNP <- snp_attach("PGS001990/ukb_imp_bct_vars_merged.rds")


## conditional to run imputation -- slow for large genotypes
if(is.null(obj.bigSNP$genotypes_imputed)){
  obj.bigSNP$genotypes_imputed <- snp_fastImputeSimple(obj.bigSNP$genotypes,
                                                       method = "mean2", 
                                                       ncores = NCORES)
  obj.bigSNP <- snp_save(obj.bigSNP)
}

map <- obj.bigSNP$map %>%
  #  dplyr::select(chromosome, marker.ID, physical.pos, allele1, allele2) %>%
  dplyr::mutate(chromosome = as.integer(chromosome)) %>%
  dplyr::rename('chr' = 'chromosome',
                'rsid' = 'marker.ID',
                'pos' = 'physical.pos',
                'a0' = 'allele1',
                'a1' = 'allele2')

#rand_sample <- sample(nrow(obj.bigSNP$genotypes_imputed), 1000, replace = F)
#test = snp_subset(obj.bigSNP, ind.row = rand_sample)

#obj.bigSNP$genotypes <- obj.bigSNP$genotypes_imputed

G   <- obj.bigSNP$genotypes_imputed
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
RSID <- obj.bigSNP$map$marker.ID
y   <- obj.bigSNP$fam$affection

#set.seed(123)
#ind.train <- sample(nrow(G), 350000) # 350k individuals for tuning hyper param
#ind.val <- setdiff(rows_along(G), ind.train) # leaves ~150k for testing

map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map)

#nrow(df_beta[as.double(df_beta$P) <= 1e-3 & df_beta$MA_FREQ >= 0.1,])

#df_beta <- df_beta[as.double(df_beta$P) <= 1e-1 & df_beta$MA_FREQ >= 0.1,]
#ind.val <- sample(nrow(G), 350000)
#ind.test <- setdiff(rows_along(G), ind.val)

#TODO move this to the correct pgs directory
#sd <- runonce::save_run(
#  sqrt(big_colstats(G, ind.val, ncores = NCORES)$var),
#  file = "data/sd.rds"
#)

#SUMMARY <- function(info_snp) {
#  chi2 <- with(info_snp, (beta / beta_se)^2)
#  print(round(mean(chi2, na.rm = TRUE), 2))
#  S <- rep(NA, ncol(G)); S[info_snp$`_NUM_ID_`] <- chi2
#  signif <- pchisq(S, df = 1, lower.tail = FALSE) < 5e-8
#  print(sum(signif, na.rm = TRUE))
#  ind.keep <- snp_clumping(
#    G, infos.chr = map$chr, infos.pos = map$pos, S = S,
#    ind.row = ind.val, thr.r2 = 0.01, size = 10e3, ncores = NCORES,
#    exclude = which(is.na(S) | !signif))
#  print(length(ind.keep))
#}
#
#(df_beta <- tidyr::drop_na(tibble::as_tibble(df_beta)))
#SUMMARY(df_beta)
#
#sd_val <- sd[df_beta$`_NUM_ID_`]
#sd_ss <- with(df_beta, 2 / sqrt(n_eff * beta_se^2))

#is_bad <-
#  sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
#qplot(sd_val, sd_ss, color = is_bad, alpha = I(0.5)) +
#  theme_bigstatsr() +
#  coord_equal() +
#  scale_color_viridis_d(direction = -1) +
#  geom_abline(linetype = 2, color = "red") +
#  labs(x = "Standard deviations in the validation set",
#       y = "Standard deviations derived from the summary statistics",
#       color = "Removed?")
#
#df_beta[is_bad, ] %>%
#  dplyr::arrange(P) %>%
#  head(20) %>%
#  dplyr::mutate(freq2 = big_scale()(G, ind.row = ind.val, ind.col = `_NUM_ID_`)$center / 2) %>%
#  dplyr::select(-`_NUM_ID_.ss`, -`_NUM_ID_`)

#saveRDS(df_beta[!is_bad, ], "PGS001990/df_beta-filtered.rds")

#if(!file.exists('data/df_beta_sort.Rdata')){
#  source('src/compute_LD.R')
#} 

#source('src/ldpred-auto.R')