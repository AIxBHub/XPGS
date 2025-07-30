library(ggplot2)
library(bigsnpr)
library(magrittr)
library(s3)
set.seed(123)
NCORES = 24
options(bigstatsr.check.parallel.blas = FALSE)

load('data/df_beta_sort.Rdata')

(ldsc <- with(df_sort, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = neff, blocks = NULL)))
ldsc_h2_est <- ldsc[['h2']]

#G_imputed <- snp_fastImputeSimple(G, method = "mean2")

####################### Auto
coef_shrink <- 0.95  # reduce this up to 0.4 if you have some (large) mismatch with the LD ref

#set.seed(123)  # to get the same result every time
# takes less than 2 min with 4 cores
multi_auto <- snp_ldpred2_auto(
  corr,
  df_sort, 
  h2_init = ldsc_h2_est, 
  vec_p_init = seq_log(1e-4, 0.2, length.out = 30), 
  ncores = NCORES,
  use_MLE = FALSE,  # uncomment if you have convergence issues or when GWAS power is low (need v1.11.9)
  allow_jump_sign = FALSE, 
  shrink_corr = coef_shrink)

#str(multi_auto, max.level = 1)
#str(multi_auto[[2]], max.level = 1)
auto <- multi_auto[[1]]  # first chain
plot_grid(
  qplot(y = auto$path_p_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv"
)
