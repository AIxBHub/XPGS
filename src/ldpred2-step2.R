library(ggplot2)
library(bigsnpr)
library(magrittr)
library(s3)
set.seed(123)

options(bigstatsr.check.parallel.blas = FALSE)

load('data/df_beta_sort.Rdata')

(ldsc <- with(df_sort, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = 408112, blocks = NULL)))
ldsc_h2_est <- ldsc[['h2']]

G_imputed <- snp_fastImputeSimple(G, method = "mean2")
ind.val <- sample(nrow(G_imputed), 350000)
ind.test <- setdiff(rows_along(G_imputed), ind.val)
##################### Grid

(h2_seq <- round(ldsc_h2_est * c(0.1,0.2, 0.3,0.4,0.5,0.6,0.7, 1, 1.4), 4))
(p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
(params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))
# takes less than 2 min with 4 cores
beta_grid <- snp_ldpred2_grid(corr, df_sort, params, ncores = 12)

#beta_inf <- snp_ldpred2_inf(corr, df_sort, h2 = ldsc_h2_est)
#
#pred_inf <- big_prodVec(G_imputed, beta_inf, ind.row = ind.val, ind.col = df_sort[["_NUM_ID_"]])
#pcor(pred_inf, y[ind.val], NULL)
pred_grid <- big_prodMat(G_imputed, beta_grid, ind.col = df_sort[["_NUM_ID_"]])
params$score <- apply(pred_grid[ind.val, ], 2, function(x) {
  if (all(is.na(x))) return(NA)
  summary(lm(y[ind.val] ~ x))$coef["x", 3]
  # summary(glm(y[ind.val] ~ x, family = "binomial"))$coef["x", 3]  # for a binary phenotype
})

ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
  theme_bigstatsr() +
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
  facet_wrap(~ sparse, labeller = label_both) +
  labs(y = "GLM Z-Score", color = "h2") +
  theme(legend.position = "top", panel.spacing = unit(1, "lines"))
####################### Auto
coef_shrink <- 0.95  # reduce this up to 0.4 if you have some (large) mismatch with the LD ref

set.seed(123)  # to get the same result every time
# takes less than 2 min with 4 cores
multi_auto <- snp_ldpred2_auto(
  corr,
  df_sort, 
  h2_init = ldsc_h2_est, 
  vec_p_init = seq_log(1e-4, 0.2, length.out = 30), 
  ncores = 12,
  use_MLE = FALSE,  # uncomment if you have convergence issues or when GWAS power is low (need v1.11.9)
  allow_jump_sign = FALSE, 
  shrink_corr = coef_shrink)

#str(multi_auto, max.level = 1)
str(multi_auto[[1]], max.level = 1)
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
