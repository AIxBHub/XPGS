##### 
#df_beta <- readRDS('PGS001990/df_beta-filtered.rds')

POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES) 

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
                   ncores = NCORES)
  
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

save(df_sort,ld, corr, file = 'data/df_beta_sort.Rdata')
