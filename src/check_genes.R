library(magrittr)
library(GenomicRanges)
library(ChIPpeakAnno)
library(EnsDb.Hsapiens.v86)

anno_ensdb_transcript <- toGRanges(EnsDb.Hsapiens.v86, 
                                   feature = "gene")

vars = read.csv('data/union_rsID.csv')

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")


var.gr <- vars %>%
  dplyr::filter(!is.na(start_GRCh38)) %>%
  dplyr::mutate(seqnames = paste0('chr',chr), start = start_GRCh38, end = start + 1) %>%
  dplyr::select(seqnames, start, end, rsid) %>%
  GRanges()

hg38 = genes(TxDb.Hsapiens.UCSC.hg38.knownGene, single.strand.genes.only = F) %>%
  unlist()

#anno <- ChIPpeakAnno::annotatePeakInBatch(var.gr,AnnotationData = anno_ensdb_transcript, )

var.gr <- var.gr %>%
  data.frame() %>%
  dplyr::mutate(in_gene = ifelse(GRanges(.) %over% hg38, T, F))

var.gr %>%
  dplyr::group_by(in_gene) %>%
  dplyr::count()
  dplyr::summarise(n = ())
  
var.gr <- var.gr %>%
  data.frame() %>%
  dplyr::mutate(in_gene = ifelse(GRanges(.) %over% anno_ensdb_transcript, T, F))
var.gr %>%
  dplyr::group_by(in_gene) %>%
  dplyr::count()
  dplyr::summarise(n = ())
  
 