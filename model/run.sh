#! /bin/bash


godag="data/GO/go-basic.obo"
annotation_file="PGS001990/ukb_imp_bct_vars_merged_clean_annotations-NodeNorm.csv"

module load python/3.12

source .venv/bin/activate

python model/go.py -obo $godag -anno $annotation_file -test 5 -out PGS001990
