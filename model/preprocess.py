"""
drafting script to take annotated variant file and generate map that can be sued to build model
"""
import pandas as pd
import sys
#sys.path.append('model')
#import 

annotation_file = 'PGS001990/ukb_imp_bct_vars_merged_clean_annotations-NodeNorm.csv'

afile = pd.read_csv(annotation_file)

