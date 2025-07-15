# git clone https://github.com/RobokopU24/ORION.git
# python 3.11.11 used
# pip install -r ORION/requirements.txt
import sys
sys.path.append("/Users/jchung/Documents/GRANT/2023_11_RFA_RM_23_003_CommonFundUtility/gits/ORION")
from Common import normalization

import pandas as pd

class respo():
    def __init__(self):
        #self.resp = resp
        self.nn = normalization.NodeNormalizer(conflate_node_types=False)
        self.nn.node_norm_endpoint="https://nodenormalization-sri.renci.org/" # This endpoint might not be permanent. 
    
    def in_druglist(self, nodeid):
        nodes = []
        nodes.append({'id': nodeid})
        self.nn.normalize_node_data(nodes)
        if len(nodes)>0:
            return nodes[0]["equivalent_identifiers"] #any(dfdrug["single_ID"].isin(nodes[0]["equivalent_identifiers"]))
        else:
            return [nodeid]

resp = respo()

# Read variant table, split rsID with multiple EMSEMBL matches, add curie prefix to each ENSEMBL ID
df = pd.read_csv("BCT_related_variants_PGS_models/bct_variant_annotated.csv")
df["gene"] = df["gene"].str.split("-")
df = df.explode("gene").reset_index(drop=True)
df["curie"] = [f"ENSEMBL:{x}" for x in df["gene"]]

# Prepare list of dictionaries for NodeNorm function to run nn. Expand the output again just for merge purpose.
gl = df["gene"].unique()
dfgene = [{'id':f"ENSEMBL:{x}"} for x in gl]
resp.nn.normalize_node_data(dfgene)
dfgene = pd.DataFrame(dfgene)
dfgene = dfgene.explode("equivalent_identifiers").reset_index(drop=True)

merged = pd.merge(df, dfgene[["id", "equivalent_identifiers", "name"]], left_on="curie", right_on="equivalent_identifiers", how="left")
merged.to_csv("/Users/jchung/Documents/GRANT/2023_11_RFA_RM_23_003_CommonFundUtility/Data/blood_cell_trait_ukbiobank/BCT_related_variants_PGS_models/bct_variant_annotated_nodenormed.csv")