from neo4j import GraphDatabase
import pandas as pd
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag

godag = GODag('data/GO/go-basic.obo',
              optional_attrs={'relationship'})

annotation_file = "PGS001990/ukb_imp_bct_vars_merged_clean_annotations-NodeNorm.csv"

effects = ['missense_variant','splice_region_variant','stop_gained','5_prime_UTR_premature_start_codon_gain_variant']
#effects = ['splice_region_variant']
#genes = ["SSU72"]
taxon = "NCBITaxon:9606" # homo sapien tax
driver = GraphDatabase.driver("bolt://robokopkg.renci.org:7687", auth=("neo4j", "password"))

query = """
MATCH (n:`biolink:GeneOrGeneProduct`)-[r]-(d)
WHERE n.name = $gene_name
  AND n.taxon = $taxon
  AND r.primary_knowledge_source = "infores:goa"
  AND "biolink:BiologicalProcess" IN labels(d)
  AND NOT "biolink:CellularComponent" IN labels(d)
RETURN 
    n.id AS gene_id,
    n.name AS gene_name,
    type(r) AS predicate,
    r.primary_knowledge_source AS pks,
    d.id AS GO_id,
    d.name AS GO_name
"""

afile = pd.read_csv(annotation_file)
genes = afile[afile['ann'].str.contains('|'.join(effects))]['name'].sort_values().unique() ## find all partial matches to effects and grab all unique gene names

all_results = []

## for testing
#genes = genes[0:5]
gene_map = {}
with driver.session() as session:
    for gene in genes:
        result = session.run(query, {"gene_name": gene, "taxon": taxon})
        data = [r.data() for r in result]
#        print(data[])
#        exit() 
        goids = [dict['GO_id'] for dict in data]
        
        for goid in goids:
            print(goid)
            gosubdag_r0 = GoSubDag([goid], godag, prt=None)
            goids.append(gosubdag_r0.rcntobj.go2ancestors[goid])

        goids = list(set(goids))
        #print(goids)
#        df = pd.DataFrame([r.data() for r in result])
#        if not df.empty:
        if len(goids) > 0 :
#            all_results.append(df)
#            for go in df['GO_id'].values:
            for go in goids:
                if go not in gene_map.keys():
                    gene_map[go] = [gene]
                else:
                    gene_map[go].append(gene)

gene_map_df = pd.DataFrame({'GOid' : gene_map.keys(), 'Genes' : gene_map.values()})

gene_map_df.to_csv('gene_map.tsv', sep = '\t',index = False)
# Concatenate all results into one DataFrame
#gtogo = pd.concat(all_results, ignore_index=True)

#print(gtogo.head())