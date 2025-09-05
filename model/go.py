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
genes = genes[0:10]
gene_map = {}
with driver.session() as session:
    for gene in genes:
        print(gene)
        result = session.run(query, {"gene_name": gene, "taxon": taxon})
        data = [r.data() for r in result]
#        print(data[])
#        exit() 

        ## extract a list of all goids associated with specific gene
        goids = [dict['GO_id'] for dict in data]
        anc=[]
        for goid in goids:
         #   if goid == 'GO:0008150': ## skip for root term BP
         #       continue
                     # Add direct gene to GO term
            if goid not in gene_map:
                gene_map[goid] = []
            if gene not in gene_map[goid]:
                gene_map[goid].append(gene)

            gosubdag_r0 = GoSubDag([goid], godag, prt=None)
            try:
                anc += list(gosubdag_r0.rcntobj.go2ancestors[goid])
            except KeyError:
                continue

        anc = list(set(anc))
        
        ## loop through all the ancestor terms for gene
        ## if the go term is not in the gene map dictionary keys add it as a key and assign gene
        if len(anc) > 0 :
#            all_results.append(df)
#            for go in df['GO_id'].values:
            for go in anc:
                if go not in gene_map.keys():
                    gene_map[go] = [gene]
                else:
                    gene_map[go].append(gene)

with open('PGS001990/ontology.txt', 'w') as f:
    for goid in gene_map.keys():
        ngenes = len(gene_map[goid])
        go_obj = godag[goid]

        # Get direct children
        direct_children = set()
        for child in go_obj.get_goterms_lower():
        # Check if this is a direct child
            if goid in [parent.id for parent in child.parents]:
                direct_children.add(child.id) 

        f.write(f'ROOT: {goid} {ngenes}\n')
        f.write(f'GENES: {' '.join(gene_map[goid])}\n')
        f.write(f'TERMS: {' '.join(direct_children)}\n')

    # Store the direct children
#    for child in go_obj.children:
#        child_id = child.id

#        print(child_id)
#gene_map_df = pd.DataFrame({'GOid' : gene_map.keys(), 'Genes' : gene_map.values()})

#gene_map_df.to_csv('gene_map.tsv', sep = '\t',index = False)

#print(gene_map_df)

# Concatenate all results into one DataFrame
#gtogo = pd.concat(all_results, ignore_index=True)

#print(gtogo.head())