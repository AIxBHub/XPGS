from neo4j import GraphDatabase
import pandas as pd
print('starting...')
driver = GraphDatabase.driver("bolt://robokopkg.renci.org:7687", auth=("neo4j", "password"))
print(driver)
query = """
MATCH (n:`biolink:GeneOrGeneProduct`)-[r]-(d)
WHERE n.name = $gene_name
  AND n.taxon = $taxon
  AND r.primary_knowledge_source = "infores:goa"
  AND NOT "biolink:CellularComponent" IN labels(d)
RETURN 
    n.id AS gene_id,
    n.name AS gene_name,
    type(r) AS predicate,
    r.primary_knowledge_source AS pks,
    d.id AS GO_id,
    d.name AS GO_name
"""

genes = ["SSU72", "BRCA1", "TP53"]
taxon = "NCBITaxon:9606"

all_results = []

with driver.session() as session:
    for gene in genes:
        result = session.run(query, {"gene_name": gene, "taxon": taxon})
        df = pd.DataFrame([r.data() for r in result])
        print(df)
        if not df.empty:
            all_results.append(df)

# Concatenate all results into one DataFrame
gtogo = pd.concat(all_results, ignore_index=True)

print(gtogo.head())