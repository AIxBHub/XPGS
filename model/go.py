from neo4j import GraphDatabase
import pandas as pd
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag

# Load GO DAG
godag = GODag('data/GO/go-basic.obo', optional_attrs={'relationship'})

annotation_file = "PGS001990/ukb_imp_bct_vars_merged_clean_annotations-NodeNorm.csv"

effects = ['missense_variant','splice_region_variant','stop_gained','5_prime_UTR_premature_start_codon_gain_variant']
taxon = "NCBITaxon:9606"  # homo sapien tax
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

# Read annotations and get genes
afile = pd.read_csv(annotation_file)
genes = afile[afile['ann'].str.contains('|'.join(effects))]['name'].sort_values().unique()

# For testing, limit to 10 genes
genes = genes[0:10]

# Create a gene-to-id mapping file
gene_to_id = {}
for i, gene in enumerate(genes):
    gene_to_id[gene] = i

# Save gene2id.txt
with open('gene2id.txt', 'w') as f:
    for gene, gene_id in gene_to_id.items():
        f.write(f"{gene} {gene_id}\n")

# Build GO term to gene mappings
go_term_relationships = []  # Store parent-child term relationships
go_gene_associations = []   # Store term-gene associations
go_term_hierarchy = {}      # Store the hierarchy for validation

with driver.session() as session:
    for gene in genes:
        print(f"Processing gene: {gene}")
        result = session.run(query, {"gene_name": gene, "taxon": taxon})
        data = [r.data() for r in result]
        
        # Extract GO terms directly associated with this gene
        goids = [d['GO_id'] for d in data]
        
        # Add direct gene to GO term associations
        for goid in goids:
            go_gene_associations.append((goid, gene, "gene"))
            
            # Get ancestors for this GO term
            try:
                gosubdag_r0 = GoSubDag([goid], godag, prt=None)
                ancestors = list(gosubdag_r0.rcntobj.go2ancestors[goid])
                
                # For each ancestor, add the gene association
                for ancestor in ancestors:
                    go_gene_associations.append((ancestor, gene, "gene"))
                    
                    # Get parent-child relationships for the GO hierarchy
                    go_obj = godag[goid]
                    for parent in go_obj.parents:
                        parent_id = parent.id
                        if (parent_id, goid) not in go_term_relationships:
                            go_term_relationships.append((parent_id, goid, "default"))
                            
                            # Track hierarchy for validation
                            if parent_id not in go_term_hierarchy:
                                go_term_hierarchy[parent_id] = []
                            if goid not in go_term_hierarchy[parent_id]:
                                go_term_hierarchy[parent_id].append(goid)
                
            except KeyError:
                print(f"Warning: Could not process GO term {goid}")
                continue

# Find the root term(s)
all_terms = set([rel[0] for rel in go_term_relationships] + [rel[1] for rel in go_term_relationships])
child_terms = set([rel[1] for rel in go_term_relationships])
root_terms = all_terms - child_terms

# If multiple roots, create a single ROOT node
if len(root_terms) > 1:
    for root in root_terms:
        go_term_relationships.append(("ROOT", root, "default"))
elif len(root_terms) == 1:
    root = list(root_terms)[0]
    # Rename the single root to "ROOT" if needed
    go_term_relationships = [("ROOT" if rel[0] == root else rel[0], 
                             "ROOT" if rel[1] == root else rel[1], 
                             rel[2]) for rel in go_term_relationships]
else:
    print("Warning: No root terms found")

# Write the ontology file in the required format
with open('ontology.txt', 'w') as f:
    # Write term-term relationships first
    for parent, child, rel_type in go_term_relationships:
        f.write(f"{parent}\t{child}\t{rel_type}\n")
    
    # Write term-gene associations
    for term, gene, rel_type in go_gene_associations:
        # Only include genes that are in our gene_to_id mapping
        if gene in gene_to_id:
            f.write(f"{term}\t{gene}\t{rel_type}\n")

print("Ontology file generated successfully")