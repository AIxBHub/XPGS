from neo4j import GraphDatabase
import pandas as pd
from goatools.obo_parser import GODag
from goatools.gosubdag.gosubdag import GoSubDag
import argparse

def main():
    parser = argparse.ArgumentParser(description = 'Build ontology file')
    parser.add_argument('-obo', help = 'GO DAG file', type = str)
    parser.add_argument('-anno', help = 'Variant annotation file', type = str)
    parser.add_argument('-eff', nargs = '+', help = 'List of effects', type = str, default = ['missense_variant','splice_region_variant','stop_gained','5_prime_UTR_premature_start_codon_gain_variant'])
    parser.add_argument('-tax', help = 'Taxon', type = str, default = 'NCBITaxon:9606')
    parser.add_argument('-kg', help = 'URL to KG', type = str, default = "bolt://robokopkg.renci.org:7687")
    parser.add_argument('-test', help = 'Test with N genes', type = int)
    parser.add_argument('-out', help = 'path to output directory', type = str, required=True)

    args = parser.parse_args()

    # Load GO DAG
    godag = GODag(args.obo, optional_attrs={'relationship'})
    driver = GraphDatabase.driver(args.kg, auth=("neo4j", "password"))

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
    afile = pd.read_csv(args.anno)
    genes = afile[afile['ann'].str.contains('|'.join(args.eff))]['name'].sort_values().unique()

    if args.test is not None:
        genes = genes[0:args.test]

    gene_to_id = map_genes(genes, args.out)
    go_term_relationships, go_term_hierarchy, go_gene_associations = map_terms(driver, genes, query, args.tax, godag)
    go_term_relationships = find_root(go_term_relationships)
    save_ontologies(args.out, go_term_relationships, go_gene_associations, gene_to_id)

def map_genes(genes, outdir):
    # Create a gene-to-id mapping file
    gene_to_id = {}
    for i, gene in enumerate(genes):
        if pd.isna(gene):
            continue
        if ' ' in gene:
            gene = gene.split(' ')[0]
        gene_to_id[gene] = i

    # Save gene2id.txt
    with open(f'{outdir}/gene2id.txt', 'w') as f:
        for gene, gene_id in gene_to_id.items():
            f.write(f"{gene} {gene_id}\n")

    return gene_to_id

def map_terms(driver, genes, query, taxon, godag):
    # Build GO term to gene mappings using sets to avoid duplicates
    go_term_relationships = set()  # Store parent-child term relationships
    go_gene_associations = set()   # Store term-gene associations
    go_term_hierarchy = {}         # Store the hierarchy for validation
    processed_terms = set()        # Track processed GO terms to avoid duplicate hierarchy processing

    with driver.session() as session:
        for gene in genes:
            result = session.run(query, {"gene_name": gene, "taxon": taxon})
            data = [r.data() for r in result]

            # Extract GO terms directly associated with this gene
            goids = [d['GO_id'] for d in data]

            # Add direct gene to GO term associations
            for goid in goids:
                go_gene_associations.add((goid, gene, "gene"))

                # Get ancestors for this GO term
                try:
                    gosubdag_r0 = GoSubDag([goid], godag, prt=None)
                    ancestors = list(gosubdag_r0.rcntobj.go2ancestors[goid])

                    # For each ancestor, add the gene association
                    for ancestor in ancestors:
                        go_gene_associations.add((ancestor, gene, "gene"))

                    # Process parent-child relationships only once per GO term
                    if goid not in processed_terms:
                        # Get parent-child relationships for the GO hierarchy
                        go_obj = godag[goid]
                        for parent in go_obj.parents:
                            parent_id = parent.id
                            go_term_relationships.add((parent_id, goid, "default"))

                            # Track hierarchy for validation
                            if parent_id not in go_term_hierarchy:
                                go_term_hierarchy[parent_id] = []
                            if goid not in go_term_hierarchy[parent_id]:
                                go_term_hierarchy[parent_id].append(goid)
                        
                        processed_terms.add(goid)

                except KeyError:
                    print(f"Warning: Could not process GO term {goid}")
                    continue

    # Convert sets back to lists for compatibility with the rest of the code
    return list(go_term_relationships), go_term_hierarchy, list(go_gene_associations)

def find_root(go_term_relationships):
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
    
    return go_term_relationships

def save_ontologies(outdir, go_term_relationships, go_gene_associations, gene_to_id):
    # Write the ontology file in the required format
    with open(f'{outdir}/ontology.txt', 'w') as f:
        # Write term-term relationships first
        for parent, child, rel_type in go_term_relationships:
            f.write(f"{parent}\t{child}\t{rel_type}\n")
        
        # Write term-gene associations
        for term, gene, rel_type in go_gene_associations:
            # Only include genes that are in our gene_to_id mapping
            if gene in gene_to_id:
                f.write(f"{term}\t{gene}\t{rel_type}\n")

if __name__ == "__main__":
    main()  