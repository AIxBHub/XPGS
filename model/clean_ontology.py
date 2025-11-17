#!/usr/bin/env python3
"""
Clean Ontology Utility Script

This script removes empty terms (terms with no gene annotations) from an ontology file.
It's useful for preprocessing ontology files before VNN training to avoid warnings and
potential training issues.

Usage:
    python clean_ontology.py -onto ontology.txt -gene2id gene2ind.txt -output ontology_cleaned.txt

Input:
    - Ontology file with format: parent\tchild\trelationship_type
    - Gene-to-ID mapping file with format: gene_name\tgene_id

Output:
    - Cleaned ontology file with empty terms removed
    - Statistics report about the cleaning process
"""

import argparse
import sys
import logging
from data_wrapper import TrainingDataWrapper


def setup_logger():
    """Configure logging for the cleaning process."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s: %(message)s'
    )
    return logging.getLogger(__name__)


def main():
    """Main function for ontology cleaning."""
    parser = argparse.ArgumentParser(
        description='Remove empty terms from ontology file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Clean ontology file
  python clean_ontology.py -onto data/ontology.txt -gene2id data/gene2ind.txt -output data/ontology_cleaned.txt

  # Print statistics only (no output file)
  python clean_ontology.py -onto data/ontology.txt -gene2id data/gene2ind.txt --dry-run

Note:
  Empty terms are those with no gene annotations, either directly or
  through descendant terms in the ontology hierarchy.
        """
    )

    # Required arguments
    parser.add_argument('-onto', '--ontology', required=True,
                        help='Input ontology file path')
    parser.add_argument('-gene2id', '--gene-mapping', required=True,
                        help='Gene-to-ID mapping file path')
    parser.add_argument('-output', '--output-file', required=True,
                        help='Output cleaned ontology file path')

    # Optional arguments
    parser.add_argument('--dry-run', action='store_true',
                        help='Only report statistics without creating output file')
    parser.add_argument('--save-empty-list', type=str,
                        help='Save list of empty terms to this file')

    args = parser.parse_args()

    # Setup logger
    logger = setup_logger()

    logger.info("Ontology Cleaning Utility")
    logger.info("="*60)

    # Create minimal args object for TrainingDataWrapper
    # We need to provide dummy values for training-related parameters
    class MinimalArgs:
        def __init__(self, ontology_file, gene2id_file):
            self.onto = ontology_file
            self.gene2id = gene2id_file
            # Dummy training parameters (not used for cleaning)
            self.genotype_hiddens = 50
            self.min_dropout_layer = 2
            self.dropout_fraction = 0.3
            self.lr = 0.001
            self.wd = 0.001
            self.alpha = 0.1
            self.epoch = 1
            self.batchsize = 64
            self.modeldir = 'MODEL'
            self.delta = 0.01
            self.train = None
            self.test = None
            self.testsetratio = 0.2
            self.optimize = None
            self.metric_output = 'metrics.tsv'

    try:
        # Initialize data wrapper (this loads the ontology)
        logger.info(f"Loading ontology from: {args.ontology}")
        logger.info(f"Loading gene mapping from: {args.gene_mapping}")

        minimal_args = MinimalArgs(args.ontology, args.gene_mapping)
        wrapper = TrainingDataWrapper(minimal_args, logger=logger)

        # Check if there are empty terms
        if not hasattr(wrapper, 'dG') or len(wrapper.dG.nodes()) == 0:
            logger.error("Failed to load ontology. Please check file format.")
            sys.exit(1)

        # Run the cleaning process
        if args.dry_run:
            logger.info("\nDRY RUN MODE - No output file will be created")

        stats = wrapper.remove_empty_terms(args.ontology, args.output_file if not args.dry_run else '/tmp/dummy.txt')

        # Display results
        print("\n" + "="*60)
        print("CLEANING RESULTS")
        print("="*60)
        print(f"Original terms:               {stats['original_terms']}")
        print(f"Empty terms removed:          {stats['empty_terms_removed']}")
        print(f"Final terms:                  {stats['final_terms']}")
        print(f"Reduction:                    {stats['empty_terms_removed']/stats['original_terms']*100:.1f}%")
        print(f"Term edges preserved:         {stats['term_edges_preserved']}")
        print(f"Gene annotations preserved:   {stats['gene_annotations_preserved']}")

        if not args.dry_run:
            print(f"\nCleaned ontology saved to:    {args.output_file}")
        print("="*60)

        # Save empty term list if requested
        if args.save_empty_list and len(stats['empty_term_list']) > 0:
            with open(args.save_empty_list, 'w') as f:
                f.write("# Empty terms removed from ontology\n")
                f.write(f"# Total: {len(stats['empty_term_list'])}\n")
                for term in stats['empty_term_list']:
                    f.write(f"{term}\n")
            logger.info(f"Empty term list saved to: {args.save_empty_list}")

        # Print sample of empty terms if any were found
        if stats['empty_terms_removed'] > 0:
            print(f"\nSample of removed empty terms (showing up to 10):")
            for term in stats['empty_term_list'][:10]:
                print(f"  - {term}")
            if len(stats['empty_term_list']) > 10:
                print(f"  ... and {len(stats['empty_term_list']) - 10} more")

        logger.info("\nCleaning completed successfully!")

    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error during cleaning: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
