"""
OBO File Parser for GO Terms

This module provides utilities to parse Gene Ontology OBO files and extract
term names and metadata.

The OBO format is used by the Gene Ontology Consortium to distribute
ontology data in a human-readable text format.

Functions:
    parse_obo_file: Parse OBO file and extract term information
    get_term_name: Get human-readable name for a GO term ID
    create_term_name_dict: Create mapping from GO IDs to names
"""

import re
from collections import defaultdict


def parse_obo_file(obo_file):
    """
    Parse a Gene Ontology OBO file.

    The OBO format consists of stanzas (blocks) for each term, with fields like:
    - id: The GO identifier (e.g., GO:0000001)
    - name: Human-readable term name
    - namespace: Category (biological_process, molecular_function, cellular_component)
    - def: Definition of the term
    - is_a: Parent terms
    - is_obsolete: Whether the term is deprecated

    Args:
        obo_file (str): Path to OBO file

    Returns:
        dict: Maps GO IDs to term information dictionaries
              Each dict contains: name, namespace, definition, is_obsolete, parents

    Example:
        >>> terms = parse_obo_file('go-basic.obo')
        >>> print(terms['GO:0000001'])
        {'name': 'mitochondrion inheritance',
         'namespace': 'biological_process',
         'definition': 'The distribution of mitochondria...',
         'is_obsolete': False,
         'parents': ['GO:0048308', 'GO:0048311']}
    """
    terms = {}
    current_term = None
    in_term_stanza = False

    with open(obo_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()

            # Start of a new term stanza
            if line == '[Term]':
                in_term_stanza = True
                current_term = {
                    'id': None,
                    'name': None,
                    'namespace': None,
                    'definition': None,
                    'is_obsolete': False,
                    'parents': []
                }
                continue

            # End of stanza
            if line == '' and in_term_stanza:
                if current_term and current_term['id']:
                    # Skip obsolete terms
                    if not current_term['is_obsolete']:
                        terms[current_term['id']] = current_term
                in_term_stanza = False
                current_term = None
                continue

            # Skip non-term stanzas
            if not in_term_stanza:
                continue

            # Parse fields within term stanza
            if line.startswith('id:'):
                current_term['id'] = line.split('id:')[1].strip()

            elif line.startswith('name:'):
                current_term['name'] = line.split('name:')[1].strip()

            elif line.startswith('namespace:'):
                current_term['namespace'] = line.split('namespace:')[1].strip()

            elif line.startswith('def:'):
                # Definition is in quotes
                match = re.search(r'"([^"]*)"', line)
                if match:
                    current_term['definition'] = match.group(1)

            elif line.startswith('is_obsolete:'):
                if 'true' in line.lower():
                    current_term['is_obsolete'] = True

            elif line.startswith('is_a:'):
                # Extract parent GO ID
                parent_id = line.split('is_a:')[1].strip().split()[0]
                current_term['parents'].append(parent_id)

    return terms


def create_term_name_dict(obo_file):
    """
    Create a simple mapping from GO IDs to names.

    Args:
        obo_file (str): Path to OBO file

    Returns:
        dict: Maps GO ID (str) to term name (str)

    Example:
        >>> name_dict = create_term_name_dict('go-basic.obo')
        >>> print(name_dict['GO:0000001'])
        'mitochondrion inheritance'
    """
    terms = parse_obo_file(obo_file)
    return {term_id: info['name'] for term_id, info in terms.items() if info['name']}


def get_term_name(term_id, name_dict, max_length=None):
    """
    Get human-readable name for a GO term ID.

    If the term is not in the dictionary, returns the original ID.
    Optionally truncates long names for visualization.

    Args:
        term_id (str): GO term ID (e.g., 'GO:0000001')
        name_dict (dict): Mapping from GO IDs to names
        max_length (int, optional): Maximum name length. If longer, truncate with '...'

    Returns:
        str: Term name (or ID if not found)

    Example:
        >>> name_dict = {'GO:0000001': 'mitochondrion inheritance'}
        >>> get_term_name('GO:0000001', name_dict)
        'mitochondrion inheritance'
        >>> get_term_name('GO:0000001', name_dict, max_length=15)
        'mitochondrion...'
        >>> get_term_name('GO:9999999', name_dict)
        'GO:9999999'
    """
    # Get name from dictionary, or use ID if not found
    name = name_dict.get(term_id, term_id)

    # Truncate if needed
    if max_length and len(name) > max_length:
        name = name[:max_length-3] + '...'

    return name


def format_term_label(term_id, name_dict, include_id=True, max_length=50):
    """
    Format a term label for visualization.

    Creates labels like "GO:0000001: mitochondrion inheritance"
    or just "mitochondrion inheritance" depending on include_id.

    Args:
        term_id (str): GO term ID
        name_dict (dict): Mapping from GO IDs to names
        include_id (bool): Whether to include the GO ID in label
        max_length (int): Maximum name length before truncation

    Returns:
        str: Formatted label

    Example:
        >>> name_dict = {'GO:0000001': 'mitochondrion inheritance'}
        >>> format_term_label('GO:0000001', name_dict, include_id=True)
        'GO:0000001: mitochondrion inheritance'
        >>> format_term_label('GO:0000001', name_dict, include_id=False)
        'mitochondrion inheritance'
    """
    name = get_term_name(term_id, name_dict, max_length=max_length)

    if include_id and term_id != name:  # Only include ID if we found a name
        return f"{term_id}: {name}"
    else:
        return name


if __name__ == "__main__":
    # Test the parser
    import sys

    if len(sys.argv) > 1:
        obo_file = sys.argv[1]
        print(f"Parsing {obo_file}...")

        terms = parse_obo_file(obo_file)
        print(f"Found {len(terms)} non-obsolete terms")

        # Show first few terms
        print("\nFirst 5 terms:")
        for i, (term_id, info) in enumerate(list(terms.items())[:5]):
            print(f"\n{term_id}:")
            print(f"  Name: {info['name']}")
            print(f"  Namespace: {info['namespace']}")
            print(f"  Parents: {info['parents']}")

        # Create name dictionary
        name_dict = create_term_name_dict(obo_file)
        print(f"\nCreated name dictionary with {len(name_dict)} entries")
    else:
        print("Usage: python obo_parser.py <obo_file>")
        print("Example: python obo_parser.py go-basic.obo")
