#!/usr/bin/env python3
# -*- coding: utf-8 -*- 


import sys
from Bio import Phylo


if len(sys.argv) < 3:
    print("[USAGE]: python clean_leaf_names.py input.newick output.newick")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

try:
    # Read tree
    trees = list(Phylo.parse(input_file, "newick"))
    if not trees:
        print("[ERROR]: Empty input file.")
        sys.exit(1)
    
    tree = trees[0]

    for term in tree.get_terminals():
        original = term.name
        # e.g. "G01_p1" -> "G01", "G01_WP_123" -> "G01"
        if original:
            cleaned = original.split('_')[0]
            term.name = cleaned

    Phylo.write(tree, output_file, "newick")
    print(f"Saved cleaned tree to: {output_file}")

except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
