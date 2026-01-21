#!/usr/bin/env python3
# normalize_genetrees.py
# -*- coding: utf-8 -*-


import sys
from Bio import Phylo


input_file = sys.argv[1]
output_file = sys.argv[2]

print(f"Normalizacja nazw w: {input_file} -> {output_file}")

try:
    trees = list(Phylo.parse(input_file, "newick"))
    
    if not trees:
        print("BŁĄD: Plik wejściowy jest pusty.")
        sys.exit(0)

    for tree in trees:
        for term in tree.get_terminals():
            term.branch_length = None
            term.confidence = None
            term.comment = None

            if term.name:
                term.name = term.name.split('_')[0]
    
    Phylo.write(trees, output_file, "newick")

except Exception as e:
    print(f"Błąd krytyczny: {e}")
    sys.exit(1)
