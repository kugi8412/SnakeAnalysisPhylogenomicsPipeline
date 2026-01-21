#!/usr/bin/env python3
import sys
from Bio import Phylo

if len(sys.argv) < 3:
    print("Użycie: python clean_leaf_names.py input.newick output.newick")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

try:
    # Wczytujemy drzewo
    trees = list(Phylo.parse(input_file, "newick"))
    if not trees:
        print("Błąd: Pusty plik wejściowy.")
        sys.exit(1)
    
    tree = trees[0]

    # Iterujemy po liściach i czyścimy nazwy
    for term in tree.get_terminals():
        original = term.name
        # Bierzemy wszystko przed pierwszym podkreślnikiem "_"
        # Np. "G01_p1" -> "G01", "G01_WP_123" -> "G01"
        if original:
            cleaned = original.split('_')[0]
            term.name = cleaned
    
    # Zapisujemy
    Phylo.write(tree, output_file, "newick")
    print(f"Zapisano wyczyszczone drzewo do: {output_file}")

except Exception as e:
    print(f"Błąd: {e}")
    sys.exit(1)
