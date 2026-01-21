#!/usr/bin/env python3
# merge_trees.py
# -*- coding: utf-8 -*-


import sys
import os

import numpy as np
import pandas as pd

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade


try:
    output_file = sys.argv[1]
    PURPOSE = sys.argv[2]               # "for_consensus" OR "for_supertree"
    FILTER_TYPE = sys.argv[3]           # "raw" OR "filtered"
    BS_THRESHOLD = float(sys.argv[4])
    input_trees_list = sys.argv[5:]
except IndexError:
    print("Błąd argumentów.", file=sys.stderr)
    sys.exit(1)

# Config from env variables
try:
    MIN_SPECIES = int(os.environ.get("MIN_SPECIES", 4))
    SAMPLES_CSV = os.environ.get("SAMPLES_CSV", "config/samples.csv")
    OUTGROUPS_STR = os.environ.get("OUTGROUPS", "")
    OUTGROUPS = set(x.strip() for x in OUTGROUPS_STR.split(",") if x.strip())
except:
    MIN_SPECIES = 4
    OUTGROUPS = set()


def get_expected_species(csv_path: str) -> set:
    try:
        df = pd.read_csv(csv_path)
        if 'id' in df.columns: return set(df['id'].astype(str))
        return set(df.iloc[:, 0].astype(str))
    except:
        return set()

ALL_SPECIES = get_expected_species(SAMPLES_CSV)


def force_bifurcation(clade: Clade):
    """(A,B,C) -> ((A,B),C) for FASTARUC
    """
    if clade.is_terminal(): return
    while len(clade.clades) > 2:
        c1 = clade.clades.pop()
        c2 = clade.clades.pop()
        new_int = Clade()
        new_int.clades = [c1, c2]
        clade.clades.append(new_int)

    for child in clade.clades: force_bifurcation(child)

def check_bootstrap(tree: Phylo.BaseTree.Tree,
                    threshold: float
                   ) -> bool:
    """Checking if bootstrap > threshold if filter.
    """
    if FILTER_TYPE == "raw":
        return True

    if threshold <= 0.01: 
        return True
    
    supports = []
    for clade in tree.find_clades():
        val = None
        if clade.confidence is not None: val = clade.confidence
        elif clade.comment:
            try: val = float(clade.comment)
            except: pass
        
        if val is not None:
            if val <= 1.0: val *= 100.0
            supports.append(val)
            
    if not supports:
        return True

    return np.mean(supports) >= threshold

def clean_tree(tree: Phylo.BaseTree.Tree) -> Phylo.BaseTree.Tree:
    """Cleaning Tree (for Fasturec/IQ-TREE).
    """
    try:
        force_bifurcation(tree.root)
    except:
        pass
    
    for clade in tree.find_clades():
        clade.branch_length = None
        clade.confidence = None
        clade.comment = None
        if not clade.is_terminal(): clade.name = None
        else:
            if clade.name:
                clade.name = str(clade.name).strip().replace(":", "_").replace(",", "_").replace(";", "")

    return tree


print(f"--- MERGE: {PURPOSE.upper()} | {FILTER_TYPE.upper()} ---")

count_saved = 0
count_total = 0

with open(output_file, 'w') as out:
    for tree_file in input_trees_list:
        if not os.path.exists(tree_file) or os.path.getsize(tree_file) == 0: continue
        try:
            trees = list(Phylo.parse(tree_file, "newick"))
            if not trees: continue
            tree = trees[0]
            count_total += 1
            leafs = set(leaf.name for leaf in tree.get_terminals())
            species_ok = False
            
            if PURPOSE == "for_consensus":
                if ALL_SPECIES.issubset(leafs):
                    species_ok = True
            else:
                ingroup = leafs - OUTGROUPS
                if len(ingroup) >= MIN_SPECIES:
                    species_ok = True
            
            if not species_ok: continue

            if not check_bootstrap(tree, BS_THRESHOLD):
                continue

            clean_tree(tree)
            Phylo.write(tree, out, "newick")
            count_saved += 1
                
        except Exception as e:
            print(f"Error {tree_file}: {e}", file=sys.stderr)

print(f"Saved {count_saved} out of {count_total} trees.")
