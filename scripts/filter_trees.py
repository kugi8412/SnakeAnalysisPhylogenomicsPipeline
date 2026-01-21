#!/usr/bin/env python3
# merge_trees.py
# -*- coding: utf-8 -*-


import sys
import numpy as np
import pandas as pd

from Bio import Phylo

# Snakemake args
input_tree_path = sys.argv[1]
threshold = float(sys.argv[2])
samples_file = sys.argv[3]

try:
    df = pd.read_csv(samples_file)
    if 'id' in df.columns:
        expected_species = set(df['id'].astype(str))
    else:
        expected_species = set(df.iloc[:, 0].astype(str))
    
    expected_count = len(expected_species)
    trees = list(Phylo.parse(input_tree_path, "newick"))

    if not trees:
        sys.exit(0)

    tree = trees[0]

    leafs = set(leaf.name for leaf in tree.get_terminals())
    if not expected_species.issubset(leafs):
        sys.exit(0)

    # Bootstrap check
    supports = []
    for clade in tree.find_clades():
        if clade.confidence is not None:
            val = clade.confidence
            if val <= 1.0 and threshold > 1.0:
                val *= 100.0
            supports.append(val)
        elif clade.comment:
            try:
                val = float(clade.comment)
                if val <= 1.0 and threshold > 1.0:
                    val *= 100.0
                supports.append(val)
            except: pass
    
    avg_support = np.mean(supports) if supports else 0.0

    if avg_support >= threshold:
        print(input_tree_path)
    else:
        sys.exit(0)

except Exception as e:
    sys.exit(0)
