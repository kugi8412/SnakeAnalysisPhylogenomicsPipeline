#!/usr/bin/env python3
# concat_matrix.py
# -*- coding: utf-8 -*- 

import sys
import os

from Bio import SeqIO


input_files = snakemake.input
output_matrix = snakemake.output.matrix
output_partitions = snakemake.output.partitions

print(f"Concatenating {len(input_files)} alignments...")

taxa_sequences = {}
all_taxa = set()
partitions = []
current_position = 1

for aln_file in input_files:
    gene_name = os.path.basename(aln_file).replace(".trim.aln", "")
    gene_seqs = {}
    aln_length = 0
    for record in SeqIO.parse(aln_file, "fasta"):
        taxon_id = record.id.split('_')[0]
        gene_seqs[taxon_id] = str(record.seq)
        aln_length = len(record.seq)
        all_taxa.add(taxon_id)

    if aln_length == 0:
        continue

    start = current_position
    end = current_position + aln_length - 1
    partitions.append((gene_name, start, end))
    current_position = end + 1

genes_data = []
all_taxa = set()
genes_data = []

for aln_file in input_files:
    gene_name = os.path.basename(aln_file).replace(".trim.aln", "")
    gene_dict = {}
    length = 0
    
    try:
        for record in SeqIO.parse(aln_file, "fasta"):
            taxon = record.id.split('_')[0]
            gene_dict[taxon] = str(record.seq)
            length = len(record.seq)
            all_taxa.add(taxon)
            
        if length > 0:
            genes_data.append((gene_name, length, gene_dict))
            
    except Exception as e:
        print(f"[Warning]: Error reading {aln_file}: {e}", file=sys.stderr)

sorted_taxa = sorted(list(all_taxa))
final_matrix = {t: [] for t in sorted_taxa}
final_partitions = []
curr_pos = 1

for gene_name, length, seqs in genes_data:
    end_pos = curr_pos + length - 1
    # RAxML/IQ-TREE
    final_partitions.append(f"JTT+G4, {gene_name} = {curr_pos}-{end_pos}")
    curr_pos = end_pos + 1
    
    # Sequences
    for taxon in sorted_taxa:
        if taxon in seqs:
            final_matrix[taxon].append(seqs[taxon])
        else:
            final_matrix[taxon].append("-" * length)

total_len = curr_pos - 1

with open(output_matrix, "w") as f:
    f.write(f"{len(sorted_taxa)} {total_len}\n")
    for taxon in sorted_taxa:
        full_seq = "".join(final_matrix[taxon])
        f.write(f"{taxon}  {full_seq}\n")

# Partitions
with open(output_partitions, "w") as f:
    for p in final_partitions:
        f.write(f"{p}\n")

print(f"Supermatrix saved to {output_matrix}")
print(f"Partitions saved to {output_partitions}")
