#!/usr/bin/env python3
# select_orthologs.py
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import shutil
import tempfile

from Bio import SeqIO
from multiprocessing import Pool
from collections import defaultdict
from typing import Tuple, List, Dict

SEQ_DB = {}


if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")


def load_all_sequences(fasta_file: str) -> Dict[str, str]:
    print(f"Loading sequences from {fasta_file}...", file=sys.stderr)
    db = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        db[record.id] = str(record.seq)

    return db

def parse_header(header: str) -> Tuple[str, str]:
    # Format: >G01_WP_12345... (SpeciesID_ProteinID)
    parts = header.split('_', 1)
    if len(parts) > 0:
        return parts[0], header

    return "UNKNOWN", header

def get_best_protein(candidates: List[str],
                     rest_of_cluster_proteins: List[str]
                    ) -> str:
    best_hit = candidates[0]
    query_fd, query_path = tempfile.mkstemp(prefix="query_")
    subj_fd, subj_path = tempfile.mkstemp(prefix="subj_")
    try:
        with os.fdopen(query_fd, 'w') as qf:
            for prot_id in candidates:
                if prot_id in SEQ_DB:
                    qf.write(f">{prot_id}\n{SEQ_DB[prot_id]}\n")
        
        count = 0
        with os.fdopen(subj_fd, 'w') as sf:
            for prot_id in rest_of_cluster_proteins:
                if prot_id in SEQ_DB:
                    sf.write(f">{prot_id}\n{SEQ_DB[prot_id]}\n")
                    count += 1
        
        if count > 0:
            cmd = [
                "blastp", 
                "-query", query_path, 
                "-subject", subj_path, 
                "-outfmt", "6 qseqid bitscore", 
                "-max_target_seqs", "1",
                "-evalue", "1e-5"
            ]
            
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                scores = defaultdict(float)
                for line in result.stdout.splitlines():
                    if not line.strip(): continue
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        scores[parts[0]] += float(parts[1])
                
                if scores:
                    best_hit = max(scores, key=scores.get)
            except:
                pass     
    finally:
        if os.path.exists(query_path): os.remove(query_path)
        if os.path.exists(subj_path): os.remove(subj_path)
    
    return best_hit

def process_cluster_job(args: Tuple) -> Tuple[str, any]:
    rep, members_dict, min_species, outgroups, mode = args
    # Singletons
    total_proteins = sum(len(prots) for prots in members_dict.values())
    if total_proteins == 1:
        return ("SINGLETON", None)

    # Min Species Check (Ingroup only)
    ingroup_species = [sp for sp in members_dict.keys() if sp not in outgroups]
    
    if len(ingroup_species) < min_species:
        return ("LOW_SPECIES", None)

    final_sequences = [] 

    if mode == "paralogs":
        for sp, prots in members_dict.items():
            for i, prot_id in enumerate(prots):
                new_header = f"{sp}_p{i+1}" if len(prots) > 1 else sp
                final_sequences.append((new_header, prot_id))

        return ("OK", (rep, final_sequences))
    
    else:
        single_copy_proteins = []
        for sp, prots in members_dict.items():
            if len(prots) == 1:
                single_copy_proteins.append(prots[0])
        
        if not single_copy_proteins:
            return ("NO_ANCHOR", None)

        for sp, prots in members_dict.items():
            if len(prots) == 1:
                final_sequences.append((sp, prots[0]))
            else:
                best_prot = get_best_protein(prots, single_copy_proteins)
                final_sequences.append((sp, best_prot))
            
        return ("OK", (rep, final_sequences))

def main():
    cluster_tsv = snakemake.input.clusters
    all_fasta = snakemake.input.fasta
    out_dir = snakemake.output.out_dir
    stats_file = snakemake.output.stats
    min_species = int(snakemake.params.min_species)
    outgroups = set(snakemake.params.get("outgroups", []))
    mode = snakemake.params.get("mode", "sco")
    threads = snakemake.threads
    
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    os.makedirs(out_dir, exist_ok=True)

    global SEQ_DB
    SEQ_DB = load_all_sequences(all_fasta)
    print(f"Parsing clusters...", file=sys.stderr)
    clusters = defaultdict(lambda: defaultdict(list))

    ALL_PROJECT_SPECIES = set() 
    with open(cluster_tsv, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2: continue
            rep, member = parts[0], parts[1]
            species, _ = parse_header(member)
            clusters[rep][species].append(member)
            if species not in outgroups:
                ALL_PROJECT_SPECIES.add(species)
            
    TOTAL_INGROUP_COUNT = len(ALL_PROJECT_SPECIES)

    tasks = [(rep, members, min_species, outgroups, mode) for rep, members in clusters.items()]
    total_clusters = len(tasks)
    print(f"Processing {total_clusters} clusters in '{mode}' mode with {threads} threads...", file=sys.stderr)
    print(f"Detected {TOTAL_INGROUP_COUNT} unique ingroup species in input data.", file=sys.stderr)
    
    stats = {
        "Total Input": total_clusters,
        "Singletons": 0,
        "Low Species Coverage": 0,
        "No Anchors (SCO only)": 0,
        "Final Families (Total)": 0,
        "Final Families (Full Coverage)": 0,
        "Total Genes Written": 0
    }
    
    with Pool(processes=threads) as pool:
        for result in pool.imap_unordered(process_cluster_job, tasks, chunksize=50):
            status, data = result
            
            if status == "OK":
                stats["Final Families (Total)"] += 1
                rep, seqs = data
                stats["Total Genes Written"] += len(seqs)
                current_species = set()
                for h, _ in seqs:
                    sp_code = h.split('_')[0]
                    if sp_code not in outgroups:
                        current_species.add(sp_code)
                
                if len(current_species) == TOTAL_INGROUP_COUNT:
                    stats["Final Families (Full Coverage)"] += 1
                
                # Save to file
                safe_name = rep.replace("|", "_").replace(":", "_")
                out_file = os.path.join(out_dir, f"{safe_name}.fasta")
                with open(out_file, "w") as f:
                    for header, prot_id in seqs:
                        if prot_id in SEQ_DB:
                            f.write(f">{header}\n{SEQ_DB[prot_id]}\n")
            
            elif status == "SINGLETON":
                stats["Singletons"] += 1
            elif status == "LOW_SPECIES":
                stats["Low Species Coverage"] += 1
            elif status == "NO_ANCHOR":
                stats["No Anchors (SCO only)"] += 1

    print(f"Writing detailed report to {stats_file}...", file=sys.stderr)
    with open(stats_file, "w") as f:
        f.write(f"SELECTION RAPORT\n")
        f.write(f"=========================\n")
        f.write(f"Mode of work:             {mode}\n")
        f.write(f"Required minimum:       {min_species}\n")
        f.write(f"Total species:          {TOTAL_INGROUP_COUNT}\n")
        f.write(f"-------------------------\n")
        f.write(f"All families (Input):        {stats['Total Input']}\n")
        f.write(f"Rejected (Singletons):         {stats['Singletons']}\n")
        f.write(f"Rejected (Low Species Coverage):     {stats['Low Species Coverage']}\n")
        if mode == "sco":
            f.write(f"Rejected (No Anchor):       {stats['No Anchors (SCO only)']}\n")

        f.write(f"-------------------------\n")
        f.write(f"Final Families (Total):   {stats['Final Families (Total)']}\n")
        f.write(f" - Full Coverage:   {stats['Final Families (Full Coverage)']}\n")
        f.write(f" - Partial Coverage: {stats['Final Families (Total)'] - stats['Final Families (Full Coverage)']}\n")
        f.write(f"Total Genes Written:          {stats['Total Genes Written']}\n")

if __name__ == "__main__":
    main()
