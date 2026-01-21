#!/usr/bin/env python3
# run_genetree.py
# -*- coding: utf-8 -*-


import os
import sys
import shlex
import shutil
import subprocess

# Snakemake Params
input_aln = snakemake.input.aln
output_tree = snakemake.output.tree
method = snakemake.params.method
raw_flags = snakemake.params.flags
extra_flags = shlex.split(raw_flags) if raw_flags else []
threads = snakemake.threads

tmp_clean_aln = input_aln + f".clean.{os.getpid()}.tmp"


def clean_and_validate_fasta(in_path: str,
                             out_path: str
                            ) -> int:
    valid_count = 0
    with open(in_path, 'r') as fin, open(out_path, 'w') as fout:
        header = None
        seq_lines = []
        
        def process_entry(h, lines):
            if not h: return 0
            full_seq = "".join(lines).strip()
            meaningful_seq = full_seq.replace("-", "").replace("?", "").replace("X", "").replace("x", "")
            
            if len(meaningful_seq) > 0:
                fout.write(f"{h}\n{full_seq}\n")
                return 1
            return 0

        for line in fin:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if header:
                    valid_count += process_entry(header, seq_lines)

                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        
        # Last sequence
        if header:
            valid_count += process_entry(header, seq_lines)
            
    return valid_count

if not os.path.exists(input_aln) or os.path.getsize(input_aln) == 0:
    print(f"[WARNING]: Skipping empty alignment {input_aln}", file=sys.stderr)
    open(output_tree, 'w').close()
    sys.exit(0)

try:
    num_seqs = clean_and_validate_fasta(input_aln, tmp_clean_aln)
except Exception as e:
    print(f"[ERROR]: Failed to clean fasta {input_aln}: {e}", file=sys.stderr)
    open(output_tree, 'w').close()
    sys.exit(0)

if num_seqs < 4:
    print(f"[WARNING]: Too few valid sequences ({num_seqs}) in {input_aln} after cleaning. Skipping.", file=sys.stderr)
    if os.path.exists(tmp_clean_aln):
        os.remove(tmp_clean_aln)
    open(output_tree, 'w').close()
    sys.exit(0)

# Buildingtree
cmd = []
if method == "NJ":
    binary = "FastTree" 
    if shutil.which("FastTree") is None and shutil.which("fasttree") is not None:
        binary = "fasttree"

    cmd = [binary] + extra_flags + [tmp_clean_aln]

elif method == "ML":
    prefix = output_tree.replace(".treefile", "")
    cmd = ["iqtree", "-s", tmp_clean_aln, "-T", str(threads), "-pre", prefix, "-quiet"] + extra_flags

print(f"DEBUG: Running {' '.join(cmd)}")

try:
    if method == "NJ":
        with open(output_tree, "w") as f:
            subprocess.run(cmd, stdout=f, check=True)
    else:
        # IQ-TREE
        subprocess.run(cmd, check=True)
        # IQ-TREE => .treefile
        expected_out = prefix + ".treefile"
        if os.path.exists(expected_out) and expected_out != output_tree:
             shutil.move(expected_out, output_tree)

except subprocess.CalledProcessError as e:
    print(f"ERROR: Tree tool failed for {input_aln}. Cmd: {e.cmd}", file=sys.stderr)
    open(output_tree, 'w').close()

finally:
    if os.path.exists(tmp_clean_aln):
        os.remove(tmp_clean_aln)
    
    if method == "ML":
        prefix = output_tree.replace(".treefile", "")
        for ext in [".iqtree", ".log", ".bionj", ".mldist", ".ckp.gz", ".contree"]:
            f = prefix + ext
            if os.path.exists(f):
                os.remove(f)
