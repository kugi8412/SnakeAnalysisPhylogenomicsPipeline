#!/bin/bash


set -euo pipefail

INPUT="$1"
OUTPUT_TREE="$2"
METHOD="$3"
THREADS="$4"
MODEL="$5"      # Config
BOOTSTRAP="$6"  # Config

PREFIX=${OUTPUT_TREE%.treefile}

if [ "$METHOD" == "ML" ]; then
    # IQ-TREE
    iqtree -s "$INPUT" -m "$MODEL" -bb "$BOOTSTRAP" -nt "$THREADS" -pre "$PREFIX"
    
elif [ "$METHOD" == "NJ" ]; then
    # FastTree
    FastTree -lg < "$INPUT" > "$OUTPUT_TREE"
    
elif [ "$METHOD" == "MP" ]; then
    # RAxML-NG Parsimony
    raxml-ng --start --model "$MODEL" --tree pars{1} --msa "$INPUT" --prefix "$PREFIX" --threads "$THREADS"
    mv "${PREFIX}.raxml.startTree" "$OUTPUT_TREE"
fi
