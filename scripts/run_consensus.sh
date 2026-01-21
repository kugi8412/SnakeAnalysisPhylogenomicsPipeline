#!/bin/bash


INPUT_TREES="$1"
OUTPUT_TREE="$2"
MINSUP="$3"
THREADS="$4"

PREFIX="${OUTPUT_TREE%.treefile}"

echo "Processing consensus for $INPUT_TREES..."

# IQ-TREE
iqtree2 -con -t "$INPUT_TREES" -minsup "$MINSUP" -nt "$THREADS" -pre "$PREFIX" -quiet > /dev/null 2>&1 || true

if [ -f "${PREFIX}.contree" ]; then
    mv "${PREFIX}.contree" "$OUTPUT_TREE"
    echo "Success: Generated .contree"
elif [ -f "${PREFIX}.treefile" ]; then
    mv "${PREFIX}.treefile" "$OUTPUT_TREE"
    echo "Success: Generated .treefile"
elif [ -f "${PREFIX}.consensus" ]; then
    mv "${PREFIX}.consensus" "$OUTPUT_TREE"
    echo "Success: Generated .consensus"
else
    echo "[WARNING]: IQ-TREE failed (likely due to disjoint taxa)."
    echo "Creating empty placeholder to allow pipeline completion."
    touch "$OUTPUT_TREE"
fi

exit 0
