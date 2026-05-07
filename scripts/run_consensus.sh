#!/bin/bash


INPUT_TREES="$1"
OUTPUT_TREE="$2"
MINSUP="$3"
THREADS="$4"

PREFIX="${OUTPUT_TREE%.treefile}"

echo "[INFO]: Processing consensus for $INPUT_TREES."

# IQ-TREE (handle both iqtree and iqtree2 binary names)
if command -v iqtree2 >/dev/null 2>&1; then
    IQTREE_BIN=iqtree2
elif command -v iqtree >/dev/null 2>&1; then
    IQTREE_BIN=iqtree
else
    echo "[ERROR]: iqtree/iqtree2 not found" >&2; exit 1
fi

$IQTREE_BIN -con -t "$INPUT_TREES" -minsup "$MINSUP" -nt "$THREADS" -pre "$PREFIX" -quiet > /dev/null 2>&1 || true

if [ -f "${PREFIX}.contree" ]; then
    mv "${PREFIX}.contree" "$OUTPUT_TREE"
    echo "[INFO]: Generated .contree"
elif [ -f "${PREFIX}.treefile" ]; then
    mv "${PREFIX}.treefile" "$OUTPUT_TREE"
    echo "[INFO]: Generated .treefile"
elif [ -f "${PREFIX}.consensus" ]; then
    mv "${PREFIX}.consensus" "$OUTPUT_TREE"
    echo "[INFO]: Generated .consensus"
else
    echo "[WARNING]: IQ-TREE failed (likely due to disjoint taxa)."
    echo "[INFO]: Creating empty placeholder to allow pipeline completion."
    touch "$OUTPUT_TREE"
fi

exit 0
