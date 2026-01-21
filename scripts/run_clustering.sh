#!/bin/bash


set -euo pipefail

INPUT_FASTA="$1"
OUT_PREFIX="$2"
TMP_DIR="$3"
MIN_ID="$4"     # From Config
COV="$5"        # From Config
THREADS="$6"    # From Config

echo "[MMseqs2] Params: MinID=$MIN_ID, Cov=$COV, Threads=$THREADS"

mmseqs easy-cluster \
    "$INPUT_FASTA" \
    "$OUT_PREFIX" \
    "$TMP_DIR" \
    --min-seq-id "$MIN_ID" \
    -c "$COV" \
    --cov-mode 0 \
    --threads "$THREADS"

if [ -d "$TMP_DIR" ]; then
    rm -rf "$TMP_DIR"
fi
