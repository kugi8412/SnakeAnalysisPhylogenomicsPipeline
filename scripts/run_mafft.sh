#!/bin/bash


INPUT="$1"
OUTPUT="$2"
THREADS="$3"
ARGS="$4"

fallback() {
    echo "[WARNING]: Alignment failed for $INPUT. Creating fallback output."
    if [ -s "$INPUT" ]; then
        cat "$INPUT" > "$OUTPUT"
    else
        touch "$OUTPUT"
    fi
    exit 0
}

if [ ! -s "$INPUT" ]; then
    echo "WARNING: Empty input $INPUT"
    touch "$OUTPUT"
    exit 0
fi

NUM_SEQS=$(grep -c "^>" "$INPUT" || echo "0")

if [ "$NUM_SEQS" -lt 2 ]; then
    echo "WARNING: Only $NUM_SEQS seq in $INPUT. Skipping alignment."
    cat "$INPUT" > "$OUTPUT"
    exit 0
fi

# MAFFT
mafft --thread "$THREADS" $ARGS "$INPUT" > "$OUTPUT" 2>/dev/null || fallback

exit 0
