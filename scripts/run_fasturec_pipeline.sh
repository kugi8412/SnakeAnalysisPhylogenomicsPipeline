#!/bin/bash


INPUT_TREES="$1"
OUTPUT_TREE="$2"
PROCESS_SCRIPT="$3"

if [ ! -s "$INPUT_TREES" ]; then
    echo "WARNING: Input file is empty."
    touch "$OUTPUT_TREE"
    exit 0
fi


CLEAN_INPUT="clean_input_$$.newick"
cp "$INPUT_TREES" "$CLEAN_INPUT"

sed -i -E 's/:[0-9.eE-]+//g' "$CLEAN_INPUT"
sed -i -E 's/\)[0-9.]+/\)/g' "$CLEAN_INPUT"
sed -i -E 's/_p[0-9]+//g' "$CLEAN_INPUT"

# FASTUREC
rm -f fu.tsv
TMP_PREFIX="fasturec_tmp_$$"

echo "DEBUG: Executing fasturec on NORMALIZED input..."
fasturec -G "$CLEAN_INPUT" -Y -b "$TMP_PREFIX" -ot > "${TMP_PREFIX}.log" 2>&1

GENERATED_FILE=$(ls ${TMP_PREFIX}* 2>/dev/null | grep -v ".log" | head -n 1)

# Fallback
if [ -z "$GENERATED_FILE" ]; then
    if [ -f "fu.tsv" ]; then
        GENERATED_FILE="fu.tsv"
    fi
fi

if [ -n "$GENERATED_FILE" ] && [ -s "$GENERATED_FILE" ]; then
    echo "SUCCESS: Found output file: $GENERATED_FILE"
    python3 "$PROCESS_SCRIPT" --fasturec_tree "$GENERATED_FILE" --output_tree "$OUTPUT_TREE" 2>> "${TMP_PREFIX}.log"
    
    if [ ! -s "$OUTPUT_TREE" ]; then
        echo "[WARNING]: Python post-processing failed."
        tail -n 1 "$GENERATED_FILE" > "$OUTPUT_TREE"
    fi
else
    echo "WARNING: Fasturec failed to produce output."
    echo "--- Last 10 lines of log ---"
    tail -n 10 "${TMP_PREFIX}.log"
    touch "$OUTPUT_TREE"
fi

rm -f "${TMP_PREFIX}"* "fu.tsv" "$CLEAN_INPUT"

exit 0
