#!/bin/bash


set -euo pipefail

# Last argument is the output file, all preceding are input files
OUTPUT="${@: -1}"
INPUTS=("${@:1:$#-1}")

cat "${INPUTS[@]}" > "$OUTPUT"
