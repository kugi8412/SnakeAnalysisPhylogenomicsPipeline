#!/bin/bash


INPUT="$1"
OUTPUT="$2"
METHOD="$3"

fallback() {
    echo "WARNING: Trimal failed for $INPUT. Creating empty output."
    touch "$OUTPUT"
    exit 0
}

if [ ! -s "$INPUT" ]; then
    echo "WARNING: Empty input for trimal ($INPUT). Skipping."
    touch "$OUTPUT"
    exit 0
fi

trimal -in "$INPUT" -out "$OUTPUT" $METHOD 2>/dev/null || fallback

exit 0
