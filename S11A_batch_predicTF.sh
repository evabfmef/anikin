#!/bin/bash

# PredicTF Batch Processing Script
# Runs PredicTF analysis on all .fa/.faa protein sequence files in a directory.
#
# Author: Eva Stare
#
# Usage: bash batch_predicTF.sh <predictf_path> <input_dir> <output_dir>

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PREDICTF_PATH="${1:?Usage: batch_predicTF.sh <predictf_path> <input_dir> <output_dir>}"
INPUT_DIR="${2:?Provide input directory containing .fa/.faa files}"
OUTPUT_BASE_DIR="${3:?Provide output directory}"
SCRIPT_PATH="$PREDICTF_PATH/scripts/predictf_in_genome.sh"

mkdir -p "$OUTPUT_BASE_DIR"

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

if [ ! -f "$SCRIPT_PATH" ]; then
    echo "ERROR: PredicTF script not found: $SCRIPT_PATH"
    exit 1
fi

# ---------------------------------------------------------------------------
# Process each genome
# ---------------------------------------------------------------------------
processed=0
failed=0

echo "Starting PredicTF batch analysis..."
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_BASE_DIR"
echo "=================================="

for input_file in "$INPUT_DIR"/*.fa "$INPUT_DIR"/*.faa; do
    [ ! -f "$input_file" ] && continue

    filename=$(basename "$input_file")
    basename_no_ext="${filename%.*}"
    output_dir="$OUTPUT_BASE_DIR/${basename_no_ext}_results"

    echo "Processing: $filename"
    mkdir -p "$output_dir"

    if sh "$SCRIPT_PATH" "$PREDICTF_PATH" "$input_file" "$output_dir"; then
        echo "  OK: $filename"
        ((processed++))

        if [ -f "$output_dir/file.out.mapping.TF" ]; then
            echo "  Confirmed TFs: $(wc -l < "$output_dir/file.out.mapping.TF")"
        fi
        if [ -f "$output_dir/file.out.mapping.potential.TF" ]; then
            echo "  Potential TFs: $(wc -l < "$output_dir/file.out.mapping.potential.TF")"
        fi
    else
        echo "  FAILED: $filename"
        ((failed++))
    fi

    echo "------------------------"
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo "=================================="
echo "Batch processing completed."
echo "Successfully processed: $processed files"
echo "Failed: $failed files"
echo "Results saved in: $OUTPUT_BASE_DIR"
