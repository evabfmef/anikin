#!/bin/bash
# =============================================================================
# S1_Prokka_annotation.sh
# Genome annotation with Prokka using a custom reference database.
# Author: Eva Stare
# =============================================================================

# Output directory
mkdir -p annot_output_e06

# Loop through all .fasta files in input
for file in input_directory/*.fasta; do
  # Extract base filename without path or extension
  base=$(basename "$file" .fasta)

  # Output directory and prefix
  outdir="annot_output_e06/$base"
  prefix="$base"

  # Run Prokka with manual database and evalue 1e-06 (default setting with manual database)
  prokka "$file" \
    --cpus 16 \
    --outdir "$outdir" \
    --prefix "$prefix" \
    --increment 1 \
    --gffver 3 \
    --kingdom Bacteria \
    --gcode 11 \
    --proteins input_manual_db/manual_database.fasta \
    --evalue 1e-06 \
    --addgenes \
    --locustag CUSTOM
done
