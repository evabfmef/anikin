#!/bin/bash
# =============================================================================
# S1_Prokka_annotation.sh
# Genome annotation with Prokka using a custom reference database.
# Author: Eva Stare
# =============================================================================

set -euo pipefail
shopt -s nullglob

files=(data/raw/fasta_genomes/*.fasta)
(( ${#files[@]} )) || { echo "no FASTA files found"; exit 1; }

mkdir -p annot_output

for file in "${files[@]}"; do
  base=$(basename "$file" .fasta)
  outdir="annot_output/$base"
  prokka "$file" \
    --cpus 16 \
    --outdir "$outdir" \
    --prefix "$base" \
    --increment 1 \
    --gffver 3 \
    --kingdom Bacteria \
    --gcode 11 \
    --proteins data/misc/manual_database.fasta \
    --evalue 1e-06 \
    --addgenes \
    --locustag "$base" \
    --force
done
