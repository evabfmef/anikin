"""
S10B_GEM_analysis_part2.py — Draft Strain-Specific Metabolic Model Generation
--------------------------------------------------------------------------
Author: Eva Stare

Adapted from the multi-strain GEM generation workflow by Norsigian et al. (2020),
originally developed for E. coli. We applied this workflow to B. subtilis using
the reference genome-scale metabolic model iBB1018 (Blázquez et al., 2023),
which is distributed in SBML (.xml) format rather than JSON. The script was
modified accordingly. Additional minor modifications were made to improve
code organization and readability.

References:
  Norsigian CJ, Fang X, Seif Y, Monk JM, Palsson BO. A workflow for generating
  multi-strain genome-scale metabolic models of prokaryotes. Nat Protoc. 2020;
  15(1):1-14. doi:10.1038/s41596-019-0254-3

  Blázquez B, San León D, Rojas A, Tortajada M, Nogales J. New Insights on
  Metabolic Features of Bacillus subtilis Based on Multistrain Genome-Scale
  Metabolic Modeling. Int J Mol Sci. 2023;24:7091. doi:10.3390/ijms24087091

"""

import os
import sys
import cobra
import pandas as pd
import numpy as np
from glob import glob
from cobra.manipulation.delete import remove_genes
from cobra.manipulation.modify import rename_genes

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
WORK_DIR = "./GEM"
MODELS_DIR = os.path.join(WORK_DIR, "Models")
REFERENCE_MODEL = os.path.join(WORK_DIR, "iBB1018.xml")
ORTHO_MATRIX = os.path.join(WORK_DIR, "ortho_matrix.csv")
GENEIDS_MATRIX = os.path.join(WORK_DIR, "geneIDs_matrix.csv")

# Artificial genes to keep (essential for spontaneous reactions and gap-filling)
ARTIFICIAL_GENES = ["GROWMATCH", "GAPFILLING"]

os.makedirs(MODELS_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load reference model and homology matrix
# ---------------------------------------------------------------------------
print("Loading reference model...")
model = cobra.io.read_sbml_model(REFERENCE_MODEL)
print(f"Reference: {len(model.genes)} genes, {len(model.reactions)} reactions")

hom_matrix = pd.read_csv(ORTHO_MATRIX, index_col=0)
print(f"Homology matrix: {hom_matrix.shape}")

# ---------------------------------------------------------------------------
# 2. Generate draft strain-specific models
# ---------------------------------------------------------------------------
successful, failed = [], []

for strain in hom_matrix.columns:
    try:
        # Create model copy via JSON intermediate (avoids deepcopy issues)
        temp_file = os.path.join(WORK_DIR, f"temp_{strain}.json")
        cobra.io.json.save_json_model(model, temp_file, pretty=False)
        model_copy = cobra.io.load_json_model(temp_file)
        os.remove(temp_file)

        # Identify genes absent in this strain
        absent = hom_matrix[strain][hom_matrix[strain] == 0.0].index.tolist()
        genes_to_remove = [g for g in absent if g not in ARTIFICIAL_GENES]

        # Get Gene objects and remove
        to_delete = []
        for gene_id in genes_to_remove:
            try:
                to_delete.append(model_copy.genes.get_by_id(gene_id))
            except KeyError:
                pass

        if to_delete:
            remove_genes(model_copy, to_delete, remove_reactions=True)

        model_copy.id = str(strain)
        output_path = os.path.join(MODELS_DIR, f"{strain}.json")
        cobra.io.json.save_json_model(model_copy, output_path, pretty=False)

        print(f"{strain}: {len(model_copy.genes)} genes, "
              f"{len(model_copy.reactions)} reactions")
        successful.append(strain)

    except Exception as e:
        print(f"Error: {strain} — {e}")
        failed.append(strain)

print(f"\nModels created: {len(successful)}, failed: {len(failed)}")

# ---------------------------------------------------------------------------
# 3. Update gene-protein-reaction rules with strain-specific locus tags
# ---------------------------------------------------------------------------
geneIDs_matrix = pd.read_csv(GENEIDS_MATRIX, index_col=0)
model_files = glob(os.path.join(MODELS_DIR, "*.json"))

for mod_path in model_files:
    mod = cobra.io.load_json_model(mod_path)

    # Match model file to strain column
    current_strain = None
    for col in geneIDs_matrix.columns:
        if col in mod_path:
            current_strain = col
            break

    if current_strain is None:
        continue

    id_mapping = geneIDs_matrix[current_strain].to_dict()
    id_mapping = {k: v for k, v in id_mapping.items() if v != "None"}

    rename_genes(mod, id_mapping)
    cobra.io.json.save_json_model(mod, mod_path, pretty=False)

print("\nGPR rules updated for all models.")

# ---------------------------------------------------------------------------
# 4. Summary of draft models
# ---------------------------------------------------------------------------
print("\nDraft model summary:")
for strain in hom_matrix.columns:
    path = os.path.join(MODELS_DIR, f"{strain}.json")
    if os.path.exists(path):
        m = cobra.io.load_json_model(path)
        print(f"  {m.id}: {len(m.genes)} genes, {len(m.reactions)} reactions")
