"""
S10C_GEM_analysis_part3.py — Reaction Presence/Absence Matrix Generation
--------------------------------------------------------------------
Author: Eva Stare

Generates a binary reaction presence/absence matrix from strain-specific
genome-scale metabolic models (GEMs) produced by S10B_GEM_analysis_part2.py.

"""

import cobra
from glob import glob
import pandas as pd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
MODELS_DIR = "./Models"
OUTPUT_FILE = "reaction_matrix.xlsx"

# ---------------------------------------------------------------------------
# 1. Load models and collect reactions
# ---------------------------------------------------------------------------
model_files = glob(f"{MODELS_DIR}/*.json")
print(f"Found {len(model_files)} model files")

all_reactions = set()
model_reactions = {}

for path in model_files:
    model = cobra.io.load_json_model(path)
    strain_name = path.split("/")[-1].replace(".json", "")
    rxns = {r.id for r in model.reactions}
    model_reactions[strain_name] = rxns
    all_reactions.update(rxns)

unique_reactions = sorted(all_reactions)
strain_names = sorted(model_reactions.keys())

print(f"Strains: {len(strain_names)}")
print(f"Unique reactions: {len(unique_reactions)}")

# ---------------------------------------------------------------------------
# 2. Build binary matrix
# ---------------------------------------------------------------------------
reaction_matrix = pd.DataFrame(0, index=strain_names, columns=unique_reactions)

for strain, rxns in model_reactions.items():
    for rxn in rxns:
        reaction_matrix.loc[strain, rxn] = 1

print(f"Reaction matrix shape: {reaction_matrix.shape}")

# ---------------------------------------------------------------------------
# 3. Save output
# ---------------------------------------------------------------------------
reaction_matrix.to_excel(OUTPUT_FILE)
print(f"Saved to {OUTPUT_FILE}")
