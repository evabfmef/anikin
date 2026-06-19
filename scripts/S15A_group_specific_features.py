"""
S15A_group_specific_features.py — Identify KD Group-Specific Genomic Features
-----------------------------------------------------------------------------
Author: Eva Stare

For each kin discrimination (KD) group, identifies:
  - Features uniquely present in all group members and absent in all non-members
  - Features uniquely absent in all group members but present in all non-members

"""

import pandas as pd
import os

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
INPUT_FILE = "input/megamatrix_filtered_strains.xlsx"
OUTPUT_FILE = "output/features_analysis_results.xlsx"

STRAIN_GROUPS = {
    "g3":    ["RO-F-3"],
    "g4":    ["RS-D-2"],
    "g5":    ["RO-DD-2"],
    "g7":    ["RO-A-4"],
    "g19":   ["NRS6085"],
    "g25_A": ["PS-65", "NRS6202", "PS-14", "NCIB_3610", "P8_B1", "PS-233",
              "MB8_B7", "NRS6118", "PS-216", "P9_B1", "NRS6121", "PS-13",
              "PS-237", "PS-168", "PS-96", "PS-30", "NRS6153", "PS-68",
              "PS-18", "PS-210", "PS-31"],
    "g29":   ["BS16045"],
    "g30":   ["NRS6105", "NRS6145"],
    "g31":   ["PS-52", "PS-53"],
    "g34":   ["NRS6128"],
    "g39":   ["PS-209"],
    "g44":   ["PS-196"],
    "g50":   ["PS-108", "PS-109", "PS-131", "PS-119", "PS-130"],
    "g52":   ["NRS6160"],
    "g53":   ["PS-24", "PS-25", "PS-20"],
    "g54":   ["PS-160"],
    "g55":   ["PS-93", "PS-95"],
    "g56":   ["PS-217", "PS-218"],
    "g57":   ["PS-15"],
    "g60_A": ["NRS6181", "PS-194", "NRS6186", "PS-64"],
    "g60_B": ["NRS6132", "NRS6099", "NRS6187"],
    "g63_A": ["PS-55", "PS-149", "MB9_B1", "MB8_B1", "NRS6127", "MB9_B6",
              "NRS6103"],
    "g82":   ["RO-FF-1"],
    "g92":   ["KF24"],
    "g95":   ["PS-263"],
    "g96":   ["73"],
}

os.makedirs("output", exist_ok=True)

# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
df = pd.read_excel(INPUT_FILE, index_col=0)
print(f"Loaded: {df.shape[0]} strains x {df.shape[1]} features")

features_unique_to_group = {}
features_absent_in_group = {}

for group_name, strains in STRAIN_GROUPS.items():
    df_group = df.loc[strains]
    df_not_group = df.drop(strains)

    unique_features = []
    absent_features = []

    for col in df.columns:
        # Present in all group members, absent in all non-members
        if df_group[col].all() and not df_not_group[col].any():
            unique_features.append(col)
        # Absent in all group members, present in all non-members
        if not df_group[col].any() and df_not_group[col].all():
            absent_features.append(col)

    features_unique_to_group[group_name] = unique_features
    features_absent_in_group[group_name] = absent_features

    print(f"{group_name}: {len(unique_features)} unique, {len(absent_features)} absent")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
df_unique = pd.DataFrame.from_dict(features_unique_to_group, orient="index")
df_absent = pd.DataFrame.from_dict(features_absent_in_group, orient="index")

with pd.ExcelWriter(OUTPUT_FILE, engine="xlsxwriter") as writer:
    df_unique.to_excel(writer, sheet_name="Unique_Features")
    df_absent.to_excel(writer, sheet_name="Absent_Features")

print(f"\nSaved to: {OUTPUT_FILE}")
