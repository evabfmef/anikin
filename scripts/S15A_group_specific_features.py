"""
S15A_group_specific_features.py — Identify KD Group-Specific Genomic Features
-----------------------------------------------------------------------------
Author: Eva Stare

For each kin discrimination (KD) group, identifies:
  - Features uniquely present in all group members and absent in all non-members
  - Features uniquely absent in all group members but present in all non-members

This is a descriptive exact-pattern criterion, not a statistical test: a feature
qualifies only under perfect co-segregation with group membership. Because the
criterion involves no test of significance, its interpretation depends on group
size. In a multi-strain group, perfect co-segregation across all members is
evidence of a shared determinant. In a single-strain group the criterion reduces
to identifying that strain's private accessory content, and in two-strain groups
it is satisfied readily by chance; features from such groups are reported for
completeness but are strain-specific rather than group-associated, and are not
treated as equivalent evidence. Groups are therefore annotated by size class
(SINGLETON, PAIR, MULTI) in a separate summary file (group_size_summary.xlsx /
group_summary.tsv); the feature workbook consumed by S15B-S15D is unchanged. 
Statistically supported associations are assessed separately by the pairwise 
enrichment analysis, which operates at the interaction level rather than the group level.

The input matrix is validated to contain only 0 and 1 before analysis. This
matters because pandas treats NaN as truthy, so an all-missing column would
otherwise satisfy .all() and be reported as uniquely present.
"""

import os

import pandas as pd

from kd_config import load_config, load_kd_groups, describe

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
#   python S15A_group_specific_features.py [--config=FILE] [--analysis=39|67]
cfg = load_config()
ANALYSIS = cfg["analysis"]

INPUT_FILE = (cfg["s15a"]["input_file_67"] if ANALYSIS == "67"
              else cfg["s15a"]["input_file_39"])
OUTPUT_DIR   = cfg["s15a"]["output_dir"]
OUTPUT_FILE  = os.path.join(OUTPUT_DIR, cfg["s15a"]["output_file"])
# Written separately: S15B and S15C parse every sheet of OUTPUT_FILE as a KD group.
SUMMARY_FILE = os.path.join(OUTPUT_DIR, cfg["s15a"]["summary_file"])

# Group memberships from metadata/ (version-controlled TSV).
STRAIN_GROUPS = load_kd_groups(cfg)

# Size-class boundaries; the exact-pattern criterion is descriptive and its
# interpretation depends on group size (see module docstring).
SINGLETON_MAX = cfg["s15a"]["size_classes"]["singleton_max"]
PAIR_MAX      = cfg["s15a"]["size_classes"]["pair_max"]

os.makedirs(OUTPUT_DIR, exist_ok=True)
print(describe(cfg))

# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
df = pd.read_csv(INPUT_FILE, sep="\t", index_col=0)
print(f"Loaded: {df.shape[0]} strains x {df.shape[1]} features")

# --- Validate that the matrix is strictly binary -------------------------
# Fail immediately on missing or non-binary values: pandas evaluates NaN as
# truthy, so .all() on an all-missing column would return True and the feature
# would be misreported as uniquely present.
n_missing = int(df.isna().sum().sum())
if n_missing:
    bad_cols = df.columns[df.isna().any()].tolist()
    raise ValueError(
        f"Matrix contains {n_missing} missing value(s) in {len(bad_cols)} "
        f"column(s), e.g. {bad_cols[:5]}. Presence/absence must be complete."
    )

observed = set(pd.unique(df.values.ravel()))
non_binary = observed - {0, 1}
if non_binary:
    raise ValueError(
        f"Matrix contains non-binary value(s): {sorted(non_binary)[:10]}. "
        f"Only 0 and 1 are permitted."
    )

print(f"Matrix validated: strictly binary, no missing values")

# --- Group size classes --------------------------------------------------
# The exact-pattern criterion is descriptive and its interpretation depends on
# group size; see module docstring.
def size_class(n):
    """Classify a KD group by the number of member strains."""
    if n <= SINGLETON_MAX:
        return "SINGLETON"
    if n <= PAIR_MAX:
        return "PAIR"
    return "MULTI"


group_sizes = {g: len(s) for g, s in STRAIN_GROUPS.items()}
group_classes = {g: size_class(n) for g, n in group_sizes.items()}

n_by_class = pd.Series(group_classes).value_counts().to_dict()
print(f"KD groups by size class: {n_by_class}")

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

    print(f"{group_name} [{group_classes[group_name]}, n={group_sizes[group_name]}]: "
          f"{len(unique_features)} unique, {len(absent_features)} absent")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
df_unique = pd.DataFrame.from_dict(features_unique_to_group, orient="index")
df_absent = pd.DataFrame.from_dict(features_absent_in_group, orient="index")

# Per-group summary annotated by size class. Written to a separate file rather
# than added as a sheet to OUTPUT_FILE: scripts S15B and S15C iterate over every
# sheet of that workbook and treat each as a KD group, so an extra sheet would
# be parsed as a phantom group.
df_summary = pd.DataFrame({
    "kd_group": list(STRAIN_GROUPS),
    "n_strains": [group_sizes[g] for g in STRAIN_GROUPS],
    "size_class": [group_classes[g] for g in STRAIN_GROUPS],
    "n_unique_features": [len(features_unique_to_group[g]) for g in STRAIN_GROUPS],
    "n_absent_features": [len(features_absent_in_group[g]) for g in STRAIN_GROUPS],
}).sort_values(["size_class", "kd_group"]).reset_index(drop=True)

with pd.ExcelWriter(OUTPUT_FILE, engine="xlsxwriter") as writer:
    df_unique.to_excel(writer, sheet_name="Unique_Features")
    df_absent.to_excel(writer, sheet_name="Absent_Features")

with pd.ExcelWriter(SUMMARY_FILE, engine="xlsxwriter") as writer:
    df_summary.to_excel(writer, sheet_name="Group_Summary", index=False)

# Plain-TSV copies (canonical; xlsx is the human copy)
df_unique.to_csv(os.path.join(OUTPUT_DIR, "features_unique_to_group.tsv"), sep="\t")
df_absent.to_csv(os.path.join(OUTPUT_DIR, "features_absent_in_group.tsv"), sep="\t")
df_summary.to_csv(os.path.join(OUTPUT_DIR, "group_summary.tsv"), sep="\t", index=False)

print("\nFeatures by group size class:")
for cls in ["MULTI", "PAIR", "SINGLETON"]:
    sub = df_summary[df_summary["size_class"] == cls]
    if len(sub):
        print(f"  {cls}: {len(sub)} group(s), "
              f"{int(sub['n_unique_features'].sum())} unique, "
              f"{int(sub['n_absent_features'].sum())} absent")

print(f"\nSaved to: {OUTPUT_FILE}")
print(f"Group size summary: {SUMMARY_FILE}")