"""
S14_KD_group_validation.py — Kin Discrimination Group Validation and Subgroup Generation
--------------------------------------------------------------------------------------
Author: Eva Stare

This script validates kin discrimination (KD) group assignments against
experimentally determined swarming interaction phenotypes and resolves
within-group violations through minimum-exclusion and subgroup reassignment.

Phenotype coding:
  1.0  = strong kin (merge)
  0.5  = weak kin (merge)
  0.25 = weak non-kin (discrimination)
  0.0  = strong non-kin (discrimination)
"""

import pandas as pd
import numpy as np
from collections import defaultdict
import networkx as nx
import os

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
INPUT_FILE = "input/Strain_interactions.xlsx"
OUTPUT_FILE = "output/KD_Validation_with_Subgroups.xlsx"

COL_S1 = "S1"
COL_S2 = "S2"
COL_PHENOTYPE = "Phenotype_interaction"
COL_GROUP_S1 = "Group_S1"
COL_GROUP_S2 = "Group_S2"

KIN_VALUES = [1.0, 0.5]
DISC_VALUES = [0.25, 0.0]

os.makedirs("output", exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
df = pd.read_excel(INPUT_FILE)
print(f"Loaded {len(df)} interactions")

# ---------------------------------------------------------------------------
# 2. Identify within-group violations
# ---------------------------------------------------------------------------
within_group = df[df[COL_GROUP_S1] == df[COL_GROUP_S2]].copy()
within_violations = within_group[within_group[COL_PHENOTYPE].isin(DISC_VALUES)].copy()

print(f"Within-group kin interactions: {len(within_group) - len(within_violations)}")
print(f"Within-group violations: {len(within_violations)}")

# Count violations per strain
strain_violation_counts = defaultdict(lambda: {"total": 0, "strong": 0, "mild": 0,
                                                "partners": [], "group": None})
for _, row in within_violations.iterrows():
    for s in [row[COL_S1], row[COL_S2]]:
        severity = "strong" if row[COL_PHENOTYPE] == 0 else "mild"
        strain_violation_counts[s]["total"] += 1
        strain_violation_counts[s][severity] += 1
        strain_violation_counts[s]["group"] = row[COL_GROUP_S1]

    strain_violation_counts[row[COL_S1]]["partners"].append(row[COL_S2])
    strain_violation_counts[row[COL_S2]]["partners"].append(row[COL_S1])

# ---------------------------------------------------------------------------
# 3. Greedy minimum-exclusion algorithm
# ---------------------------------------------------------------------------
def find_minimum_exclusion(violations_df, col_s1, col_s2, col_group):
    """
    For each group with violations, iteratively remove the strain with the
    most violations until no violations remain. Returns a dict of
    {group: [excluded strains]}.
    """
    results = {}

    for group in violations_df[col_group].unique():
        group_viols = violations_df[violations_df[col_group] == group]
        remaining = set()
        for _, row in group_viols.iterrows():
            remaining.add((row[col_s1], row[col_s2]))

        to_remove = []

        while remaining:
            strain_scores = defaultdict(int)
            for s1, s2 in remaining:
                strain_scores[s1] += 1
                strain_scores[s2] += 1

            worst = max(strain_scores, key=strain_scores.get)
            to_remove.append(worst)
            remaining = {(s1, s2) for s1, s2 in remaining
                         if s1 != worst and s2 != worst}

        results[group] = to_remove
        if to_remove:
            print(f"  {group}: exclude {to_remove}")

    return results


print("\nMinimum exclusion:")
excluded_strains = find_minimum_exclusion(within_violations, COL_S1, COL_S2, COL_GROUP_S1)

# ---------------------------------------------------------------------------
# 4. Reassign excluded strains via maximal clique detection
# ---------------------------------------------------------------------------
# Build interaction lookup
interactions = {}
for _, row in df.iterrows():
    key = tuple(sorted([row[COL_S1], row[COL_S2]]))
    interactions[key] = row[COL_PHENOTYPE]


def get_interaction(s1, s2):
    return interactions.get(tuple(sorted([s1, s2])), None)


def find_subgroups_for_excluded(strains, group_name):
    """
    Among excluded strains, build a kin-compatibility graph and find maximal
    cliques (size >= 2) as new subgroups. Remaining strains are singletons.
    """
    if len(strains) < 2:
        return [], list(strains)

    G = nx.Graph()
    G.add_nodes_from(strains)

    for i, s1 in enumerate(strains):
        for s2 in strains[i + 1:]:
            pheno = get_interaction(s1, s2)
            if pheno is not None and pheno in KIN_VALUES:
                G.add_edge(s1, s2)

    # Find all maximal cliques of size >= 2
    cliques = [c for c in nx.find_cliques(G) if len(c) >= 2]

    # Greedy assignment: largest cliques first, no strain assigned twice
    cliques.sort(key=len, reverse=True)
    assigned = set()
    subgroups = []

    for clique in cliques:
        available = [s for s in clique if s not in assigned]
        if len(available) >= 2:
            subgroups.append(available)
            assigned.update(available)

    singletons = [s for s in strains if s not in assigned]
    return subgroups, singletons


# Build strain-to-group mapping
strain_to_group = {}
for _, row in df.iterrows():
    strain_to_group[row[COL_S1]] = row[COL_GROUP_S1]
    strain_to_group[row[COL_S2]] = row[COL_GROUP_S2]

# All groups (including those without violations)
all_groups = sorted(set(strain_to_group.values()))

# Generate new assignments
new_assignments = []

print("\nSubgroup generation:")
for group in all_groups:
    all_in_group = sorted(s for s, g in strain_to_group.items() if g == group)

    if group in excluded_strains and excluded_strains[group]:
        excluded = set(excluded_strains[group])
        retained = [s for s in all_in_group if s not in excluded]

        # Retained strains keep the group as _A
        for s in retained:
            new_assignments.append({
                "Strain": s,
                "Original_Group": group,
                "New_Group": f"{group}_A",
                "Status": "Retained"
            })

        # Find subgroups among excluded strains
        subgroups, singletons = find_subgroups_for_excluded(list(excluded), group)

        for idx, sg in enumerate(subgroups):
            suffix = chr(ord("B") + idx)
            print(f"  {group}_{suffix}: {sorted(sg)}")
            for s in sg:
                new_assignments.append({
                    "Strain": s,
                    "Original_Group": group,
                    "New_Group": f"{group}_{suffix}",
                    "Status": "New_subgroup"
                })

        for s in singletons:
            print(f"  {group} singleton: {s}")
            new_assignments.append({
                "Strain": s,
                "Original_Group": group,
                "New_Group": f"{group}_{s}",
                "Status": "Singleton"
            })
    else:
        # No violations — group unchanged
        for s in all_in_group:
            new_assignments.append({
                "Strain": s,
                "Original_Group": group,
                "New_Group": group,
                "Status": "Unchanged"
            })

assignments_df = pd.DataFrame(new_assignments)

# ---------------------------------------------------------------------------
# 5. Validate new assignments
# ---------------------------------------------------------------------------
new_group_map = dict(zip(assignments_df["Strain"], assignments_df["New_Group"]))

remaining_violations = []
for _, row in df.iterrows():
    s1, s2 = row[COL_S1], row[COL_S2]
    pheno = row[COL_PHENOTYPE]
    new_g1 = new_group_map.get(s1)
    new_g2 = new_group_map.get(s2)

    if new_g1 and new_g2 and new_g1 == new_g2:
        if pheno in DISC_VALUES:
            remaining_violations.append((new_g1, s1, s2, pheno))

print(f"\nValidation: {len(remaining_violations)} violations remaining")
if remaining_violations:
    for group, s1, s2, pheno in remaining_violations:
        print(f"  {group}: {s1} vs {s2} = {pheno}")

# ---------------------------------------------------------------------------
# 6. Within-group interaction heatmaps (post-validation)
# ---------------------------------------------------------------------------
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap, BoundaryNorm

# Apply new group assignments to the dataframe for visualization
for idx_row, row in df.iterrows():
    if row[COL_S1] in new_group_map:
        df.loc[idx_row, COL_GROUP_S1] = new_group_map[row[COL_S1]]
    if row[COL_S2] in new_group_map:
        df.loc[idx_row, COL_GROUP_S2] = new_group_map[row[COL_S2]]


def extract_pheno(val):
    if pd.isna(val):
        return np.nan
    val_str = str(val).strip()
    if val_str in ["?"]:
        return np.nan
    try:
        return float(val_str.split()[0].replace(",", "."))
    except ValueError:
        return np.nan


same_group = df[df[COL_GROUP_S1] == df[COL_GROUP_S2]].copy()
same_group["pheno_numeric"] = same_group[COL_PHENOTYPE].apply(extract_pheno)

cmap = ListedColormap(["#d62728", "#ff7f0e", "#90EE90", "#2ca02c"])
norm = BoundaryNorm(boundaries=[0, 0.125, 0.375, 0.75, 1.0], ncolors=4)


def plot_group_matrix(grp, ax=None):
    """Plot within-group phenotype interaction matrix as a heatmap."""
    grp_data = same_group[same_group[COL_GROUP_S1] == grp].copy()
    if len(grp_data) == 0:
        return None
    strains = sorted(set(grp_data[COL_S1].astype(str)) | set(grp_data[COL_S2].astype(str)))
    matrix = pd.DataFrame(index=strains, columns=strains, dtype=float)
    for _, row in grp_data.iterrows():
        s1, s2, pheno = str(row[COL_S1]), str(row[COL_S2]), row["pheno_numeric"]
        matrix.loc[s1, s2] = pheno
        matrix.loc[s2, s1] = pheno
    for s in strains:
        matrix.loc[s, s] = 1.0
    pheno_vals = grp_data["pheno_numeric"].dropna()
    n_kin = (pheno_vals >= 0.5).sum()
    n_nonkin = (pheno_vals < 0.5).sum()
    if ax is None:
        fig, ax = plt.subplots(figsize=(max(8, len(strains) * 0.6),) * 2)
    fontsize = 10 if len(strains) <= 8 else (8 if len(strains) <= 15 else 6)
    sns.heatmap(matrix.astype(float), ax=ax, cmap=cmap, norm=norm,
                annot=True, fmt=".2f", annot_kws={"size": fontsize},
                cbar_kws={"label": "Phenotype (0,0.25=Nonkin | 0.5,1=Kin)"},
                linewidths=0.5, linecolor="gray")
    ax.set_title(f"{grp} ({len(strains)} strains)\nKin: {n_kin}, Nonkin: {n_nonkin}",
                 fontsize=12, fontweight="bold")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    plt.setp(ax.get_yticklabels(), rotation=0)
    return matrix


# Plot all groups (excluding singletons)
groups = sorted(same_group.groupby(COL_GROUP_S1).size()
                .pipe(lambda s: s[s >= 1]).index)
groups = [g for g in groups if "_singleton" not in g]

n_cols = 3
n_rows = (len(groups) + n_cols - 1) // n_cols
fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 6 * n_rows))
axes = axes.flatten()

for i, grp in enumerate(groups):
    plot_group_matrix(grp, ax=axes[i])
for i in range(len(groups), len(axes)):
    axes[i].axis("off")

plt.suptitle("Within-group Phenotype Interactions\n"
             "(Green = Kin [1.0, 0.5], Red/Orange = Nonkin [0.25, 0])",
             fontsize=14, fontweight="bold", y=1.02)
plt.tight_layout()
plt.savefig(os.path.join("output", "all_kd_groups_matrices_new.png"),
            dpi=150, bbox_inches="tight")
print("Saved: output/all_kd_groups_matrices_new.png")

# ---------------------------------------------------------------------------
# 7. Summary and export
# ---------------------------------------------------------------------------
print("\nSummary:")
for status, count in assignments_df["Status"].value_counts().items():
    print(f"  {status}: {count}")

print(f"\nOriginal groups: {len(all_groups)}")
print(f"New groups: {assignments_df['New_Group'].nunique()}")

# Save
summary_data = pd.DataFrame({
    "Metric": [
        "Total Strains", "Original Groups", "New Groups",
        "Retained", "New Subgroups", "Singletons",
        "Original Violations", "Remaining Violations"
    ],
    "Value": [
        len(assignments_df), len(all_groups),
        assignments_df["New_Group"].nunique(),
        len(assignments_df[assignments_df["Status"] == "Retained"]),
        len(assignments_df[assignments_df["Status"] == "New_subgroup"]),
        len(assignments_df[assignments_df["Status"] == "Singleton"]),
        len(within_violations), len(remaining_violations)
    ]
})

with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
    assignments_df.to_excel(writer, sheet_name="Group_Assignments", index=False)
    summary_data.to_excel(writer, sheet_name="Summary", index=False)

print(f"\nSaved to: {OUTPUT_FILE}")
