"""
S14_KD_group_validation.py - Kin discrimination group validation and subgroup generation
----------------------------------------------------------------------------------------
Author: Eva Stare

Validates kin discrimination (KD) group assignments against experimentally
determined swarming interaction phenotypes and resolves within-group
discrimination events by exclusion and subgroup reassignment.

Strains are treated as vertices and within-group discrimination events as
edges, so identifying the fewest strains whose exclusion removes all
within-group discrimination is a minimum vertex cover problem. Exact minimum
covers are obtained by exhaustive enumeration, which is tractable at the group
sizes encountered here (n <= 27). All equally optimal covers are enumerated and
exported. A greedy vertex-removal heuristic is run alongside for comparison but
does not determine the reported assignments.

Where the minimum vertex cover was not unique, one solution was fixed a priori so 
that group assignments are deterministic and reproducible across runs. The alternatives 
are equivalent by every criterion applied here — identical cover size, identical 
retained group size, and no within-group conflicts under either — so the choice among 
them is arbitrary and the sensitivity of downstream results to this choice is 
reported in Supplementary Results 3.

Excluded strains are partitioned into kin-consistent subgroups by maximal
clique enumeration under an explicit deterministic objective.

The solver comparison, the full set of equally optimal covers, and any
alternative subgroup partitions are reported to stdout at run time; the
written outputs are the group assignments and summary only.

Phenotype coding:
  1.0  = strong kin (merge)
  0.5  = weak kin (merge)
  0.25 = weak non-kin (discrimination)
  0.0  = strong non-kin (discrimination)
"""

import os
from collections import defaultdict
from itertools import combinations

import networkx as nx
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
INPUT_FILE = "data/processed/strain_interactions.tsv"
OUTPUT_DIR = "data/kd/output"
OUTPUT_FILE = OUTPUT_DIR + "/KD_Validation_with_Subgroups.xlsx"

COL_S1 = "S1"
COL_S2 = "S2"
COL_PHENOTYPE = "Phenotype_interaction"
COL_GROUP_S1 = "Group_S1"
COL_GROUP_S2 = "Group_S2"

KIN_VALUES = [1.0, 0.5]
DISC_VALUES = [0.25, 0.0]

# Published exclusion sets, retained where the minimum vertex cover is not
# unique. Both alternatives in g25 are provably minimum covers of equal size;
# this one is fixed by the original analysis and by Supplementary Results 3.
# Keys are group labels, values are strain names as they appear in INPUT_FILE.
PINNED_EXCLUSIONS = {
    "g25": ["NRS6107", "NRS6116", "NRS6190", "P8B3", "PS11", "PS51"],
}

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
df = pd.read_csv(INPUT_FILE, sep="\t")
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
# 3. Exact minimum exclusion (minimum vertex cover)
# ---------------------------------------------------------------------------
# Interaction lookup, used by the exclusion solver and the subgroup assignment.
interactions = {}
for _, row in df.iterrows():
    interactions[tuple(sorted([str(row[COL_S1]), str(row[COL_S2])]))] = row[COL_PHENOTYPE]


def get_interaction(s1, s2):
    """Return the phenotype score for a strain pair, or None if untested."""
    return interactions.get(tuple(sorted([str(s1), str(s2)])), None)


def edges_for_group(violations_df, group, col_s1, col_s2, col_group):
    """Return the set of within-group violation edges for one group."""
    g = violations_df[violations_df[col_group] == group]
    return {tuple(sorted((str(r[col_s1]), str(r[col_s2])))) for _, r in g.iterrows()}


def greedy_exclusion(edges):
    """Greedy vertex-removal heuristic.

    Iteratively removes the strain involved in the largest number of remaining
    violations, breaking ties alphabetically. Retained for comparison only:
    this is a heuristic and is not guaranteed to return a minimum cover.
    """
    remaining = set(edges)
    removed = []
    while remaining:
        scores = defaultdict(int)
        for s1, s2 in remaining:
            scores[s1] += 1
            scores[s2] += 1
        best = min(scores, key=lambda s: (-scores[s], s))
        removed.append(best)
        remaining = {(a, b) for a, b in remaining if a != best and b != best}
    return sorted(removed)


def exact_minimum_vertex_covers(edges):
    """Enumerate all minimum vertex covers by exhaustive search.

    Searches subsets of increasing size and returns the first size at which a
    valid cover exists, together with every cover of that size.
    Returns (cover_size, [frozenset, ...]).
    """
    if not edges:
        return 0, [frozenset()]
    vertices = sorted({v for e in edges for v in e})
    for k in range(1, len(vertices) + 1):
        covers = [frozenset(sub) for sub in combinations(vertices, k)
                  if all(a in set(sub) or b in set(sub) for a, b in edges)]
        if covers:
            return k, covers
    return len(vertices), [frozenset(vertices)]


def retained_kin_edges(cover, edges):
    """Number of kin interactions among strains retained under a cover.

    Used as the tie-break criterion when no pinned solution is defined: the
    cover leaving the most internally connected group is preferred.
    """
    members = {v for e in edges for v in e} - set(cover)
    return sum(1 for a, b in combinations(sorted(members), 2)
               if get_interaction(a, b) in KIN_VALUES)


def solve_exclusions(violations_df, col_s1, col_s2, col_group):
    """Solve the exact minimum exclusion problem for every group with violations.

    Returns (excluded_strains, alternatives) where excluded_strains maps
    group -> selected exclusion set and alternatives maps group -> all equally
    optimal exclusion sets. The greedy heuristic is evaluated for comparison
    and reported to stdout only.
    """
    excluded, alternatives = {}, {}
    n_optimal = 0

    for group in sorted(violations_df[col_group].unique()):
        edges = edges_for_group(violations_df, group, col_s1, col_s2, col_group)
        greedy = greedy_exclusion(edges)
        k_exact, covers = exact_minimum_vertex_covers(edges)

        if group in PINNED_EXCLUSIONS:
            chosen = sorted(PINNED_EXCLUSIONS[group])
            assert frozenset(chosen) in covers, (
                f"{group}: pinned exclusion set {chosen} is not a minimum "
                f"vertex cover of size {k_exact}. Input data may have changed."
            )
            selection = "pinned (published solution)"
        else:
            chosen = sorted(max(
                covers,
                key=lambda c: (retained_kin_edges(c, edges), tuple(sorted(c)))
            ))
            selection = "tie-break" if len(covers) > 1 else "unique"

        excluded[group] = chosen
        alternatives[group] = sorted(sorted(c) for c in covers)

        n_optimal += int(len(greedy) == k_exact)

        flag = "optimal" if len(greedy) == k_exact else "SUBOPTIMAL"
        print(f"  {group}: {len(edges)} violation edges")
        print(f"    exact minimum    = {k_exact}: {chosen}  [{selection}]")
        print(f"    greedy heuristic = {len(greedy)}: {greedy}  [{flag}]")
        if len(covers) > 1:
            print(f"    {len(covers)} equally optimal solutions:")
            for c in alternatives[group]:
                print(f"      - {c}")

    n_groups = len(excluded)
    n_nonunique = sum(1 for v in alternatives.values() if len(v) > 1)
    print(f"\n  Greedy matched the exact minimum in {n_optimal}/{n_groups} groups")
    print(f"  Groups with a non-unique minimum cover: {n_nonunique}")
    print(f"  Total strains excluded: {sum(len(v) for v in excluded.values())}")

    return excluded, alternatives


print("\nExact minimum exclusion (minimum vertex cover):")
excluded_strains, exclusion_alternatives = solve_exclusions(
    within_violations, COL_S1, COL_S2, COL_GROUP_S1)

# ---------------------------------------------------------------------------
# 4. Deterministic clique-based subgroup assignment
# ---------------------------------------------------------------------------
def find_subgroups_for_excluded(strains, group_name, max_partitions=200):
    """Partition excluded strains into non-overlapping kin-consistent subgroups.

    Maximal cliques of the kin graph may overlap, so greedy first-come
    assignment is order-dependent. Assignment is instead made under an explicit
    objective: maximise the number of strains placed into subgroups of size
    >= 2, breaking ties by largest-subgroup-first then alphabetical order.
    Returns (subgroups, singletons, alternative_partitions).
    """
    strains = sorted(set(map(str, strains)))
    if len(strains) < 2:
        return [], strains, []

    G = nx.Graph()
    G.add_nodes_from(strains)
    for s1, s2 in combinations(strains, 2):
        if get_interaction(s1, s2) in KIN_VALUES:
            G.add_edge(s1, s2)

    cliques = sorted([sorted(c) for c in nx.find_cliques(G) if len(c) >= 2],
                     key=lambda c: (-len(c), c))
    if not cliques:
        return [], strains, []

    best = {"covered": -1, "partitions": [], "keys": set()}

    def _search(idx, chosen, used):
        if idx == len(cliques):
            covered = len(used)
            key = tuple(sorted(tuple(c) for c in chosen))
            if covered > best["covered"]:
                best.update(covered=covered, partitions=[chosen[:]], keys={key})
            elif covered == best["covered"] and key not in best["keys"]:
                if len(best["partitions"]) < max_partitions:
                    best["partitions"].append(chosen[:])
                    best["keys"].add(key)
            return
        c = cliques[idx]
        if not set(c) & used:
            _search(idx + 1, chosen + [c], used | set(c))
        _search(idx + 1, chosen, used)

    _search(0, [], set())

    partitions = sorted(
        (sorted(p, key=lambda c: (-len(c), c)) for p in best["partitions"]),
        key=lambda p: [(-len(c), c) for c in p]
    )
    subgroups = partitions[0] if partitions else []
    assigned = {s for c in subgroups for s in c}
    singletons = [s for s in strains if s not in assigned]
    alternatives = partitions[1:]

    if alternatives:
        print(f"  {group_name}: {len(partitions)} equally optimal partitions "
              f"of the excluded strains")
        for p in partitions:
            print(f"    - {p}")

    return subgroups, singletons, alternatives


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
        subgroups, singletons, subgroup_alternatives = find_subgroups_for_excluded(
            sorted(excluded), group)

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
plt.savefig(os.path.join(OUTPUT_DIR, "all_kd_groups_matrices_new.png"),
            dpi=150, bbox_inches="tight")
print(f"Saved: {OUTPUT_DIR}/all_kd_groups_matrices_new.png")

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

# Plain-TSV copies of the two output sheets (canonical; xlsx is human copy)
assignments_df.to_csv(os.path.join(OUTPUT_DIR, "KD_Validation_Group_Assignments.tsv"),
                      sep="\t", index=False)
summary_data.to_csv(os.path.join(OUTPUT_DIR, "KD_Validation_Summary.tsv"),
                    sep="\t", index=False)

print(f"\nSaved to: {OUTPUT_FILE}")