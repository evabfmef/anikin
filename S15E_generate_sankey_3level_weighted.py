"""
S15E_generate_sankey_3level_weighted.py

Generate three-level weighted Sankey diagram for group-specific gene features.

Visualizes the flow: KD groups → functional categories → feature category
(present and unique / present but additional / absent).

Genes with multiple functional categories are weight-adjusted so each gene
contributes a total weight of 1, split equally across its categories.

Author: Eva Stare
"""

import pandas as pd
import plotly.graph_objects as go

# =============================================================================
# 1. Read and validate data
# =============================================================================

INPUT_FILE = "39_KD_features_summary_9kat_p_PRESENT.csv"
REQUIRED_COLS = ["Kin group", "Gene name", "Functional category", "Category"]

df = pd.read_csv(INPUT_FILE, sep=";")
df = df.dropna(subset=REQUIRED_COLS).copy()

for col in REQUIRED_COLS:
    df[col] = df[col].astype(str).str.strip()

# Remove exact duplicate rows
df = df.drop_duplicates(subset=REQUIRED_COLS).copy()

print(f"Loaded {len(df)} gene–group–category entries from {INPUT_FILE}")

# =============================================================================
# 2. Split multi-category genes and compute weights
# =============================================================================

# A gene can belong to multiple functional categories (semicolon-separated).
# Each gene contributes total weight = 1, divided equally across its categories.
df["Functional category"] = df["Functional category"].str.split(r"\s*;\s*")
df["n_cat"] = df["Functional category"].apply(len)
df = df.explode("Functional category").copy()

df["Functional category"] = df["Functional category"].astype(str).str.strip()
df["Category"] = df["Category"].astype(str).str.strip()
df = df[(df["Functional category"] != "") & (df["Category"] != "")].copy()

df["weight"] = 1 / df["n_cat"]

# =============================================================================
# 3. Aggregate weighted links for both Sankey levels
# =============================================================================

# Left → Middle: KD group → Functional category
agg_left_mid = (
    df.groupby(["Kin group", "Functional category"], as_index=False)["weight"]
      .sum()
)

# Middle → Right: Functional category → Feature category
agg_mid_right = (
    df.groupby(["Functional category", "Category"], as_index=False)["weight"]
      .sum()
)

# =============================================================================
# 4. Define node ordering
# =============================================================================

preferred_kin_order = [
    "g1", "g2", "g3", "g4", "g5", "g6",
    "g7", "g8", "g9", "g10", "g11", "g12",
    "g13", "g14"
]

preferred_func_order = [
    "Toxin-antitoxin & competition",
    "Cell wall & membrane",
    "Metabolism",
    "Stress response & resistance",
    "DNA-related",
    "Regulation",
    "Transport",
    "Quorum sensing",
    "Other"
]

preferred_cat_order = [
    "present and unique",
    "present but additional",
    "absent"
]

# Build ordered node lists (preferred order first, then any extras alphabetically)
kin_present = set(agg_left_mid["Kin group"])
func_present = set(agg_left_mid["Functional category"]).union(
    set(agg_mid_right["Functional category"])
)
cat_present = set(agg_mid_right["Category"])

left_nodes = [x for x in preferred_kin_order if x in kin_present]
left_nodes += sorted(x for x in kin_present if x not in left_nodes)

mid_nodes = [x for x in preferred_func_order if x in func_present]
mid_nodes += sorted(x for x in func_present if x not in mid_nodes)

right_nodes = [x for x in preferred_cat_order if x in cat_present]
right_nodes += sorted(x for x in cat_present if x not in right_nodes)

labels = left_nodes + mid_nodes + right_nodes
label_to_index = {label: i for i, label in enumerate(labels)}

print(f"Nodes: {len(left_nodes)} KD groups, {len(mid_nodes)} functional categories, "
      f"{len(right_nodes)} feature categories")

# =============================================================================
# 5. Define color palettes
# =============================================================================

kin_colors = {
    "g1":  "rgba(99,110,250,0.85)",
    "g2":  "rgba(239,85,59,0.85)",
    "g3":  "rgba(0,204,150,0.85)",
    "g4":  "rgba(171,99,250,0.85)",
    "g5":  "rgba(255,161,90,0.85)",
    "g6":  "rgba(25,211,243,0.85)",
    "g7":  "rgba(255,102,146,0.85)",
    "g8":  "rgba(182,232,128,0.85)",
    "g9":  "rgba(255,151,255,0.85)",
    "g10": "rgba(254,203,82,0.85)",
    "g11": "rgba(120,120,120,0.85)",
    "g12": "rgba(50,50,50,0.85)",
}

func_colors = {
    "Toxin-antitoxin & competition": "rgba(240,200,80,0.85)",
    "Cell wall & membrane":          "rgba(70,200,230,0.85)",
    "Metabolism":                     "rgba(245,120,160,0.85)",
    "Stress response & resistance":   "rgba(255,140,100,0.85)",
    "DNA-related":                    "rgba(140,140,255,0.85)",
    "Regulation":                     "rgba(180,180,180,0.85)",
    "Transport":                      "rgba(120,180,255,0.85)",
    "Quorum sensing":                 "rgba(220,150,230,0.85)",
    "Other":                          "rgba(200,200,200,0.85)",
}

cat_colors = {
    "present and unique":    "rgba(60,180,75,0.85)",
    "present but additional": "rgba(255,140,0,0.85)",
    "absent":                "rgba(180,180,180,0.85)",
}

# Assign node colors
node_colors = []
for label in labels:
    if label in kin_colors:
        node_colors.append(kin_colors[label])
    elif label in func_colors:
        node_colors.append(func_colors[label])
    elif label in cat_colors:
        node_colors.append(cat_colors[label])
    else:
        node_colors.append("rgba(180,180,180,0.85)")

# =============================================================================
# 6. Build Sankey links
# =============================================================================

# Left → Middle links (colored by KD group, semi-transparent)
source_1 = agg_left_mid["Kin group"].map(label_to_index).tolist()
target_1 = agg_left_mid["Functional category"].map(label_to_index).tolist()
value_1 = agg_left_mid["weight"].tolist()
link_colors_1 = [
    kin_colors.get(kg, "rgba(160,160,160,0.45)").replace("0.85", "0.45")
    for kg in agg_left_mid["Kin group"]
]

# Middle → Right links (colored by functional category, semi-transparent)
source_2 = agg_mid_right["Functional category"].map(label_to_index).tolist()
target_2 = agg_mid_right["Category"].map(label_to_index).tolist()
value_2 = agg_mid_right["weight"].tolist()
link_colors_2 = [
    func_colors.get(fc, "rgba(160,160,160,0.45)").replace("0.85", "0.45")
    for fc in agg_mid_right["Functional category"]
]

# Merge all links
source = source_1 + source_2
target = target_1 + target_2
value = value_1 + value_2
link_colors = link_colors_1 + link_colors_2

print(f"Links: {len(source_1)} (groups→functions) + {len(source_2)} (functions→categories) "
      f"= {len(source)} total")

# =============================================================================
# 7. Create Sankey diagram
# =============================================================================

fig = go.Figure(data=[go.Sankey(
    arrangement="snap",
    node=dict(
        pad=20,
        thickness=18,
        line=dict(color="rgba(50,50,50,0.6)", width=0.4),
        label=labels,
        color=node_colors,
    ),
    link=dict(
        source=source,
        target=target,
        value=value,
        color=link_colors,
    ),
)])

fig.update_layout(
    title_text="Group-specific genes: kin groups → functional categories → category",
    font=dict(size=13),
    width=1500,
    height=900,
    paper_bgcolor="white",
    plot_bgcolor="white",
)

# =============================================================================
# 8. Export
# =============================================================================

fig.write_html("sankey_3level_weighted39.html")
print("Saved: sankey_3level_weighted39.html")

fig.write_image("sankey_3level_weighted39.png", width=1800, height=1100, scale=3)
print("Saved: sankey_3level_weighted39.png (5400×3300 px)")

fig.write_image("sankey_3level_weighted39.svg")
print("Saved: sankey_3level_weighted39.svg")

fig.write_image("sankey_3level_weighted39.pdf")
print("Saved: sankey_3level_weighted39.pdf")

fig.show()
