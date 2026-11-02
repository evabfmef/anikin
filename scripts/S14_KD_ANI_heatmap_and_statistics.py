"""
S14_KD_ANI_heatmap_and_statistics.py - KD phenotype heatmap and ANI-KD statistics

Relates genome-wide relatedness (average nucleotide identity, ANI) to kin
discrimination (KD) phenotypes for the 850 experimentally tested strain pairs
(78 strains). Produces the KD phenotype heatmap (Fig. 3B) and the statistical
analyses described in Supplementary Methods 8.

Analyses
--------
1. Heatmap. Strains are ordered by UPGMA clustering of ANI distances
   (1 - ANI/100). Cells show the KD phenotype score of each tested pair
   (0, 0.25, 0.5, 1); untested pairs are grey. Diagonal cells show
   self-pairings (kin positive control, KD = 1), which are excluded from all
   statistical analyses. Side bars show candidate KD group membership
   (Supplementary Table S9) and the dashed line on the dendrogram marks the
   99.5% ANI threshold.

2. Pearson correlation between pairwise ANI and KD score (descriptive).

3. Strain-label permutation test (exploratory). Strain labels are permuted
   on the ANI matrix while the 850 tested pairs and their phenotypes are kept
   fixed (9,999 permutations). Because the tested pairs are not a random
   subset of all 3,003 possible pairs, exchangeability may be violated; this
   is therefore not a conventional Mantel test and is reported as exploratory.

4. Bayesian ordinal mixed-effects model (primary inference). KD category
   (0 < 0.25 < 0.5 < 1) is modelled as an ordered response with a cumulative
   logit link, mean-centred ANI (percentage points) as predictor, and one
   random intercept per strain that enters every pair in which the strain
   participates, irrespective of its position in the pair:
       eta = beta * ANI + u_i + u_j,   u ~ N(0, sigma^2)
   Priors: beta ~ N(0, 2); sigma ~ Half-Normal(1);
           ordered cutpoints ~ N((-2, 0, 2), 2).
   Sampling: 2 chains x 1,000 draws after 1,000 tuning steps,
   target acceptance 0.95. Convergence is checked by R-hat (< 1.01), bulk
   effective sample size (> 400) and divergent transitions. A posterior
   predictive check compares observed and simulated KD category frequencies
   (marginal frequencies only; it does not test the proportional-odds
   assumption). Posterior summaries use 95% highest-density intervals.

Inputs
------
Supplementary Tables workbook (--input):
  sheet "5"  - pairwise averaged reciprocal ANI matrix, 313 strains
  sheet "9"  - candidate KD group composition (one group per column)
  sheet "10" - pairwise swarming phenotypes and ANI values (850 pairs)

Outputs (--outdir)
------------------
  KD_heatmap_ANI_candidate_groups.png / .pdf / .svg
  ANI_KD_permutation_9999_null.csv
  ANI_KD_ordinal_model_summary.csv
  ANI_KD_posterior_predictive_frequencies.csv
  ANI_KD_ordinal_model_trace.nc

Usage
-----
  python S14_KD_ANI_heatmap_and_statistics.py \
      --input "Supplementary Tables isme.xlsx" --outdir output/S14

Environment
-----------
Python 3 with pandas, numpy, scipy, matplotlib, openpyxl, pymc and arviz.
Run in Google Colab for the published analysis.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="KD heatmap and ANI-KD statistics (S14).")
parser.add_argument("--input", default="Supplementary_Tables_isme.xlsx",
                    help="Supplementary Tables workbook (sheets 5, 9 and 10).")
parser.add_argument("--outdir", default="output/S14", help="Output directory.")
parser.add_argument("--skip-model", action="store_true",
                    help="Skip the Bayesian model (heatmap and permutation test only).")
args = parser.parse_args()

INPUT_FILE = args.input
OUTPUT_DIR = args.outdir
RUN_ORDINAL_MODEL = not args.skip_model

if not os.path.isfile(INPUT_FILE):
    raise FileNotFoundError(INPUT_FILE)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Expected panel size (fails loudly if the input does not match the manuscript)
N_STRAINS_EXPECTED = 78
N_PAIRS_EXPECTED = 850
N_GROUPS_EXPECTED = 25
N_ANI_STRAINS = 313

# ANI threshold shown on the dendrogram
ANI_THRESHOLD = 99.5
DISTANCE_THRESHOLD = 1 - ANI_THRESHOLD / 100

# Statistics
SEED = 20261020
N_PERM = 9999

# Figure
FIG_WIDTH, FIG_HEIGHT = 21, 19
STRAIN_FONT = 7.5
GROUP_FONT = 7
SHOW_GROUP_LABELS = False   # print group names above the colour strip
SPECIES = r"$\it{Bacillus\ subtilis}$"

KD_SCORES = [0, 0.25, 0.5, 1]
KD_COLORS = ["#B2182B", "#EF8A62", "#92C5DE", "#2166AC"]
KD_LABELS = ["0 – Strong non-kin", "0.25 – Weak non-kin",
             "0.5 – Weak kin", "1 – Strong kin"]
UNTESTED_COLOR = "#E6E6E6"
CONTROL_COLOR = "#08306B"      # diagonal: self-pairing, kin positive control

print(f"Input:  {INPUT_FILE}")
print(f"Output: {OUTPUT_DIR}")


def normalise_name(name):
    """Harmonise strain identifiers across sheets."""
    name = str(name).strip()
    return "NCIB 3610" if name in ("NCIB 3610", "NCIB_3610") else name


# ---------------------------------------------------------------------------
# 1. KD phenotypes (Supplementary Table S10)
# ---------------------------------------------------------------------------
df = pd.read_excel(INPUT_FILE, sheet_name="10", header=2)
df.columns = df.columns.astype(str).str.strip()

required = ["S1", "S2", "Phenotype", "ANI value"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns in sheet 10: {missing}")

df["Phenotype"] = pd.to_numeric(df["Phenotype"], errors="coerce")
df["ANI value"] = pd.to_numeric(df["ANI value"], errors="coerce")
df = df.dropna(subset=required).copy()
df["S1"] = df["S1"].apply(normalise_name)
df["S2"] = df["S2"].apply(normalise_name)

strains = sorted(set(df["S1"]) | set(df["S2"]))
n = len(strains)
print(f"Strains: {n}; tested pairs: {len(df)}")

if n != N_STRAINS_EXPECTED:
    raise ValueError(f"Expected {N_STRAINS_EXPECTED} strains, found {n}")
if len(df) != N_PAIRS_EXPECTED:
    raise ValueError(f"Expected {N_PAIRS_EXPECTED} interactions, found {len(df)}")
if not set(df["Phenotype"].unique()).issubset(set(KD_SCORES)):
    raise ValueError("Unexpected KD phenotype values")

kd_matrix = pd.DataFrame(np.nan, index=strains, columns=strains)
seen_pairs = set()
for _, row in df.iterrows():
    s1, s2, score = row["S1"], row["S2"], float(row["Phenotype"])
    if s1 == s2:
        raise ValueError(f"Self-comparison: {s1}")
    pair = tuple(sorted([s1, s2]))
    if pair in seen_pairs:
        raise ValueError(f"Duplicate pair: {pair}")
    seen_pairs.add(pair)
    kd_matrix.loc[s1, s2] = score
    kd_matrix.loc[s2, s1] = score

# ---------------------------------------------------------------------------
# 2. ANI matrix (Supplementary Table S5)
# ---------------------------------------------------------------------------
ani_raw = pd.read_excel(INPUT_FILE, sheet_name="5", header=None, skiprows=2)
ani_raw = ani_raw.iloc[:N_ANI_STRAINS].reset_index(drop=True)

ani_strains = [normalise_name(x) for x in ani_raw.iloc[:, 0]]
if len(set(ani_strains)) != len(ani_strains):
    raise ValueError("Duplicate strain names in sheet 5")

# Lower-triangular matrix; decimal commas are converted if present
ani_values = ani_raw.iloc[:, 1:].map(
    lambda x: x.replace(",", ".") if isinstance(x, str) else x)
ani_values = ani_values.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

# Rebuild the full symmetric matrix from the lower triangle
ani_array = np.full((len(ani_strains), len(ani_strains)), np.nan)
np.fill_diagonal(ani_array, 100.0)
for i in range(len(ani_strains)):
    for j in range(i):
        value = ani_values[i, j]
        if not np.isnan(value):
            ani_array[i, j] = value
            ani_array[j, i] = value
ani_complete = pd.DataFrame(ani_array, index=ani_strains, columns=ani_strains)

missing = [s for s in strains if s not in ani_complete.index]
if missing:
    raise ValueError(f"Strains missing from sheet 5: {missing}")

ani_panel = ani_complete.loc[strains, strains].copy()
if ani_panel.isna().any().any():
    raise ValueError("Missing ANI values in the extracted matrix")
if not np.allclose(ani_panel.values, ani_panel.values.T):
    raise ValueError("ANI matrix is not symmetric")
if ani_panel.values.min() < 0 or ani_panel.values.max() > 100:
    raise ValueError("ANI values outside 0–100")

# ---------------------------------------------------------------------------
# 3. Candidate KD groups (Supplementary Table S9)
# ---------------------------------------------------------------------------
# Group labels are in Excel row 3; each column lists the strains of one group.
groups_raw = pd.read_excel(INPUT_FILE, sheet_name="9", header=2)

group_assignment = {}
for column in groups_raw.columns:
    group_name = str(column).strip()
    if not group_name.startswith("G"):
        continue
    for value in groups_raw[column].dropna():
        strain = normalise_name(value)
        if strain in ("", "nan"):
            continue
        if strain in group_assignment and group_assignment[strain] != group_name:
            raise ValueError(f"{strain} assigned to more than one group")
        group_assignment[strain] = group_name

missing = [s for s in strains if s not in group_assignment]
if missing:
    raise ValueError(f"Missing candidate KD group assignments: {missing}")

groups = sorted({group_assignment[s] for s in strains}, key=lambda g: int(g[1:]))
print(f"Candidate KD groups: {len(groups)}")
if len(groups) != N_GROUPS_EXPECTED:
    raise ValueError(f"Expected {N_GROUPS_EXPECTED} groups, found {len(groups)}")

# ---------------------------------------------------------------------------
# 4. UPGMA clustering on ANI distance
# ---------------------------------------------------------------------------
ani_distance = 1 - ani_panel.to_numpy(dtype=float) / 100
np.fill_diagonal(ani_distance, 0)
Z = linkage(squareform(ani_distance, checks=True), method="average")

order = dendrogram(Z, no_plot=True)["leaves"]
ordered_strains = [strains[i] for i in order]
kd_ordered = kd_matrix.loc[ordered_strains, ordered_strains]
ordered_groups = [group_assignment[s] for s in ordered_strains]

# ---------------------------------------------------------------------------
# 5. Heatmap (Fig. 3B)
# ---------------------------------------------------------------------------
palette = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
group_colors = {g: palette[i] for i, g in enumerate(groups)}
group_cmap = ListedColormap([group_colors[g] for g in groups])
group_numbers = np.array([groups.index(g) for g in ordered_groups])

kd_cmap = ListedColormap(KD_COLORS)
kd_cmap.set_bad(UNTESTED_COLOR)
kd_norm = BoundaryNorm([-0.1, 0.125, 0.375, 0.75, 1.1], kd_cmap.N)

fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT), facecolor="white")
# Rows: dendrogram, group colour strip, heatmap.
# Columns: left group strip, heatmap, legends.
gs = GridSpec(3, 3, figure=fig,
              width_ratios=[0.025, 1, 0.28], height_ratios=[0.20, 0.025, 1],
              left=0.16, right=0.97, bottom=0.09, top=0.94,
              wspace=0.015, hspace=0.015)
ax_dendro = fig.add_subplot(gs[0, 1])
ax_col_colors = fig.add_subplot(gs[1, 1])
ax_row_colors = fig.add_subplot(gs[2, 0])
ax_heatmap = fig.add_subplot(gs[2, 1])
ax_legend = fig.add_subplot(gs[2, 2])
ax_legend.axis("off")

# Dendrogram with the 99.5% ANI reference line
dendrogram(Z, ax=ax_dendro, no_labels=True, color_threshold=0,
           above_threshold_color="#444444", link_color_func=lambda k: "#444444")
ax_dendro.set_xticks([])
ax_dendro.set_ylabel("ANI distance\n(1 − ANI/100)", fontsize=9)
ax_dendro.tick_params(axis="y", labelsize=8)
ax_dendro.axhline(DISTANCE_THRESHOLD, color="#B2182B", linestyle="--", linewidth=1)
ax_dendro.set_xlim(0, n * 10)

# Candidate KD group colour strips (top and left)
strip_kwargs = dict(aspect="auto", interpolation="nearest", cmap=group_cmap,
                    vmin=-0.5, vmax=len(groups) - 0.5)
ax_col_colors.imshow(group_numbers[np.newaxis, :], **strip_kwargs)
ax_col_colors.set_xlim(-0.5, n - 0.5)
ax_row_colors.imshow(group_numbers[:, np.newaxis], **strip_kwargs)
ax_row_colors.set_ylim(n - 0.5, -0.5)
for ax in (ax_col_colors, ax_row_colors):
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

# KD phenotype matrix
ax_heatmap.imshow(np.ma.masked_invalid(kd_ordered.values.astype(float)),
                  cmap=kd_cmap, norm=kd_norm, interpolation="nearest", aspect="equal")
# Diagonal: self-pairings (kin positive control), not part of the 850 pairs
diagonal = np.ma.masked_where(~np.eye(n, dtype=bool), np.ones((n, n)))
ax_heatmap.imshow(diagonal, cmap=ListedColormap([CONTROL_COLOR]),
                  interpolation="nearest", aspect="equal")
ax_heatmap.set_facecolor(UNTESTED_COLOR)
ax_heatmap.set_xlim(-0.5, n - 0.5)
ax_heatmap.set_ylim(n - 0.5, -0.5)
ax_heatmap.set_xticks(np.arange(n))
ax_heatmap.set_yticks(np.arange(n))
ax_heatmap.set_xticklabels(ordered_strains, rotation=90, fontsize=STRAIN_FONT)
ax_heatmap.set_yticklabels(ordered_strains, rotation=0, fontsize=STRAIN_FONT)
ax_heatmap.yaxis.tick_left()
ax_heatmap.tick_params(axis="both", length=2, pad=3)
ax_heatmap.set_xlabel(f"{SPECIES} strains", fontsize=12, labelpad=10)
ax_heatmap.set_ylabel(f"{SPECIES} strains", fontsize=12, labelpad=10)
for spine in ax_heatmap.spines.values():
    spine.set_visible(False)

# Optional group labels above the colour strip (blocks of >= 2 strains)
if SHOW_GROUP_LABELS:
    start = 0
    for i in range(1, n + 1):
        if i == n or ordered_groups[i] != ordered_groups[start]:
            if i - start >= 2:
                ax_col_colors.text((start + i - 1) / 2, -1.1, ordered_groups[start],
                                   ha="center", va="bottom", fontsize=GROUP_FONT,
                                   rotation=90, clip_on=False)
            start = i

# Legends
kd_legend = [Patch(facecolor=c, label=l) for c, l in zip(KD_COLORS, KD_LABELS)]
kd_legend.append(Patch(facecolor=CONTROL_COLOR, label="1 – Kin positive control"))
kd_legend.append(Patch(facecolor=UNTESTED_COLOR, label="Not tested"))
legend_kd = ax_legend.legend(handles=kd_legend, title="KD phenotype",
                             loc="upper left", bbox_to_anchor=(0, 1),
                             frameon=False, fontsize=10, title_fontsize=11)
ax_legend.add_artist(legend_kd)
ax_legend.legend(handles=[Patch(facecolor=group_colors[g], label=g) for g in groups],
                 title="Candidate KD group", loc="upper left",
                 bbox_to_anchor=(0, 0.76), frameon=False, fontsize=9,
                 title_fontsize=11, ncol=1, columnspacing=0.8, handlelength=1.2)

fig.suptitle(f"Kin discrimination phenotypes across ANI-clustered {SPECIES} strains",
             fontsize=17, y=0.985)

prefix = os.path.join(OUTPUT_DIR, "KD_heatmap_ANI_candidate_groups")
fig.savefig(prefix + ".png", dpi=400, bbox_inches="tight", facecolor="white")
fig.savefig(prefix + ".pdf", bbox_inches="tight", facecolor="white")
fig.savefig(prefix + ".svg", bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Saved: {prefix}.png / .pdf / .svg")

# ---------------------------------------------------------------------------
# 6. Pearson correlation (descriptive)
# ---------------------------------------------------------------------------
A = ani_panel.to_numpy(dtype=float)
idx = {s: i for i, s in enumerate(strains)}
ii = df["S1"].map(idx).to_numpy(dtype=int)
jj = df["S2"].map(idx).to_numpy(dtype=int)
kd = df["Phenotype"].to_numpy(dtype=float)
ani = A[ii, jj]
ani_s10 = df["ANI value"].to_numpy(dtype=float)

print("\nInput cross-check")
print(f"  Tested pairs: {len(df)} of {n * (n - 1) // 2} possible")
print(f"  Maximum S5 vs S10 ANI discrepancy: {np.max(np.abs(ani - ani_s10))}")
print(f"  Pearson r (ANI from S10): {pearsonr(ani_s10, kd).statistic:.4f}")
print(f"  Pearson r (ANI from S5):  {pearsonr(ani, kd).statistic:.4f}")

# ---------------------------------------------------------------------------
# 7. Strain-label permutation test (exploratory)
# ---------------------------------------------------------------------------
rng = np.random.default_rng(SEED)
r_obs = np.corrcoef(ani, kd)[0, 1]
r_perm = np.empty(N_PERM)
for t in range(N_PERM):
    q = rng.permutation(n)
    r_perm[t] = np.corrcoef(A[q[ii], q[jj]], kd)[0, 1]
p_two = (1 + np.count_nonzero(np.abs(r_perm) >= abs(r_obs))) / (N_PERM + 1)
p_one = (1 + np.count_nonzero(r_perm >= r_obs)) / (N_PERM + 1)

print(f"\nStrain-label permutation test ({N_PERM} permutations, exploratory)")
print(f"  Observed r = {r_obs:.6f}")
print(f"  Two-sided p = {p_two:.6f}; one-sided p = {p_one:.6f}")
print(f"  Null mean = {r_perm.mean():.6f}; null SD = {r_perm.std():.6f}")
pd.DataFrame({"permuted_r": r_perm}).to_csv(
    os.path.join(OUTPUT_DIR, "ANI_KD_permutation_9999_null.csv"), index=False)

# ---------------------------------------------------------------------------
# 8. Bayesian ordinal mixed-effects model (primary inference)
# ---------------------------------------------------------------------------
if RUN_ORDINAL_MODEL:
    import arviz as az
    import pymc as pm

    categories = {0.0: 0, 0.25: 1, 0.5: 2, 1.0: 3}
    y = np.array([categories[float(v)] for v in kd], dtype=int)
    # Centred for sampler stability; beta remains per 1 percentage point of ANI.
    x = ani - ani.mean()

    with pm.Model() as model:
        beta = pm.Normal("beta_ANI_per_1pct", mu=0, sigma=2)
        sigma_strain = pm.HalfNormal("sigma_strain", sigma=1)
        z_strain = pm.Normal("z_strain", mu=0, sigma=1, shape=n)
        u = pm.Deterministic("strain_effect", z_strain * sigma_strain)
        cutpoints = pm.Normal("cutpoints", mu=[-2, 0, 2], sigma=2, shape=3,
                              transform=pm.distributions.transforms.ordered,
                              initval=np.array([-2.0, 0.0, 2.0]))
        eta = beta * x + u[ii] + u[jj]
        pm.OrderedLogistic("KD", eta=eta, cutpoints=cutpoints, observed=y)
        trace = pm.sample(draws=1000, tune=1000, chains=2, cores=2,
                          random_seed=SEED, target_accept=0.95, progressbar=True)

    summary = az.summary(trace, var_names=["beta_ANI_per_1pct", "sigma_strain"],
                         hdi_prob=0.95)
    beta_samples = trace.posterior["beta_ANI_per_1pct"].values.flatten()
    n_divergences = int(trace.sample_stats["diverging"].sum().values)

    print("\nBayesian ordinal mixed-effects model")
    print(summary.to_string())
    print(f"  Posterior P(beta > 0): {float(np.mean(beta_samples > 0))}")
    print(f"  Posterior P(beta < 0): {float(np.mean(beta_samples < 0))}")
    print("\nDiagnostics")
    print(f"  Divergent transitions: {n_divergences}")
    print(f"  All R-hat < 1.01:      {bool((summary['r_hat'] < 1.01).all())}")
    print(f"  All bulk ESS > 400:    {bool((summary['ess_bulk'] > 400).all())}")

    # Posterior predictive check of marginal KD category frequencies
    with model:
        ppc = pm.sample_posterior_predictive(trace, var_names=["KD"], random_seed=SEED,
                                             progressbar=True, return_inferencedata=True)
    simulated = ppc.posterior_predictive["KD"].values
    simulated = simulated.reshape(-1, simulated.shape[-1])
    obs_prop = np.array([np.mean(y == k) for k in range(4)])
    sim_prop = np.array([np.mean(simulated == k, axis=1) for k in range(4)])
    ppc_table = pd.DataFrame({
        "KD score": ["0", "0.25", "0.5", "1"],
        "Observed proportion": obs_prop,
        "Predicted mean": sim_prop.mean(axis=1),
        "Predicted 2.5%": np.quantile(sim_prop, 0.025, axis=1),
        "Predicted 97.5%": np.quantile(sim_prop, 0.975, axis=1),
    })
    print("\nPosterior predictive check")
    print(ppc_table.to_string(index=False))

    ppc_table.to_csv(os.path.join(OUTPUT_DIR, "ANI_KD_posterior_predictive_frequencies.csv"),
                     index=False)
    summary.to_csv(os.path.join(OUTPUT_DIR, "ANI_KD_ordinal_model_summary.csv"))
    az.to_netcdf(trace, os.path.join(OUTPUT_DIR, "ANI_KD_ordinal_model_trace.nc"))

print(f"\nSaved to: {OUTPUT_DIR}")