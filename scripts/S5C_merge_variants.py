#!/usr/bin/env python3
"""
S5C_merge_variants.py: Merge Gene Variants Into a Reduced Presence–Absence Matrix

Author: Eva Stare

Integrates the merge decisions from the conservative clustering step
(Script S5A) and the assembly-artifact detection step (Script S5B) and applies
them to the original Roary presence–absence matrix.

"""

import logging
import os
from collections import defaultdict
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
os.makedirs("results", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("results/variant_merging.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)


# ===========================================================================
class VariantMerger:
    """Merge gene variants from conservative clustering and artifact detection."""

    def __init__(self):
        self.config = {
            # Input
            "mega_matrix_file": "data/mega_matrix.xlsx",
            "clustering_summary_file": (
                "results/conservative_clustering_output/"
                "01_conservative_clustering_summary.xlsx"
            ),
            "artifact_summary_file": (
                "results/artifact_detection_output/"
                "02_assembly_artifact_detection_summary.xlsx"
            ),
            # Output
            "output_dir": "results/merged_matrix_output",
            "merged_matrix_file": "results/merged_matrix_output/mega_matrix_merged.xlsx",
            "merge_mapping_file": "results/merged_matrix_output/merge_mapping.xlsx",
            "merge_report_file": "results/merged_matrix_output/merge_report.xlsx",
            # Merge options
            "rejection_rate_threshold": 0.0,
            "merge_artifacts": True,
            "artifact_confidence_levels": ["HIGH", "VERY_HIGH", "MEDIUM"],
        }
        self.mega_matrix: pd.DataFrame | None = None
        self.gene_summary: pd.DataFrame | None = None
        self.clusters_df: pd.DataFrame | None = None
        self.artifacts_df: pd.DataFrame | None = None
        self.merge_groups: dict[str, list[str]] = {}
        self.merge_mapping: dict[str, str] = {}
        self.merge_sources: dict[str, set[str]] = {}

    # ---- data loading -----------------------------------------------------
    def load_data(self) -> None:
        self.mega_matrix = pd.read_excel(self.config["mega_matrix_file"], sheet_name=0)
        logger.info("Matrix: %d features × %d columns", *self.mega_matrix.shape)

        self.gene_summary = pd.read_excel(
            self.config["clustering_summary_file"], sheet_name="Gene_Summary"
        )
        self.clusters_df = pd.read_excel(
            self.config["clustering_summary_file"], sheet_name="Conservative_Clusters"
        )
        logger.info("Conservative clusters: %d", len(self.clusters_df))

        if self.config["merge_artifacts"]:
            path = self.config["artifact_summary_file"]
            if os.path.isfile(path):
                df = pd.read_excel(path, sheet_name="Artifact_Candidates")
                levels = self.config["artifact_confidence_levels"]
                self.artifacts_df = df[df["Confidence_Level"].isin(levels)]
                logger.info("Artifact candidates: %d", len(self.artifacts_df))
            else:
                logger.warning("Artifact file not found — continuing without.")
                self.config["merge_artifacts"] = False

    # ---- graph construction -----------------------------------------------
    def build_merge_groups(self) -> None:
        """Build merge groups via connected-component analysis."""
        adjacency: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
        all_variants: dict[str, set[str]] = defaultdict(set)
        edge_src: dict[tuple, str] = {}

        # Conservative edges
        eligible = self.gene_summary[
            self.gene_summary["Rejection_Rate"] <= self.config["rejection_rate_threshold"]
        ]
        for base_gene in eligible["Base_Gene"]:
            for _, row in self.clusters_df[
                self.clusters_df["Base_Gene"] == base_gene
            ].iterrows():
                variants = self._parse_variants(row["Variants"], base_gene)
                for i, v1 in enumerate(variants):
                    for v2 in variants[i + 1:]:
                        adjacency[base_gene][v1].add(v2)
                        adjacency[base_gene][v2].add(v1)
                        all_variants[base_gene].update([v1, v2])
                        edge_src.setdefault(tuple(sorted([v1, v2])), "CONSERVATIVE")

        # Artifact edges
        if self.config["merge_artifacts"] and self.artifacts_df is not None:
            for _, row in self.artifacts_df.iterrows():
                bg, v1, v2 = row["Base_Gene"], row["Variant1"], row["Variant2"]
                adjacency[bg][v1].add(v2)
                adjacency[bg][v2].add(v1)
                all_variants[bg].update([v1, v2])
                edge_src.setdefault(tuple(sorted([v1, v2])), "ARTIFACT")

        # Connected components (DFS)
        for bg, adj in adjacency.items():
            visited: set[str] = set()
            for node in all_variants[bg]:
                if node in visited:
                    continue
                comp: set[str] = set()
                stack = [node]
                while stack:
                    n = stack.pop()
                    if n in visited:
                        continue
                    visited.add(n)
                    comp.add(n)
                    stack.extend(adj[n] - visited)
                if len(comp) > 1:
                    rep = sorted(comp)[0]
                    self.merge_groups[rep] = sorted(comp)
                    for v in comp:
                        self.merge_mapping[v] = rep
                    sources = set()
                    for v1 in comp:
                        for v2 in comp:
                            if v1 < v2:
                                ek = (v1, v2)
                                if ek in edge_src:
                                    sources.add(edge_src[ek])
                    self.merge_sources[rep] = sources

        logger.info("Merge groups: %d", len(self.merge_groups))

    @staticmethod
    def _parse_variants(variants_str: str, base_gene: str) -> list[str]:
        """Parse comma-separated variant string, handling multi-name genes."""
        if "," in base_gene:
            n_parts = base_gene.count(",") + 1
            parts = [p.strip() for p in variants_str.split(",")]
            return [", ".join(parts[i:i + n_parts])
                    for i in range(0, len(parts), n_parts)
                    if i + n_parts <= len(parts)]
        return [v.strip() for v in variants_str.split(",")]

    # ---- applying merges --------------------------------------------------
    def apply_merges(self) -> tuple[pd.DataFrame, list[dict], list[dict]]:
        """Merge variant rows in the mega-matrix."""
        mat = self.mega_matrix.copy()
        strain_cols = [c for c in mat.columns if c != "Feature_Name"]
        dropped, details = [], []

        for rep, members in self.merge_groups.items():
            if len(members) <= 1:
                continue
            idxs = {}
            for v in members:
                rows = mat[mat["Feature_Name"] == v]
                if rows.empty and "." in v:
                    short_base, suf = v.rsplit(".", 1)
                    esc = short_base.replace(" ", r"\s*")
                    rows = mat[mat["Feature_Name"].str.contains(
                        f"^{esc},.*\\.{suf}$", regex=True, na=False
                    )]
                if not rows.empty:
                    idxs[v] = rows.index[0]
            if not idxs:
                continue
            rep_idx = idxs.get(rep, next(iter(idxs.values())))
            for col in strain_cols:
                mat.loc[rep_idx, col] = max(
                    mat.loc[idx, col] for idx in idxs.values()
                )
            src = "+".join(sorted(self.merge_sources.get(rep, {"UNKNOWN"})))
            for v, idx in idxs.items():
                if v != rep:
                    dropped.append({
                        "Dropped_Variant": v, "Representative": rep,
                        "Merge_Source": src,
                    })
                    mat = mat.drop(idx)
            details.append({
                "Representative": rep,
                "Variants_Merged": ", ".join(members),
                "Num_Variants": len(members),
                "Features_Saved": len(members) - 1,
                "Merge_Source": src,
            })

        mat = mat.reset_index(drop=True)
        logger.info("Original: %d → Merged: %d (removed %d)",
                     len(self.mega_matrix), len(mat), len(dropped))
        return mat, dropped, details

    # ---- output -----------------------------------------------------------
    def save(self, merged: pd.DataFrame, dropped: list, details: list) -> None:
        os.makedirs(self.config["output_dir"], exist_ok=True)
        merged.to_excel(self.config["merged_matrix_file"], index=False)
        pd.DataFrame([
            {"Variant": v, "Representative": r,
             "Merge_Source": "+".join(sorted(self.merge_sources.get(r, set())))}
            for v, r in self.merge_mapping.items()
        ]).to_excel(self.config["merge_mapping_file"], index=False)

        with pd.ExcelWriter(self.config["merge_report_file"], engine="openpyxl") as w:
            summary = {
                "Metric": ["Original_Features", "Merged_Features",
                           "Features_Removed", "Reduction_Pct",
                           "Total_Merge_Groups"],
                "Value": [len(self.mega_matrix), len(merged), len(dropped),
                          f"{len(dropped)/len(self.mega_matrix)*100:.2f}%",
                          len(self.merge_groups)],
            }
            pd.DataFrame(summary).to_excel(w, sheet_name="Summary", index=False)
            if details:
                pd.DataFrame(details).to_excel(w, sheet_name="Merge_Details", index=False)
            if dropped:
                pd.DataFrame(dropped).to_excel(w, sheet_name="Dropped_Variants", index=False)
        logger.info("Results saved to %s", self.config["output_dir"])

    # ---- main pipeline ----------------------------------------------------
    def run(self) -> tuple[pd.DataFrame, list, list]:
        self.load_data()
        self.build_merge_groups()
        merged, dropped, details = self.apply_merges()
        self.save(merged, dropped, details)
        return merged, dropped, details


# ===========================================================================
def main():
    merger = VariantMerger()
    merged, dropped, details = merger.run()
    pct = len(dropped) / len(merger.mega_matrix) * 100
    logger.info("Final matrix: %d features (%.1f %% reduction)", len(merged), pct)


if __name__ == "__main__":
    main()
