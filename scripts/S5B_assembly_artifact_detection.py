#!/usr/bin/env python3
"""
S5B_assembly_artifact_detection.py: Assembly Artifact Detection

Author: Eva Stare

Identifies potential assembly artifacts among variant pairs that were *not*
merged by the conservative clustering step (Script S5A).  Assembly artifacts
arise when fragmented contigs split a single gene into two or more partial
open reading frames that Roary treats as separate clusters.  These fragments
typically show very high sequence identity (≥90–98 %) but limited alignment
coverage (20–60 %), because one copy represents only a truncated portion of
the full-length gene.

"""

import glob
import json
import logging
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
os.makedirs("results", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("results/assembly_artifact_detection.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)


# ===========================================================================
# Artifact detection
# ===========================================================================
class AssemblyArtifactDetector:
    """Detect gene-variant pairs that are likely assembly artifacts."""

    # Confidence-tier criteria
    CRITERIA = {
        "very_high_identity_fragment": {
            "min_identity": 98.0, "min_coverage": 20.0,
            "max_coverage": 55.0, "confidence": "VERY_HIGH",
            "description": "Very high identity fragments (almost certainly artifacts)",
        },
        "high_identity_fragment": {
            "min_identity": 95.0, "min_coverage": 25.0,
            "max_coverage": 60.0, "confidence": "HIGH",
            "description": "High identity fragments (likely assembly artifacts)",
        },
        "moderate_identity_fragment": {
            "min_identity": 90.0, "min_coverage": 30.0,
            "max_coverage": 50.0, "max_length_ratio": 0.6,
            "confidence": "MEDIUM",
            "description": "Moderate identity fragments with size bias",
        },
    }

    def __init__(self):
        self.config = {
            "conservative_results_dir": "results/conservative_clustering_output",
            "output_dir": "results/artifact_detection_output",
            "save_individual_results": True,
        }
        self.conservative_results: list[dict] = []

    # ---- data loading -----------------------------------------------------
    def load_conservative_results(self) -> bool:
        """Load per-gene JSON files from the conservative clustering step."""
        src = self.config["conservative_results_dir"]
        if not os.path.isdir(src):
            logger.error("Directory not found: %s — run Script S5A first.", src)
            return False
        files = glob.glob(str(Path(src) / "*_conservative_results.json"))
        if not files:
            logger.error("No result files in %s.", src)
            return False
        for fp in files:
            try:
                with open(fp) as fh:
                    self.conservative_results.append(json.load(fh))
            except Exception as exc:
                logger.warning("Could not load %s: %s", fp, exc)
        logger.info("Loaded %d gene results.", len(self.conservative_results))
        return True

    # ---- per-gene analysis ------------------------------------------------
    def analyze_gene(self, gene_result: dict) -> dict:
        """Scan un-clustered variant pairs for artifact signatures."""
        base_gene = gene_result["base_gene"]
        comparisons = gene_result.get("variant_comparisons", [])
        clusters = gene_result.get("conservative_clusters", [])

        # Already-merged pairs
        merged = {
            tuple(sorted(c["variants"]))
            for c in clusters if len(c["variants"]) == 2
        }

        candidates = []
        for comp in comparisons:
            pair = tuple(sorted([comp["variant1"], comp["variant2"]]))
            if pair in merged:
                continue
            result = self._classify(comp)
            if result["is_artifact"]:
                candidates.append(result)

        return {
            "base_gene": base_gene,
            "conservative_clusters": len(clusters),
            "additional_artifact_candidates": len(candidates),
            "artifact_candidates": candidates,
            "total_additional_features_saved": sum(
                c["features_saved"] for c in candidates
            ),
        }

    def _classify(self, comp: dict) -> dict:
        """Test a comparison against each artifact criterion."""
        identity = comp["identity_mean"]
        coverage = comp["coverage_mean"]

        # Length ratio from pairwise data
        pw = comp.get("pairwise_comparisons", [])
        ratios = []
        for p in pw:
            lo, hi = sorted([p["seq1_length"], p["seq2_length"]])
            if hi > 0:
                ratios.append(lo / hi)
        avg_ratio = float(np.mean(ratios)) if ratios else 1.0

        for name, crit in self.CRITERIA.items():
            if not (identity >= crit["min_identity"]
                    and coverage >= crit["min_coverage"]
                    and coverage <= crit["max_coverage"]):
                continue
            if "max_length_ratio" in crit and avg_ratio > crit["max_length_ratio"]:
                continue
            score = self._confidence_score(identity, coverage, avg_ratio, crit)
            return {
                "is_artifact": True,
                "artifact_type": name,
                "criterion_description": crit["description"],
                "confidence_level": crit["confidence"],
                "confidence_score": score,
                "recommendation": f"MERGE_{crit['confidence']}_CONFIDENCE",
                "variant1": comp["variant1"],
                "variant2": comp["variant2"],
                "features_saved": 1,
                "identity_mean": identity,
                "coverage_mean": coverage,
                "avg_length_ratio": avg_ratio,
            }

        return {
            "is_artifact": False,
            "variant1": comp["variant1"],
            "variant2": comp["variant2"],
            "identity_mean": identity,
            "coverage_mean": coverage,
            "avg_length_ratio": avg_ratio,
        }

    @staticmethod
    def _confidence_score(identity, coverage, ratio, crit) -> float:
        """Weighted score combining identity, coverage, and length ratio."""
        id_s = min(1.0, (identity - crit["min_identity"]) / (100 - crit["min_identity"]))
        cov_range = crit["max_coverage"] - crit["min_coverage"]
        cov_center = (crit["max_coverage"] + crit["min_coverage"]) / 2
        cov_s = 1.0 - abs(coverage - cov_center) / cov_range
        len_s = (1.0 - ratio / crit["max_length_ratio"]
                 if "max_length_ratio" in crit else 1.0)
        return min(1.0, 0.5 * id_s + 0.3 * cov_s + 0.2 * len_s)

    # ---- pipeline ---------------------------------------------------------
    def run(self) -> list[dict] | None:
        """Execute the full artifact detection pipeline."""
        os.makedirs(self.config["output_dir"], exist_ok=True)
        if not self.load_conservative_results():
            return None

        all_results = []
        for i, gr in enumerate(self.conservative_results, 1):
            logger.info("[%d/%d] %s", i, len(self.conservative_results), gr["base_gene"])
            ar = self.analyze_gene(gr)
            all_results.append(ar)
            if self.config["save_individual_results"]:
                out = Path(self.config["output_dir"]) / f"{gr['base_gene']}_artifact_analysis.json"
                with open(out, "w") as fh:
                    json.dump(ar, fh, indent=2, default=str)

        self._write_summary(all_results)
        logger.info("Artifact detection complete.")
        return all_results

    def _write_summary(self, results: list[dict]) -> None:
        """Write the summary Excel workbook."""
        gene_rows, artifact_rows = [], []
        for r in results:
            gene_rows.append({
                "Base_Gene": r["base_gene"],
                "Conservative_Clusters": r["conservative_clusters"],
                "Additional_Artifact_Candidates": r["additional_artifact_candidates"],
                "Additional_Features_Saved": r["total_additional_features_saved"],
            })
            for c in r["artifact_candidates"]:
                artifact_rows.append({
                    "Base_Gene": r["base_gene"],
                    "Variant1": c["variant1"],
                    "Variant2": c["variant2"],
                    "Artifact_Type": c["artifact_type"],
                    "Confidence_Level": c["confidence_level"],
                    "Confidence_Score": c["confidence_score"],
                    "Identity_Mean": c["identity_mean"],
                    "Coverage_Mean": c["coverage_mean"],
                    "Avg_Length_Ratio": c["avg_length_ratio"],
                })
        out = Path(self.config["output_dir"]) / "02_assembly_artifact_detection_summary.xlsx"
        with pd.ExcelWriter(out, engine="openpyxl") as w:
            pd.DataFrame(gene_rows).to_excel(w, sheet_name="Gene_Summary", index=False)
            if artifact_rows:
                df = pd.DataFrame(artifact_rows)
                df.to_excel(w, sheet_name="Artifact_Candidates", index=False)
        logger.info("Summary saved to %s", out)


# ===========================================================================
def main():
    detector = AssemblyArtifactDetector()
    detector.run()


if __name__ == "__main__":
    main()
