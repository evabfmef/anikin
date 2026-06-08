#!/usr/bin/env python3
"""
S5A_conservative_gene_clustering.py: Conservative Gene Variant Clustering

Author: Eva Stare

Identifies and clusters redundant gene variants in a pangenome presence–absence
matrix using pairwise amino acid sequence similarity. Variant pairs meeting
conservative thresholds (≥80% identity, ≥50% coverage) are merged.

"""

import json
import glob
import logging
import multiprocessing as mp
import os
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
os.makedirs("results", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("results/conservative_gene_clustering.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)


# ===========================================================================
# Main class
# ===========================================================================
class ConservativeGeneVariantClustering:
    """Cluster Roary gene variants by amino acid sequence similarity.

    Thresholds
    ----------
    min_aa_identity : float
        Minimum mean pairwise amino acid identity (default 80 %).
    min_coverage : float
        Minimum mean alignment coverage relative to the shorter sequence
        (default 50 %).
    """

    def __init__(self):
        self.config = self._setup_config()
        self.aligner = self._setup_local_aligner()
        self.mega_matrix: pd.DataFrame | None = None
        self.locus_mapping: pd.DataFrame | None = None
        self.strain_mapping: pd.DataFrame | None = None
        self.gene_variants: dict[str, list[str]] = {}
        self.strain_to_file: dict[str, str] = {}
        self.gene_to_locus: dict[str, str] = {}
        self.gbk_file_cache: dict[str, str] = {}

    # ---- configuration ----------------------------------------------------
    @staticmethod
    def _setup_config() -> dict:
        """Return default configuration dictionary."""
        return {
            # Input paths (adjust to local directory structure)
            "mega_matrix_file": "data/mega_matrix.xlsx",
            "locus_tag_file": "data/locus_tag_mapping.xlsx",
            "strain_list_file": "data/strain_list.xlsx",
            "gbk_folder_path": "data/gbk_files",
            # Output paths
            "output_dir": "results/conservative_clustering_output",
            "sequences_dir": "results/sequences_conservative",
            # Clustering thresholds
            "min_aa_identity": 80.0,
            "min_coverage": 50.0,
            # Sampling settings
            "max_sequences_per_variant": 10,
            "min_sequences_for_sampling": 10,
            # Performance
            "max_workers": min(32, mp.cpu_count()),
            "save_detailed_comparisons": True,
            "save_sequences_to_disk": True,
            "gbk_extensions": [".gbk", ".gb", ".gbff"],
        }

    def _setup_local_aligner(self) -> PairwiseAligner:
        """Configure a Smith–Waterman local aligner for protein sequences."""
        aligner = PairwiseAligner()
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5
        aligner.mode = "local"
        return aligner

    # ---- data loading -----------------------------------------------------
    def load_data(self) -> None:
        """Load the mega-matrix, locus-tag mapping, and strain list."""
        logger.info("Loading data files …")
        self.mega_matrix = pd.read_excel(
            self.config["mega_matrix_file"], sheet_name=0
        )
        self.locus_mapping = pd.read_excel(
            self.config["locus_tag_file"], sheet_name=0
        )
        self.strain_mapping = pd.read_excel(
            self.config["strain_list_file"], sheet_name=0
        )
        self._process_data()

    def _process_data(self) -> None:
        """Build internal look-up tables and identify multi-variant genes."""
        self.strain_to_file = dict(
            zip(self.strain_mapping["StrainName"],
                self.strain_mapping["FileName"])
        )
        self.gene_to_locus = {}
        for _, row in self.locus_mapping.iterrows():
            names = row.get("All Gene Names & Synonyms")
            tag = row.get("Locus Tag")
            if pd.notna(names) and pd.notna(tag):
                for name in str(names).split(", "):
                    name = name.strip()
                    if name:
                        self.gene_to_locus[name] = tag
        self._identify_gene_variants()
        self._index_gbk_files()

    def _identify_gene_variants(self) -> None:
        """Parse variant suffixes (e.g. .0, .1) and retain multi-variant genes."""
        variant_pattern = re.compile(r"^(.+)\.(\d+)$")
        for gene_name in self.mega_matrix["Feature_Name"].dropna():
            m = variant_pattern.match(str(gene_name))
            if m:
                base = m.group(1)
                self.gene_variants.setdefault(base, []).append(gene_name)
        self.gene_variants = {
            k: v for k, v in self.gene_variants.items() if len(v) > 1
        }
        logger.info(
            "Found %d genes with multiple variants", len(self.gene_variants)
        )

    def _index_gbk_files(self) -> None:
        """Build a filename → path cache for GenBank files."""
        self.gbk_file_cache = {}
        for ext in self.config["gbk_extensions"]:
            for path in glob.glob(
                str(Path(self.config["gbk_folder_path"]) / f"*{ext}")
            ):
                base = os.path.basename(path)
                stem = os.path.splitext(base)[0]
                self.gbk_file_cache[stem] = path
                self.gbk_file_cache[base] = path
                for suffix in ("_edit", "_edited", "_edit2"):
                    if stem.endswith(suffix):
                        self.gbk_file_cache[stem[: -len(suffix)]] = path

    # ---- sequence extraction ----------------------------------------------
    def _find_gbk_file(self, strain_name: str) -> str | None:
        filename = self.strain_to_file.get(strain_name, strain_name)
        for candidate in (
            filename, f"{filename}.gbk", f"{filename}.gb",
            f"{filename}_edit2.gbk", f"{filename}_edit.gbk",
            f"{filename}_edited.gbk",
        ):
            if candidate in self.gbk_file_cache:
                return self.gbk_file_cache[candidate]
        return None

    def _extract_by_locus_tag(self, gbk_path: str, locus_tag: str) -> dict | None:
        """Search a GenBank file for a CDS matching *locus_tag*."""
        try:
            for record in SeqIO.parse(gbk_path, "genbank"):
                for feat in record.features:
                    if feat.type != "CDS":
                        continue
                    # Primary: locus tag embedded in product annotation
                    product = feat.qualifiers.get("product", [""])[0]
                    if f"[locus_tag={locus_tag}]" in product:
                        return self._seq_from_feature(feat, record)
                    if locus_tag.startswith("BSU_"):
                        alt = locus_tag.replace("BSU_", "")
                        if (f"[locus_tag={alt}]" in product
                                or f"[locus_tag=BSU{alt}]" in product):
                            return self._seq_from_feature(feat, record)
                    # Fallback: standard locus_tag qualifier
                    ft_tag = feat.qualifiers.get("locus_tag", [""])[0]
                    if ft_tag == locus_tag:
                        return self._seq_from_feature(feat, record)
                    if locus_tag.startswith("BSU_"):
                        custom = locus_tag.replace("BSU_", "CUSTOM_")
                        if ft_tag == custom:
                            return self._seq_from_feature(feat, record)
                    if locus_tag.replace("_", "") == ft_tag.replace("_", ""):
                        return self._seq_from_feature(feat, record)
        except Exception as exc:
            logger.warning("Error reading %s: %s", os.path.basename(gbk_path), exc)
        return None

    @staticmethod
    def _seq_from_feature(feature, record) -> dict:
        """Extract the nucleotide sequence for a CDS feature."""
        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = feature.location.strand
        seq = record.seq[start:end]
        if strand == -1:
            seq = seq.reverse_complement()
        return {
            "sequence": str(seq), "start": start,
            "end": end, "strand": strand, "length": len(seq),
        }

    @staticmethod
    def _translate(dna_sequence: str) -> str:
        """Translate a CDS in the longest reading frame (forward + reverse)."""
        seq = Seq(dna_sequence)
        best = ""
        for frame in range(3):
            for s in (seq, seq.reverse_complement()):
                try:
                    t = str(s[frame:].translate(to_stop=False))
                    if len(t) > len(best):
                        best = t
                except Exception:
                    pass
        return best

    def _get_strains_with_variant(self, variant_name: str) -> list[str]:
        """Return strain names where *variant_name* is present (value == 1)."""
        row = self.mega_matrix.loc[
            self.mega_matrix["Feature_Name"] == variant_name
        ]
        if row.empty:
            return []
        row = row.iloc[0]
        strain_cols = [c for c in self.mega_matrix.columns if c != "Feature_Name"]
        return [s for s in strain_cols if pd.notna(row.get(s)) and row[s] == 1]

    def extract_sequences(self, base_gene: str, variants: list[str]) -> dict:
        """Extract amino acid sequences for all variants of a gene."""
        # Resolve locus tag
        locus_tag = None
        for name in [base_gene] + (base_gene.split(", ") if ", " in base_gene else []):
            if name in self.gene_to_locus:
                locus_tag = self.gene_to_locus[name]
                break
        if locus_tag is None:
            return {}

        result: dict[str, list[dict]] = {}
        for variant in variants:
            strains = self._get_strains_with_variant(variant)
            if not strains:
                continue
            # Sampling strategy: use all strains up to threshold, then subsample
            target = (
                len(strains)
                if len(strains) <= self.config["min_sequences_for_sampling"]
                else self.config["max_sequences_per_variant"]
            )
            step = max(1, len(strains) // target)
            ordered = [strains[i * step] for i in range(min(target, len(strains)))]
            backup = [s for s in strains if s not in ordered]
            seqs: list[dict] = []
            for strain in ordered + backup:
                if len(seqs) >= target:
                    break
                gbk = self._find_gbk_file(strain)
                if gbk is None:
                    continue
                info = self._extract_by_locus_tag(gbk, locus_tag)
                if info is None:
                    continue
                aa = self._translate(info["sequence"])
                seqs.append({
                    "strain": strain, "variant": variant,
                    "sequence": info["sequence"], "aa_sequence": aa,
                    "aa_length": len(aa), **info,
                })
            if seqs:
                result[variant] = seqs
        return result

    # ---- alignment --------------------------------------------------------
    def _align(self, seq1: str, seq2: str) -> dict | None:
        """Compute identity and coverage from a local alignment."""
        s1 = seq1.replace("*", "").strip()
        s2 = seq2.replace("*", "").strip()
        if not s1 or not s2:
            return None
        if s1 == s2:
            return {
                "identity": 100.0, "coverage": 100.0,
                "alignment_length": len(s1),
                "seq1_length": len(s1), "seq2_length": len(s2),
                "alignment_type": "identical",
            }
        # Fast path for near-identical equal-length sequences
        if len(s1) == len(s2):
            matches = sum(a == b for a, b in zip(s1, s2))
            if matches / len(s1) > 0.95:
                return {
                    "identity": matches / len(s1) * 100,
                    "coverage": 100.0,
                    "alignment_length": len(s1),
                    "seq1_length": len(s1), "seq2_length": len(s2),
                    "alignment_type": "high_similarity_direct",
                }
        # Full local alignment with overflow fallback
        try:
            aln = self.aligner.align(s1, s2)[0]
            a1, a2 = str(aln[0]), str(aln[1])
            matches = aligned = 0
            for a, b in zip(a1, a2):
                if a != "-" and b != "-":
                    aligned += 1
                    matches += a == b
            if aligned == 0:
                return None
            shorter = min(len(s1), len(s2))
            return {
                "identity": matches / aligned * 100,
                "coverage": aligned / shorter * 100,
                "alignment_length": aligned,
                "seq1_length": len(s1), "seq2_length": len(s2),
                "alignment_type": "local",
            }
        except (OverflowError, ValueError):
            return self._fallback_align(s1, s2)

    def _fallback_align(self, seq1: str, seq2: str) -> dict | None:
        """Sliding-window identity for sequences that cause alignment overflow."""
        short, long = (seq1, seq2) if len(seq1) <= len(seq2) else (seq2, seq1)
        if not short:
            return None
        best_id = 0.0
        w = len(short)
        for i in range(len(long) - w + 1):
            m = sum(short[j] == long[i + j] for j in range(w))
            best_id = max(best_id, m / w * 100)
        return {
            "identity": best_id,
            "coverage": len(short) / len(long) * 100,
            "alignment_length": len(short),
            "seq1_length": len(seq1), "seq2_length": len(seq2),
            "alignment_type": "fallback_sliding_window",
        }

    # ---- comparison and clustering ----------------------------------------
    def compare_variants(self, variant_seqs: dict) -> list[dict]:
        """All-versus-all pairwise comparison between variant groups."""
        variants = list(variant_seqs)
        results = []
        for i, v1 in enumerate(variants):
            for v2 in variants[i + 1:]:
                pairwise = []
                for s1 in variant_seqs[v1]:
                    for s2 in variant_seqs[v2]:
                        m = self._align(s1["aa_sequence"], s2["aa_sequence"])
                        if m:
                            pairwise.append({
                                "variant1": v1, "variant2": v2,
                                "strain1": s1["strain"], "strain2": s2["strain"],
                                **m,
                            })
                if pairwise:
                    ids = [p["identity"] for p in pairwise]
                    covs = [p["coverage"] for p in pairwise]
                    results.append({
                        "variant1": v1, "variant2": v2,
                        "num_comparisons": len(pairwise),
                        "identity_mean": np.mean(ids),
                        "identity_std": np.std(ids),
                        "identity_min": np.min(ids),
                        "identity_max": np.max(ids),
                        "coverage_mean": np.mean(covs),
                        "coverage_std": np.std(covs),
                        "coverage_min": np.min(covs),
                        "coverage_max": np.max(covs),
                        "pairwise_comparisons": pairwise,
                    })
        return results

    def cluster(self, comparisons: list[dict]) -> list[dict]:
        """Apply identity and coverage thresholds to form merge clusters."""
        clusters = []
        for comp in comparisons:
            if (comp["identity_mean"] >= self.config["min_aa_identity"]
                    and comp["coverage_mean"] >= self.config["min_coverage"]):
                clusters.append({
                    "variants": [comp["variant1"], comp["variant2"]],
                    "cluster_type": "CONSERVATIVE",
                    "recommendation": "MERGE_CONSERVATIVE",
                    "confidence_level": "MEDIUM",
                    "features_saved": 1,
                    "identity_mean": comp["identity_mean"],
                    "identity_std": comp["identity_std"],
                    "identity_min": comp["identity_min"],
                    "identity_max": comp["identity_max"],
                    "coverage_mean": comp["coverage_mean"],
                    "coverage_std": comp["coverage_std"],
                    "coverage_min": comp["coverage_min"],
                    "coverage_max": comp["coverage_max"],
                    "num_comparisons": comp["num_comparisons"],
                })
        return clusters

    # ---- gene-level processing --------------------------------------------
    def process_gene(self, base_gene: str, variants: list[str]) -> dict | None:
        """Extract, compare, and cluster variants for a single gene."""
        try:
            seqs = self.extract_sequences(base_gene, variants)
            if len(seqs) < 2:
                return None
            comparisons = self.compare_variants(seqs)
            if not comparisons:
                return None
            clusters = self.cluster(comparisons)
            result = {
                "base_gene": base_gene,
                "variants": variants,
                "variant_count": len(variants),
                "sequences_extracted": {v: len(s) for v, s in seqs.items()},
                "total_sequences": sum(len(s) for s in seqs.values()),
                "variant_comparisons": comparisons,
                "conservative_clusters": clusters,
                "features_saved": sum(c["features_saved"] for c in clusters),
                "clustering_summary": {
                    "total_variant_pairs": len(comparisons),
                    "clustered_pairs": len(clusters),
                    "rejection_rate": (
                        (len(comparisons) - len(clusters)) / len(comparisons) * 100
                        if comparisons else 0
                    ),
                },
            }
            if self.config["save_sequences_to_disk"]:
                out = Path(self.config["output_dir"]) / f"{base_gene}_conservative_results.json"
                with open(out, "w") as fh:
                    json.dump(result, fh, indent=2, default=str)
            return result
        except Exception as exc:
            logger.error("Error processing %s: %s", base_gene, exc)
            return None

    # ---- main entry point -------------------------------------------------
    def run(self) -> list[dict]:
        """Execute the full conservative clustering pipeline."""
        os.makedirs(self.config["output_dir"], exist_ok=True)
        os.makedirs(self.config["sequences_dir"], exist_ok=True)
        self.load_data()

        results = []
        total = len(self.gene_variants)
        for i, (gene, variants) in enumerate(self.gene_variants.items(), 1):
            logger.info("[%d/%d] %s", i, total, gene)
            r = self.process_gene(gene, variants)
            if r:
                results.append(r)

        self._write_summary(results)
        logger.info("Conservative clustering complete — %d/%d genes processed.", len(results), total)
        return results

    def _write_summary(self, results: list[dict]) -> None:
        """Write an Excel summary of all clustering results."""
        rows, cluster_rows = [], []
        for r in results:
            rows.append({
                "Base_Gene": r["base_gene"],
                "Variant_Count": r["variant_count"],
                "Total_Sequences": r["total_sequences"],
                "Variant_Pairs_Compared": r["clustering_summary"]["total_variant_pairs"],
                "Conservative_Clusters": len(r["conservative_clusters"]),
                "Features_Saved": r["features_saved"],
                "Rejection_Rate": r["clustering_summary"]["rejection_rate"],
            })
            for c in r["conservative_clusters"]:
                cluster_rows.append({
                    "Base_Gene": r["base_gene"],
                    "Variants": ", ".join(c["variants"]),
                    "Cluster_Type": c["cluster_type"],
                    "Confidence_Level": c["confidence_level"],
                    "Features_Saved": c["features_saved"],
                    "Identity_Mean": c["identity_mean"],
                    "Identity_Min": c["identity_min"],
                    "Coverage_Mean": c["coverage_mean"],
                    "Coverage_Min": c["coverage_min"],
                    "Num_Comparisons": c["num_comparisons"],
                })
        out = Path(self.config["output_dir"]) / "01_conservative_clustering_summary.xlsx"
        with pd.ExcelWriter(out, engine="openpyxl") as w:
            pd.DataFrame(rows).to_excel(w, sheet_name="Gene_Summary", index=False)
            if cluster_rows:
                pd.DataFrame(cluster_rows).to_excel(w, sheet_name="Conservative_Clusters", index=False)
        logger.info("Summary saved to %s", out)


# ===========================================================================
# Entry point
# ===========================================================================
def main():
    clustering = ConservativeGeneVariantClustering()
    clustering.run()


if __name__ == "__main__":
    main()
