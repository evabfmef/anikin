#!/usr/bin/env python3
"""
S7_phaster_to_votu_mapping.py

Author: Eva Stare

Cross-reference PHASTER prophage predictions against the Bacillus subtilis
prophage vOTU classification of Stefanic et al. (2025).

Each PHASTER-predicted prophage region is assigned to the viral operational
taxonomic unit (vOTU) yielding the highest aggregate BLASTn hit score,
subject to user-specified coverage and identity thresholds. Only PHASTER
regions classified as 'intact' contribute to the output tables;
'questionable' and 'incomplete' predictions are discarded after BLAST
search and before downstream summarisation.

"""

import argparse
import csv
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# PHASTER output parsing
# ---------------------------------------------------------------------------

_COORD_RE = re.compile(r"(\d{4,})-(\d{4,})")
_HIT_RE = re.compile(r"(PHAGE_[A-Za-z0-9_.]+)\(\d+\)")
_TRAIL_COORD_RE = re.compile(r"(\d+)-(\d+)\s*$")
_LEAD_INT_RE = re.compile(r">\s*(\d+)\s")


def parse_phaster_summary(path):
    """Parse a PHASTER summary.txt file into a list of region records.

    Handles both single-contig PHASTER outputs (region position recorded as
    a bare 'START-END' token) and multi-contig outputs (region position
    embedded within a comma-separated contig description that contains
    embedded whitespace and additional metadata).
    """
    regions = []
    with open(path) as fh:
        for line in fh:
            tokens = line.strip().split()
            if not tokens or not tokens[0].isdigit():
                continue
            try:
                region_n = int(tokens[0])
                length_kb = float(tokens[1].replace("Kb", ""))
                completeness_full = tokens[2]
                completeness = completeness_full.split("(")[0]
                score = int(completeness_full.split("(")[1].rstrip(")"))
                coord_match = _COORD_RE.search(line)
                if not coord_match:
                    continue
                start, end = int(coord_match.group(1)), int(coord_match.group(2))
                hit_match = _HIT_RE.search(line)
                top_hit = hit_match.group(1) if hit_match else ""
            except (ValueError, IndexError):
                continue
            regions.append({
                "region": region_n,
                "length_kb": length_kb,
                "completeness": completeness,
                "score": score,
                "start": start,
                "end": end,
                "length_bp": end - start + 1,
                "top_phaster_hit": top_hit,
            })
    return regions


def collect_phaster_metadata(phaster_dir):
    """Collect PHASTER region metadata across all strains in a directory."""
    metadata = {}
    for summary in sorted(Path(phaster_dir).glob("*_summary.txt")):
        strain = summary.name.replace("_summary.txt", "")
        for record in parse_phaster_summary(summary):
            metadata[(strain, record["region"])] = record
    return metadata


def assemble_query_fasta(phaster_dir, out_path):
    """Concatenate all per-strain prophage FASTA files into one query FASTA.

    Output headers follow the format '>{STRAIN}__region{N}__{START}-{END}'
    so that downstream BLAST results can be unambiguously attributed.
    Returns a dictionary mapping each generated header to a metadata tuple
    (strain, region_n, start, end, length_bp).
    """
    regions_meta = {}

    def _flush(out_fh, header, chunks, meta):
        if header is None:
            return
        seq = "".join(chunks)
        out_fh.write(f">{header}\n{seq}\n")
        regions_meta[header] = (*meta, len(seq))

    with open(out_path, "w") as out_fh:
        for fna in sorted(Path(phaster_dir).glob("*_phage_regions.fna")):
            strain = fna.name.replace("_phage_regions.fna", "")
            current_header = None
            current_meta = None
            chunks = []
            region_counter = 0
            with open(fna) as fh:
                for line in fh:
                    line = line.rstrip()
                    if line.startswith(">"):
                        _flush(out_fh, current_header, chunks, current_meta)
                        region_counter += 1
                        coord_match = _TRAIL_COORD_RE.search(line[1:])
                        if not coord_match:
                            raise ValueError(
                                f"No coordinates in FASTA header: {line!r}"
                            )
                        start = int(coord_match.group(1))
                        end = int(coord_match.group(2))
                        lead_match = _LEAD_INT_RE.match(line)
                        region_n = (
                            int(lead_match.group(1)) if lead_match else region_counter
                        )
                        current_header = f"{strain}__region{region_n}__{start}-{end}"
                        current_meta = (strain, region_n, start, end)
                        chunks = []
                    else:
                        chunks.append(line)
                _flush(out_fh, current_header, chunks, current_meta)
    return regions_meta


# ---------------------------------------------------------------------------
# vOTU reference assembly
# ---------------------------------------------------------------------------

def assemble_target_fasta(votu_dir, out_path):
    """Concatenate all vOTU reference FASTAs into a single target FASTA.

    Each sequence header is prefixed with the originating vOTU identifier
    (derived from the source filename stem) so the vOTU label can be
    recovered from BLAST output.
    """
    votu_counts = defaultdict(int)
    with open(out_path, "w") as out_fh:
        for pattern in ("*.fasta", "*.fa", "*.fna"):
            for fa in sorted(Path(votu_dir).glob(pattern)):
                stem = fa.stem
                votu_id = stem if stem.lower().startswith("votu") else f"vOTU_{stem}"
                with open(fa) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            original = line[1:].rstrip()
                            out_fh.write(f">{votu_id}__{original}\n")
                            votu_counts[votu_id] += 1
                        else:
                            out_fh.write(line)
    return dict(votu_counts)


# ---------------------------------------------------------------------------
# BLAST invocation
# ---------------------------------------------------------------------------

def run_makeblastdb(target_fasta, db_prefix):
    """Build a nucleotide BLAST database from the merged vOTU FASTA."""
    subprocess.run(
        ["makeblastdb", "-in", str(target_fasta),
         "-dbtype", "nucl", "-out", str(db_prefix)],
        check=True, capture_output=True,
    )


def run_blastn(query_fasta, db_prefix, out_tsv, task, threads):
    """Run BLASTn (default task: dc-megablast) and write tabular output."""
    subprocess.run(
        ["blastn",
         "-task", task,
         "-query", str(query_fasta),
         "-db", str(db_prefix),
         "-out", str(out_tsv),
         "-outfmt",
         "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
         "-evalue", "1e-20",
         "-num_threads", str(threads),
         "-max_target_seqs", "500"],
        check=True,
    )


# ---------------------------------------------------------------------------
# vOTU assignment
# ---------------------------------------------------------------------------

def assign_votus(blast_tsv, regions_meta, min_coverage, min_identity):
    """Assign each query region to its best-matching vOTU.

    For each (query, vOTU) pair, BLASTn HSPs are summed by alignment length
    and length-weighted percent identity. The aggregate score for a vOTU is
    aligned_bp * mean_identity / 100. The vOTU with the highest aggregate
    score is selected, subject to minimum coverage and identity thresholds;
    queries failing either threshold are marked 'unassigned'.
    """
    aggregates = defaultdict(
        lambda: defaultdict(lambda: {"aln": 0, "id_sum": 0.0})
    )
    with open(blast_tsv) as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            qid, sid = parts[0], parts[1]
            pident = float(parts[2])
            length = int(parts[3])
            votu = sid.split("__", 1)[0]
            aggregates[qid][votu]["aln"] += length
            aggregates[qid][votu]["id_sum"] += pident * length

    assignments = {}
    for qid, meta in regions_meta.items():
        region_len = meta[4]
        per_votu = aggregates.get(qid, {})
        if not per_votu:
            assignments[qid] = {
                "assigned_votu": "unassigned",
                "coverage": 0.0,
                "mean_identity": 0.0,
                "score": 0.0,
            }
            continue
        scored = []
        for votu, stats in per_votu.items():
            aln = stats["aln"]
            mean_id = stats["id_sum"] / aln if aln else 0.0
            coverage = aln / region_len
            agg_score = aln * mean_id / 100.0
            scored.append((agg_score, votu, coverage, mean_id))
        scored.sort(reverse=True)
        best_score, best_votu, best_cov, best_id = scored[0]
        if best_cov < min_coverage or best_id < min_identity * 100:
            best_votu = "unassigned"
        assignments[qid] = {
            "assigned_votu": best_votu,
            "coverage": round(best_cov, 3),
            "mean_identity": round(best_id, 2),
            "score": round(best_score, 1),
        }
    return assignments


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_mapping_table(path, regions_meta, assignments, phaster_meta):
    """Write the per-region mapping table."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "strain", "region", "start", "end", "length_bp",
            "phaster_completeness", "phaster_score", "phaster_top_hit",
            "assigned_votu", "coverage", "mean_identity", "score",
        ])
        for qid in sorted(regions_meta):
            strain, region_n, start, end, length_bp = regions_meta[qid]
            pm = phaster_meta.get((strain, region_n), {})
            a = assignments[qid]
            writer.writerow([
                strain, region_n, start, end, length_bp,
                pm.get("completeness", ""), pm.get("score", ""),
                pm.get("top_phaster_hit", ""),
                a["assigned_votu"], a["coverage"],
                a["mean_identity"], a["score"],
            ])


def write_phaster_crosstab(path, regions_meta, assignments, phaster_meta):
    """Write the PHASTER top-hit x assigned vOTU contingency table."""
    counts = defaultdict(lambda: defaultdict(int))
    for qid in regions_meta:
        strain, region_n, *_ = regions_meta[qid]
        pm = phaster_meta.get((strain, region_n), {})
        phaster_hit = pm.get("top_phaster_hit", "NA")
        votu = assignments[qid]["assigned_votu"]
        counts[phaster_hit][votu] += 1
    all_votus = sorted({v for row in counts.values() for v in row})
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["phaster_top_hit"] + all_votus + ["TOTAL"])
        for hit in sorted(counts):
            row = [hit] + [counts[hit].get(v, 0) for v in all_votus]
            row.append(sum(counts[hit].values()))
            writer.writerow(row)


def write_strain_x_votu_matrix(path, regions_meta, assignments):
    """Write the binary strain x vOTU presence/absence matrix."""
    strains = sorted({m[0] for m in regions_meta.values()})
    votus = sorted({
        a["assigned_votu"] for a in assignments.values()
        if a["assigned_votu"] != "unassigned"
    })
    matrix = {s: {v: 0 for v in votus} for s in strains}
    for qid, meta in regions_meta.items():
        strain = meta[0]
        votu = assignments[qid]["assigned_votu"]
        if votu in matrix[strain]:
            matrix[strain][votu] = 1
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["strain"] + votus)
        for strain in strains:
            writer.writerow([strain] + [matrix[strain][v] for v in votus])


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--phaster-dir", required=True,
                        help="Directory containing per-strain PHASTER output files.")
    parser.add_argument("--votu-dir", required=True,
                        help="Directory containing vOTU reference FASTAs.")
    parser.add_argument("--out-dir", required=True,
                        help="Output directory (created if absent).")
    parser.add_argument("--min-coverage", type=float, default=0.5,
                        help="Minimum aligned query fraction (default: 0.5).")
    parser.add_argument("--min-identity", type=float, default=0.75,
                        help="Minimum mean nucleotide identity (default: 0.75).")
    parser.add_argument("--task", default="dc-megablast",
                        help="BLASTn task (default: dc-megablast).")
    parser.add_argument("--threads", type=int, default=4,
                        help="BLAST thread count (default: 4).")
    args = parser.parse_args()

    for tool in ("makeblastdb", "blastn"):
        if shutil.which(tool) is None:
            sys.exit(f"ERROR: '{tool}' not found on PATH. Install NCBI BLAST+.")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    query_fasta = out_dir / "queries.fasta"
    regions_meta = assemble_query_fasta(args.phaster_dir, query_fasta)
    print(f"Assembled {len(regions_meta)} prophage regions.")

    target_fasta = out_dir / "votu_references.fasta"
    votu_counts = assemble_target_fasta(args.votu_dir, target_fasta)
    print(f"Assembled {len(votu_counts)} vOTU references "
          f"({sum(votu_counts.values())} sequences).")

    db_prefix = out_dir / "votu_db"
    run_makeblastdb(target_fasta, db_prefix)

    blast_tsv = out_dir / "blast_hits.tsv"
    run_blastn(query_fasta, db_prefix, blast_tsv,
               task=args.task, threads=args.threads)

    assignments = assign_votus(blast_tsv, regions_meta,
                               args.min_coverage, args.min_identity)
    phaster_meta = collect_phaster_metadata(args.phaster_dir)

    # Restrict outputs to intact PHASTER predictions, matching the standard
    # PHASTER best practice of excluding questionable and incomplete regions
    # from downstream presence/absence analyses.
    intact_keys = {
        qid for qid, meta in regions_meta.items()
        if phaster_meta.get((meta[0], meta[1]), {}).get("completeness") == "intact"
    }
    n_total = len(regions_meta)
    regions_meta = {k: v for k, v in regions_meta.items() if k in intact_keys}
    assignments = {k: v for k, v in assignments.items() if k in intact_keys}
    print(f"Retained {len(regions_meta)}/{n_total} intact prophage regions "
          f"for downstream vOTU mapping.")

    write_mapping_table(out_dir / "prophage_votu_mapping.tsv",
                        regions_meta, assignments, phaster_meta)
    write_phaster_crosstab(out_dir / "phaster_to_votu_crosstab.tsv",
                           regions_meta, assignments, phaster_meta)
    write_strain_x_votu_matrix(out_dir / "strain_x_votu_matrix.tsv",
                               regions_meta, assignments)

    n_assigned = sum(1 for a in assignments.values()
                     if a["assigned_votu"] != "unassigned")
    print(f"Assigned {n_assigned}/{len(assignments)} intact regions to vOTUs.")
    print(f"Outputs written to {out_dir}/")


if __name__ == "__main__":
    main()
