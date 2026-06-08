"""
S10A_GEM_analysis_part1.py — Bidirectional BLASTp and Homology Matrix Generation
-----------------------------------------------------------------------------
Author: Eva Stare

Adapted from the multi-strain GEM generation workflow by Norsigian et al. (2020),
originally developed for E. coli. We applied this workflow to B. subtilis using
the reference genome-scale metabolic model iBB1018 (Blázquez et al., 2023),
which is distributed in SBML (.xml) format rather than JSON. The script was
modified accordingly. Steps for downloading genomes from NCBI (step 1 of the
original notebook) and BLASTn-based curation of unannotated ORFs (step 4) were
omitted, as genome assemblies were obtained separately and annotation was
performed with Prokka. Additional minor modifications were made to improve
code organization and readability.

References:
  Norsigian CJ, Fang X, Seif Y, Monk JM, Palsson BO. A workflow for generating
  multi-strain genome-scale metabolic models of prokaryotes. Nat Protoc. 2020;
  15(1):1-14. doi:10.1038/s41596-019-0254-3

  Blázquez B, San León D, Rojas A, Tortajada M, Nogales J. New Insights on
  Metabolic Features of Bacillus subtilis Based on Multistrain Genome-Scale
  Metabolic Modeling. Int J Mol Sci. 2023;24:7091. doi:10.3390/ijms24087091

"""

import os
import sys
import pandas as pd
from glob import glob
from Bio import SeqIO
import cobra
import subprocess

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
WORK_DIR = "./GEM"
GENOMES_DIR = os.path.join(WORK_DIR, "genomes")
PROTS_DIR = os.path.join(WORK_DIR, "genomes", "prots")
NUCL_DIR = os.path.join(WORK_DIR, "genomes", "nucl")
BBH_DIR = os.path.join(WORK_DIR, "bbh")
REF_DIR = os.path.join(WORK_DIR, "genomes", "reference_strain")

REFERENCE_ID = "iBB1018"
REFERENCE_MODEL = os.path.join(WORK_DIR, "iBB1018.xml")
STRAIN_INFO = os.path.join(GENOMES_DIR, "StrainInformation.xlsx")

# BLAST parameters
EVALUE = 0.001
COVERAGE_THRESHOLD = 0.25
IDENTITY_THRESHOLD > 80.0

# Artificial genes to exclude from homology analysis
ARTIFICIAL_GENES = ["GROWMATCH", "GAPFILLING"]

os.makedirs(WORK_DIR, exist_ok=True)
for d in [GENOMES_DIR, PROTS_DIR, NUCL_DIR, BBH_DIR]:
    os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load strain information
# ---------------------------------------------------------------------------
strain_info = pd.read_excel(STRAIN_INFO)
target_strains = list(strain_info["NCBI ID"])
print(f"Target strains: {len(target_strains)}")

# ---------------------------------------------------------------------------
# 2. Extract protein and nucleotide sequences from GenBank files
# ---------------------------------------------------------------------------
def parse_genome(strain_id, seq_type="prot", in_folder=GENOMES_DIR,
                 out_folder=PROTS_DIR, ext="gbk"):
    """Parse GenBank file to produce protein or nucleotide FASTA."""
    os.makedirs(out_folder, exist_ok=True)
    in_file = os.path.join(in_folder, f"{strain_id}.{ext}")
    out_file = os.path.join(out_folder, f"{strain_id}.fa")

    if os.path.exists(out_file):
        return

    try:
        with open(in_file) as handle, open(out_file, "w") as fout:
            x = 0
            for record in SeqIO.parse(handle, "genbank"):
                for f in record.features:
                    if f.type == "CDS":
                        seq = f.extract(record.seq)
                        seq = str(seq) if seq_type == "nucl" else str(seq.translate())
                        if "locus_tag" in f.qualifiers:
                            locus = f.qualifiers["locus_tag"][0]
                        elif "gene" in f.qualifiers:
                            locus = f.qualifiers["gene"][0]
                        else:
                            locus = f"gene_{x}"
                            x += 1
                        fout.write(f">{locus}\n{seq}\n")
        print(f"Parsed {strain_id}")
    except Exception as e:
        print(f"Error parsing {strain_id}: {e}")


# Parse target strains (.gbk)
for strain in target_strains:
    parse_genome(strain, "prot", GENOMES_DIR, PROTS_DIR, "gbk")
    parse_genome(strain, "nucl", GENOMES_DIR, NUCL_DIR, "gbk")

# Parse reference strain (.gbff format)
parse_genome(REFERENCE_ID, "prot", REF_DIR, PROTS_DIR, "gbff")
parse_genome(REFERENCE_ID, "nucl", REF_DIR, NUCL_DIR, "gbff")

# ---------------------------------------------------------------------------
# 3. Build BLAST databases
# ---------------------------------------------------------------------------
def make_blast_db(strain_id, folder=PROTS_DIR, db_type="prot"):
    """Create BLAST database if it does not already exist."""
    ext = "fna" if db_type == "nucl" else "fa"
    db_check = os.path.join(folder, f"{strain_id}.fa.pin")
    if os.path.exists(db_check):
        return
    cmd = f"makeblastdb -in {folder}/{strain_id}.{ext} -dbtype {db_type}"
    os.system(cmd)


for strain in target_strains:
    make_blast_db(strain, PROTS_DIR, "prot")
make_blast_db(REFERENCE_ID, PROTS_DIR, "prot")

# ---------------------------------------------------------------------------
# 4. Run bidirectional BLASTp
# ---------------------------------------------------------------------------
def run_blastp(seq, db, in_folder=PROTS_DIR, out_folder=BBH_DIR,
               evalue=EVALUE, threads=1):
    """Run BLASTp if output does not already exist."""
    os.makedirs(out_folder, exist_ok=True)
    out_file = os.path.join(out_folder, f"{seq}_vs_{db}.txt")
    if os.path.exists(out_file):
        return out_file
    cmd = (f"blastp -db {in_folder}/{db}.fa -query {in_folder}/{seq}.fa "
           f"-out {out_file} -evalue {evalue} -outfmt 6 -num_threads {threads}")
    os.system(cmd)
    return out_file


def get_gene_lens(query, in_folder=PROTS_DIR):
    """Return DataFrame of gene lengths from a FASTA file."""
    records = SeqIO.parse(os.path.join(in_folder, f"{query}.fa"), "fasta")
    return pd.DataFrame([{"gene": r.name, "gene_length": len(r.seq)} for r in records])


def get_bbh(query, subject, in_folder=BBH_DIR, prot_folder=PROTS_DIR):
    """Compute bidirectional best BLAST hits between query and subject."""
    run_blastp(query, subject, prot_folder, in_folder)
    run_blastp(subject, query, prot_folder, in_folder)

    query_lengths = get_gene_lens(query, prot_folder)
    subject_lengths = get_gene_lens(subject, prot_folder)

    cols = ["gene", "subject", "PID", "alnLength", "mismatchCount",
            "gapOpenCount", "queryStart", "queryEnd", "subjectStart",
            "subjectEnd", "eVal", "bitScore"]

    bbh_fwd = pd.read_csv(os.path.join(in_folder, f"{query}_vs_{subject}.txt"),
                          sep="\t", names=cols)
    bbh_fwd = pd.merge(bbh_fwd, query_lengths)
    bbh_fwd["COV"] = bbh_fwd["alnLength"] / bbh_fwd["gene_length"]

    bbh_rev = pd.read_csv(os.path.join(in_folder, f"{subject}_vs_{query}.txt"),
                          sep="\t", names=cols)
    bbh_rev = pd.merge(bbh_rev, subject_lengths)
    bbh_rev["COV"] = bbh_rev["alnLength"] / bbh_rev["gene_length"]

    # Filter by coverage
    bbh_fwd = bbh_fwd[bbh_fwd.COV >= COVERAGE_THRESHOLD]
    bbh_rev = bbh_rev[bbh_rev.COV >= COVERAGE_THRESHOLD]

    out = pd.DataFrame()
    for g in bbh_fwd.gene.unique():
        res = bbh_fwd[bbh_fwd.gene == g]
        if len(res) == 0:
            continue
        best_hit = res.loc[res.PID.idxmax()].copy()
        best_gene = best_hit.subject
        res2 = bbh_rev[bbh_rev.gene == best_gene]
        if len(res2) == 0:
            continue
        best_hit2 = res2.loc[res2.PID.idxmax()].copy()
        best_hit["BBH"] = "<=>" if g == best_hit2.subject else "->"
        out = pd.concat([out, pd.DataFrame(best_hit).transpose()])

    out_file = os.path.join(in_folder, f"{query}_vs_{subject}_parsed.csv")
    out.to_csv(out_file, index=False)
    print(f"BBH: {query} vs {subject} -> {len(out)} hits")


# Run BBH for all target strains
for strain in target_strains:
    get_bbh(REFERENCE_ID, strain, BBH_DIR, PROTS_DIR)

# ---------------------------------------------------------------------------
# 5. Build homology and gene ID matrices
# ---------------------------------------------------------------------------
model = cobra.io.read_sbml_model(REFERENCE_MODEL)
model_genes = [g.id for g in model.genes if g.id not in ARTIFICIAL_GENES]

ortho_matrix = pd.DataFrame(index=model_genes, columns=target_strains)
geneIDs_matrix = pd.DataFrame(index=model_genes, columns=target_strains)

blast_files = glob(os.path.join(BBH_DIR, "*_parsed.csv"))

for blast_file in blast_files:
    bbh = pd.read_csv(blast_file)
    ids, pids = [], []
    for gene in model_genes:
        match = bbh[bbh["gene"] == gene]
        if len(match) > 0:
            row = match.iloc[0]
            ids.append(row["subject"])
            pids.append(row["PID"])
        else:
            ids.append("None")
            pids.append(0)
    for col in ortho_matrix.columns:
        if col in blast_file:
            ortho_matrix[col] = pids
            geneIDs_matrix[col] = ids

# ---------------------------------------------------------------------------
# 6. Binarize at identity threshold
# ---------------------------------------------------------------------------
for col in ortho_matrix.columns:
    ortho_matrix.loc[ortho_matrix[col] <= IDENTITY_THRESHOLD, col] = 0
    ortho_matrix.loc[ortho_matrix[col] > IDENTITY_THRESHOLD, col] = 1

# ---------------------------------------------------------------------------
# 7. Save outputs
# ---------------------------------------------------------------------------
ortho_matrix.to_csv(os.path.join(WORK_DIR, "ortho_matrix.csv"))
geneIDs_matrix.to_csv(os.path.join(WORK_DIR, "geneIDs_matrix.csv"))

print(f"\nOrtho matrix: {ortho_matrix.shape}")
print(f"Gene IDs matrix: {geneIDs_matrix.shape}")
print("Done.")
