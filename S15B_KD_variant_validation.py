"""
# KD Group Variant Validation Pipeline
## Pairwise Protein Alignment of Gene Variants vs Reference

**Author:** Eva Stare

This notebook validates group-specific gene variants identified in the kin discrimination
analysis by performing pairwise local protein alignments (Smith-Waterman, BLOSUM62) against
reference variants.

**Goal:** For each group-specific variant (e.g., `tagF.1` in g50), extract the
protein sequence from the group strain's GFF + genome, compare it against the
reference variant (`.0`, typically from NCIB_3610), and report % amino acid identity.

**Flags:**
- **>95% AA identity** → Potential clustering artefact (like g30's yybE at 99%)
- **<80% AA identity** → Strong variant divergence (like g50's tagF at 49%)
- **Split/truncated** → Gene integrity issues

"""
# What the script does

### Step 1: Setup (Cells 1-2)

Installs BioPython and pandas. Mounts Google Drive and sets paths to input files: megamatrix Excel, features analysis Excel, and the folder of GBK files.

### Step 2: Map strain names to GBK filenames (Cell 3)

67 strains have different naming conventions: PS-xxx strains use their strain name directly (e.g., `PS-130_edit2.gbk`), while others use GCF accessions (e.g., `NCIB_3610` → `GCF-029027845.1_edit2.gbk`). This cell contains the hardcoded mapping and verifies all 67 GBK files are found.

### Step 3: Parse GBK files (Cell 4)

For each of the 67 strains, reads the GenBank file using BioPython's `SeqIO` and builds a dictionary: `{gene_name → [{locus_tag, protein_sequence, start, end, strand, length}, ...]}`. Gene names are extracted from the `/gene` qualifier first, falling back to `[gene=XXX]` in the `/product` field (Prokka annotation format). A gene can have multiple entries if it appears as multiple copies in the genome. Runs a verification check on known genes in NCIB_3610.

### Step 4: Build validation tasks (Cell 5)

For each group-specific variant (e.g., `tagF.1` in g50) from the features analysis Excel:

- Strips the variant number to get the base gene name (`tagF.1` → `tagF`)
- Looks up all variants of that gene in the megamatrix (e.g., `tagF.0`, `.1`, `.2`, etc.)
- Picks the reference variant to compare against, in order of preference:
  1. `.0` from NCIB_3610
  2. `.0` from whichever strain carries it
  3. Lowest-numbered variant NCIB_3610 carries
  4. **If no `.0` exists in the megamatrix:** uses the **base gene name** (e.g., `tagF`) from NCIB_3610 as reference — the absence of `.0` means Roary did not split out a separate `.0` cluster, but the gene itself is still present in genomes under its unsuffixed name. This comparison is methodologically identical to the `.0` case.
  5. If the chosen reference equals the current variant: falls back to the most common other variant, or the base gene from NCIB_3610
- For multi-strain groups (g50, g31, g55, etc.), creates separate tasks for each group strain (up to 3), so results can be averaged
- Multi-name genes like `"tagF, rodC"` are split on commas: tries `tagF` first in the GBK, stops if found, otherwise tries `rodC`

### Step 5: Extract proteins (Cell 6)

Helper function `get_protein_for_gene`: given a strain and gene name(s), finds all matching protein sequences in the parsed GBK data. Tries exact match first. Returns all copies if the gene is duplicated.

### Step 6: Pairwise alignment (Cell 7)

For each validation task:

- Computes all pairwise combinations between reference strain copies and group strain copies (e.g., if ref has 1 copy and group has 3 → 3 alignments)
- Uses BLOSUM62 local alignment (BioPython `pairwise2`)
- Calculates: % AA identity, alignment length, length ratio (shorter/longer × 100)
- From all pairs, selects the **most divergent pair** as the primary result (lowest identity), because that's most likely the actual group-specific variant rather than a conserved copy matching `.0`
- Records all pair identities in `all_pairs_identities` and `all_pairs_detail` columns, plus `identity_spread` (max − min identity across all pairs)
- Reports multi-copy cases at the end for review

### Step 7: Average and classify (Cell 8)

For multi-strain groups, per-strain results are averaged: mean, min, and max identity across the 2-3 group strains tested. Then each variant is classified into three categories:

| Classification | Criteria | Meaning |
|:---|:---|:---|
| **LIKELY ARTEFACT** | ≥95% AA identity AND ≥80% length ratio | Roary clustering error — should have been merged (Roary threshold was 80%) |
| **LIKELY TRUNCATED** | <80% length ratio | Assembly fragment or partial annotation, not a full-length variant |
| **VARIANT** | <95% AA identity AND ≥80% length ratio | Genuinely divergent protein, correctly split by Roary |

### Step 8: Show failed lookups (Cell 9)

Lists variants where the gene couldn't be found in the GBK file — these need manual checking. Common causes: gene was manually renamed in the GBK differently from the megamatrix, or the gene annotation uses a synonym.

### Step 9: Save validation results (Cell 10)

Writes `variant_validation_results.xlsx` with sheets:

- **All Results**: averaged across strains, one row per variant
- **Raw Per-Strain**: individual alignments for each group strain
- **Likely Artefacts / Likely Truncated / Variants**: filtered sheets by classification
- **Group Summary**: per-group counts of each classification

### Step 10: Visualization (Cell 11)

Histogram of all variant identities with threshold lines, plus per-group boxplot.

### Step 11: Generate corrected features analysis (Cells 12a-c)

Reads the validation results and the original features analysis, then for each group:

- Annotates every variant with its AA identity and classification
- Flags artefacts for removal (≥95% identity with ≥80% length ratio)
- Safeguards multi-copy genes: if a variant was flagged as artefact but has `identity_spread` >10% across copies, it's kept (the high identity is from comparing the wrong copy)
- Removes paired absent variants: if present variant `tagF.1` is an artefact, the corresponding absent `tagF.0` is also removed
- Non-variant features (original presence/absence genes like `kdpA`) are kept untouched — they don't have a `.0` to compare against
- Outputs `corrected_features_analysis.xlsx` with the same group-sheet structure as the original, plus two new columns: `AA% vs ref` and `status`

### Step 12: Quick overview (Cell 13)

Prints: (a) exactly which features were removed from each group, and (b) all validated variants per group with their identity percentages, by classification.
"""

"""
## 1. Setup and Installation
"""

# Install required packages

import pandas as pd
import re
import os
import glob
from collections import defaultdict, Counter
from Bio import pairwise2
from Bio.Align import substitution_matrices
import warnings
warnings.filterwarnings('ignore')

print("Setup complete!")

"""
## 2. Mount Google Drive and Set Paths

Upload your files to Google Drive and set the paths below.
"""

# ============================================================
# SET YOUR PATHS HERE
# ============================================================
BASE_DIR = 'input/'  # UPDATE THIS PATH  # Change this!

MEGAMATRIX_PATH = os.path.join(BASE_DIR, '67_only_megamatrix_filtered_strains.xlsx')
FEATURES_PATH = os.path.join(BASE_DIR, 'Features_analysis_67_output.xlsx')
GBK_DIR = os.path.join(BASE_DIR, '67_strains_gbks')    # Folder with .gbk files

# Output
OUTPUT_DIR = os.path.join(BASE_DIR, 'validation_results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"Base dir: {BASE_DIR}")
print(f"GBK dir exists: {os.path.isdir(GBK_DIR)}")
print(f"Megamatrix exists: {os.path.isfile(MEGAMATRIX_PATH)}")
print(f"Features exists: {os.path.isfile(FEATURES_PATH)}")

"""
## 3. Map Strain Names to GBK File Names

All GBK files are named `<base>_edit2.gbk`.
PS-xxx strains use their strain name; others use GCF accessions.
**Edit `STRAIN_TO_FILENAME` if any files are named differently.**
"""

# Load strain names from megamatrix
df_mega = pd.read_excel(MEGAMATRIX_PATH, sheet_name='Megamatrix', header=None)
strain_names = [str(df_mega.iloc[i, 0]) for i in range(1, len(df_mega))]
print(f"Total strains in megamatrix: {len(strain_names)}")

# ============================================================
# STRAIN NAME -> FILE BASE NAME (without _edit2.gbk)
# ============================================================
STRAIN_TO_FILENAME = {
    '73': 'GCF-009662275.1',
    'BS16045': 'GCF-001720505.1',
    'KF24': 'GCF-030123145.1',
    'MB8_B1': 'GCF-009662255.1',
    'MB8_B7': 'GCF-009662215.1',
    'MB9_B1': 'GCF-009662175.1',
    'MB9_B6': 'GCF-009662375.1',
    'NCIB_3610': 'GCF-029027845.1',
    'NRS6085': 'GCF-905310975.2',
    'NRS6099': 'GCF-905310985.2',
    'NRS6103': 'GCF-905310995.2',
    'NRS6105': 'GCF-905311425.2',
    'NRS6118': 'GCF-905311435.2',
    'NRS6121': 'GCF-905315035.2',
    'NRS6127': 'GCF-905315045.2',
    'NRS6128': 'GCF-905315055.2',
    'NRS6132': 'GCF-905315015.2',
    'NRS6145': 'GCF-905315685.2',
    'NRS6153': 'GCF-905315705.2',
    'NRS6160': 'GCF-905315715.2',
    'NRS6181': 'GCF-905318255.2',
    'NRS6186': 'GCF-905319135.2',
    'NRS6187': 'GCF-905319545.2',
    'NRS6202': 'GCF-905319535.2',
    'P8_B1': 'GCF-009662435.1',
    'P9_B1': 'GCF-009662455.1',
    'RO-A-4': 'RO-A-4',
    'RO-DD-2': 'RO-DD-2',
    'RO-F-3': 'RO-F-3',
    'RO-FF-1': 'RO-FF-1',
    'RS-D-2': 'RS-D-2',
}
# All PS-xxx strains: filename = strain name
for s in strain_names:
    if s not in STRAIN_TO_FILENAME:
        STRAIN_TO_FILENAME[s] = s

# ============================================================
# All files are named: <base>_edit2.gbk
# ============================================================
gbk_map = {}  # strain_name -> gbk_path

for strain in strain_names:
    fname = STRAIN_TO_FILENAME.get(strain, strain)
    path = os.path.join(GBK_DIR, f"{fname}_edit2.gbk")
    if os.path.isfile(path):
        gbk_map[strain] = path
    else:
        # Fallback: glob
        matches = glob.glob(os.path.join(GBK_DIR, f"{fname}*.gbk"))
        if matches:
            gbk_map[strain] = matches[0]

print(f"GBK files found: {len(gbk_map)} / {len(strain_names)}")

missing = [s for s in strain_names if s not in gbk_map]
if missing:
    print(f"\n  Missing GBK ({len(missing)}):")
    for s in missing:
        expected = STRAIN_TO_FILENAME.get(s, s)
        print(f"    {s} -> expected: {expected}_edit2.gbk")
else:
    print("\n All GBK files found!")

# Sanity check
print("\nSample file paths:")
for strain in list(gbk_map.keys())[:5]:
    print(f"  {strain:15s} -> {os.path.basename(gbk_map[strain])}")

"""
## 4. Parse GBK Files → Gene Names + Protein Sequences

Each GBK file contains gene annotations with translated protein sequences.
Gene names are extracted from the `/gene` or `/product` qualifiers.
This builds: `{strain: {gene_name: [(locus_tag, protein_seq, start, end, strand), ...]}}`
"""

from Bio import SeqIO

def parse_gbk(gbk_path):
    """
    Parse a GenBank file and return:
      gene_data: {gene_name: [{'locus_tag': ..., 'translation': ...,
                                'start': ..., 'end': ..., 'strand': ..., 'length_bp': ...}, ...]}

    Handles Prokka GBK format where gene names may be in:
      - /gene qualifier
      - /product qualifier as [gene%3DgeneName] (URL-encoded)
      - /product qualifier as plain text containing gene name
    """
    gene_data = defaultdict(list)

    for record in SeqIO.parse(gbk_path, 'genbank'):
        for feature in record.features:
            if feature.type != 'CDS':
                continue

            qualifiers = feature.qualifiers

            # Get locus_tag
            locus_tag = qualifiers.get('locus_tag', [''])[0]
            if not locus_tag:
                continue

            # Get translation
            translation = qualifiers.get('translation', [''])[0]
            if not translation:
                continue

            # Get gene name - try multiple sources
            gene_name = None

            # Source 1: /gene qualifier
            if 'gene' in qualifiers:
                gene_name = qualifiers['gene'][0]

            # Source 2: product field with [gene%3DgeneName]
            if not gene_name and 'product' in qualifiers:
                product = qualifiers['product'][0]
                match = re.search(r'gene%3D([^\]\s%]+)', product)
                if match:
                    gene_name = match.group(1)
                # Also try URL-decoded: [gene=geneName]
                if not gene_name:
                    match = re.search(r'\[gene=([^\]]+)\]', product)
                    if match:
                        gene_name = match.group(1)

            if gene_name:
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = '+' if feature.location.strand == 1 else '-'

                gene_data[gene_name].append({
                    'locus_tag': locus_tag,
                    'translation': translation,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'length_bp': end - start,
                    'length_aa': len(translation),
                })

    return dict(gene_data)


# Parse all GBK files
print("Parsing GBK files (this may take a few minutes)...")
strain_data = {}  # {strain: {gene_name: [entries]}}

for i, (strain, gbk_path) in enumerate(gbk_map.items()):
    if i % 10 == 0:
        print(f"  {i}/{len(gbk_map)}: {strain}...")
    strain_data[strain] = parse_gbk(gbk_path)

print(f"\nParsed {len(strain_data)} GBK files")

# Quick stats
for strain in list(strain_data.keys())[:5]:
    n_genes = len(strain_data[strain])
    n_prots = sum(len(v) for v in strain_data[strain].values())
    print(f"  {strain}: {n_genes} gene names, {n_prots} CDS entries")

# Verify known genes in NCIB_3610
if 'NCIB_3610' in strain_data:
    print(f"\nGene check in NCIB_3610:")
    for g in ['tagA', 'tagF', 'dnaA', 'sigV', 'conJ', 'comQ', 'rapI']:
        if g in strain_data['NCIB_3610']:
            entries = strain_data['NCIB_3610'][g]
            for e in entries:
                print(f"   {g}: {e['locus_tag']}, {e['length_aa']} AA, {e['start']}-{e['end']} ({e['strand']})")
        else:
            print(f"   {g}: NOT FOUND — check GBK annotation")

"""
## 5. Build Validation Task List

For each group-specific variant, determine:
- Reference variant (`.0` or lowest common variant)
- Reference strain (preferably NCIB_3610)
- Group strain carrying the variant
"""

# Load megamatrix
df_mega = pd.read_excel(MEGAMATRIX_PATH, sheet_name='Megamatrix', header=None)
strain_list = [str(df_mega.iloc[i, 0]) for i in range(1, len(df_mega))]
gene_list = [str(df_mega.iloc[0, j]) for j in range(1, df_mega.shape[1])]

ncib_row = strain_list.index('NCIB_3610') + 1  # +1 for header row

# Build gene->column mapping
gene_to_col = {}
for j, g in enumerate(gene_list):
    gene_to_col[g] = j + 1  # +1 for strain name column

# ============================================================
# GROUP -> STRAIN MAPPING
# Edit this if needed!
# For singleton groups, just one strain.
# For multi-strain groups, first strain is used for comparison.
# ============================================================
GROUP_STRAINS = {
    'g3': ['RO-F-3'],
    'g4': ['RS-D-2'],
    'g5': ['RO-DD-2'],
    'g7': ['RO-A-4'],
    'g19': ['NRS6085'],
    'g25_A': ['PS-216', 'MB8_B7', 'P8_B1'],  # 3 representatives of 21
    'g29': ['BS16045'],
    'g30': ['NRS6105', 'NRS6145'],  # both strains
    'g31': ['PS-52', 'PS-53'],  # both strains
    'g34': ['NRS6128'],
    'g39': ['PS-209'],
    'g44': ['PS-196'],
    'g50': ['PS-130', 'PS-108', 'PS-119'],  # 3 of 5
    'g52': ['NRS6160'],
    'g53': ['PS-24', 'PS-25', 'PS-20'],  # all 3
    'g54': ['PS-160'],
    'g55': ['PS-93', 'PS-95'],  # both strains
    'g56': ['PS-217', 'PS-218'],  # both strains
    'g57': ['PS-15'],
    'g60_A': ['PS-194', 'PS-64', 'NRS6186'],  # 3 of 4
    'g60_B': ['NRS6132', 'NRS6099', 'NRS6187'],  # all 3
    'g63_A': ['PS-149', 'NRS6127', 'PS-55'],  # 3 of 7
    'g82': ['RO-FF-1'],
    'g92': ['KF24'],
    'g95': ['PS-263'],
    'g96': ['73'],
}

# Load features
xls_feat = pd.ExcelFile(FEATURES_PATH)

validation_tasks = []

for sheet in xls_feat.sheet_names:
    if sheet.lower() == 'removed':
        continue
    df_feat = pd.read_excel(xls_feat, sheet_name=sheet, header=None)

    for i in range(2, len(df_feat)):
        gene = df_feat.iloc[i, 0]
        if not pd.notna(gene) or str(gene).strip() in ['/', '']:
            continue
        gene_str = str(gene).strip()
        if not re.search(r'\.\d+$', gene_str):
            continue  # skip non-variants

        base = re.sub(r'\.\d+$', '', gene_str)
        var_num = int(re.search(r'\.(\d+)$', gene_str).group(1))

        # Find all variants of this gene in the matrix
        all_vars = {}
        for gname, col in gene_to_col.items():
            if gname.startswith(base + '.') and re.search(r'\.\d+$', gname):
                vn = int(re.search(r'\.(\d+)$', gname).group(1))
                count = sum(1 for r in range(1, len(df_mega))
                           if str(df_mega.iloc[r, col]) == '1')
                ncib_has = str(df_mega.iloc[ncib_row, col]) == '1'
                all_vars[gname] = {'var_num': vn, 'count': count, 'ncib': ncib_has}

        # Pick reference variant
        ref_var = None
        ref_strain = None
        no_dot_zero = False  # flag: .0 absent from megamatrix

        # Prefer .0 with NCIB
        v0 = base + '.0'
        if v0 in all_vars and all_vars[v0]['ncib']:
            ref_var = v0
            ref_strain = 'NCIB_3610'
        elif v0 in all_vars:
            ref_var = v0
            col = gene_to_col[v0]
            for r in range(1, len(df_mega)):
                if str(df_mega.iloc[r, col]) == '1':
                    ref_strain = str(df_mega.iloc[r, 0])
                    break
        else:
            # ==========================================================
            # No .0 in megamatrix. The base gene (e.g. "tagF" without
            # variant suffix) was not split by Roary into a .0 cluster,
            # but the gene is still present in genomes under its base
            # name. Use NCIB_3610 as reference with the base gene name
            # — methodologically identical to the .0 case: we compare
            # all copies in the reference vs all copies in the group
            # strain using the same pairwise alignment approach.
            # ==========================================================
            ncib_vars = {k: v for k, v in all_vars.items() if v['ncib']}
            if ncib_vars:
                ref_var = min(ncib_vars.keys(), key=lambda x: all_vars[x]['var_num'])
                ref_strain = 'NCIB_3610'
            else:
                # NCIB doesn't carry any numbered variant -> use base
                # gene name from NCIB_3610 as reference
                ref_var = base + ' (base gene, no .0 in matrix)'
                ref_strain = 'NCIB_3610'
                no_dot_zero = True

            # If ref_var is the current variant itself, fall back
            if not no_dot_zero and ref_var == gene_str:
                other_vars = {k: v for k, v in all_vars.items() if k != gene_str}
                if other_vars:
                    ref_var = max(other_vars.keys(), key=lambda x: other_vars[x]['count'])
                    col = gene_to_col[ref_var]
                    for r in range(1, len(df_mega)):
                        if str(df_mega.iloc[r, col]) == '1':
                            ref_strain = str(df_mega.iloc[r, 0])
                            break
                else:
                    ref_var = base + ' (base gene, no .0 in matrix)'
                    ref_strain = 'NCIB_3610'
                    no_dot_zero = True

        if not no_dot_zero and (ref_var is None or ref_var == gene_str):
            continue

        # Extract gene names for GFF lookup
        # Handle multi-name genes like "tagF, rodC"
        gene_names_to_try = [n.strip() for n in base.split(',')]

        grp_strains = GROUP_STRAINS.get(sheet, [])
        if not grp_strains:
            continue

        func = str(df_feat.iloc[i, 1]).strip() if pd.notna(df_feat.iloc[i, 1]) else '?'

        # Create a task for EACH group strain (up to 3)
        for grp_strain in grp_strains:
            validation_tasks.append({
                'group': sheet,
                'variant': gene_str,
                'ref_variant': ref_var,
                'ref_strain': ref_strain,
                'group_strain': grp_strain,
                'gene_names': gene_names_to_try,
                'base': base,
                'function': func[:80],
                'no_dot_zero': no_dot_zero,
            })

print(f"Total validation tasks: {len(validation_tasks)}")
n_no0 = sum(1 for t in validation_tasks if t.get('no_dot_zero', False))
print(f"  Standard (.0 reference): {len(validation_tasks) - n_no0}")
print(f"  No .0 in matrix (base gene from NCIB_3610): {n_no0}")
print(f"\nTasks per group:")
for g, c in sorted(Counter(v['group'] for v in validation_tasks).items()):
    n0 = sum(1 for t in validation_tasks if t['group'] == g and t.get('no_dot_zero', False))
    suffix = f"  ({n0} via base gene)" if n0 > 0 else ""
    print(f"  {g}: {c}{suffix}")

"""
## 6. Helper: Extract Protein for a Gene

Simple lookup: gene name → protein sequence from the parsed GBK data.
"""

def get_protein_for_gene(strain, gene_names, strain_data):
    """
    Given a strain and list of possible gene names, find protein sequence(s).
    Returns list of dicts with locus_tag, sequence, length_aa, start, end, strand.
    """
    if strain not in strain_data:
        return []

    gdata = strain_data[strain]
    results = []
    seen_tags = set()

    for gname in gene_names:
        gname = gname.strip()
        if not gname:
            continue

        # Exact match — stop at first name that hits
        if gname in gdata:
            for entry in gdata[gname]:
                if entry['locus_tag'] not in seen_tags:
                    seen_tags.add(entry['locus_tag'])
                    results.append({
                        'locus_tag': entry['locus_tag'],
                        'sequence': entry['translation'],
                        'length_aa': entry['length_aa'],
                        'length_bp': entry['length_bp'],
                        'start': entry['start'],
                        'end': entry['end'],
                        'strand': entry['strand'],
                        'gene_name': gname,
                    })
            break
    return results


# Test with a known example
if 'NCIB_3610' in strain_data:
    test = get_protein_for_gene('NCIB_3610', ['tagA'], strain_data)
    if test:
        print(f"Test: tagA in NCIB_3610")
        for t in test:
            print(f"  {t['locus_tag']}: {t['length_aa']} AA, {t['length_bp']} bp")
            print(f"  Sequence: {t['sequence'][:50]}...")
    else:
        print("WARNING: Could not find tagA in NCIB_3610!")
        print("Available gene names (first 20):")
        print(list(strain_data['NCIB_3610'].keys())[:20])

"""
## 7. Run Pairwise Local Alignments

For each validation task, align the reference protein against the group-specific variant protein
and calculate % amino acid identity and % coverage.
"""

def calculate_aa_identity(seq1, seq2):
    """
    Calculate % AA identity and % coverage using LOCAL alignment.
    seq1 = reference (.0 variant)
    seq2 = group variant

    Returns dict with:
      - pct_identity: % identical positions over aligned region
      - coverage_ref: % of reference covered by alignment
      - coverage_query: % of query covered by alignment
      - aligned_length: length of aligned region (excluding terminal gaps)
      - identical_positions: number of identical AAs
      - ref_length: full length of reference
      - query_length: full length of query
    """
    if not seq1 or not seq2:
        return None

    try:
        matrix = substitution_matrices.load("BLOSUM62")
        alignments = pairwise2.align.localds(
            seq1, seq2, matrix, -10, -0.5,
            one_alignment_only=True
        )
    except Exception:
        alignments = pairwise2.align.localxx(
            seq1, seq2,
            one_alignment_only=True
        )

    if not alignments:
        return None

    aln = alignments[0]
    aligned_seq1 = aln.seqA
    aligned_seq2 = aln.seqB

    # Find the aligned region (local alignment pads with gaps outside)
    # Find first and last position where both are non-gap
    first_pos = None
    last_pos = None
    for idx, (a, b) in enumerate(zip(aligned_seq1, aligned_seq2)):
        if a != '-' and b != '-':
            if first_pos is None:
                first_pos = idx
            last_pos = idx

    if first_pos is None:
        return None

    # Count matches and alignment length in the aligned region only
    matches = 0
    aligned_length = 0
    ref_aligned = 0    # how many ref residues are in the alignment
    query_aligned = 0  # how many query residues are in the alignment

    for idx in range(first_pos, last_pos + 1):
        a = aligned_seq1[idx]
        b = aligned_seq2[idx]
        if a != '-' or b != '-':
            aligned_length += 1
        if a != '-':
            ref_aligned += 1
        if b != '-':
            query_aligned += 1
        if a == b and a != '-':
            matches += 1

    pct_identity = (matches / aligned_length * 100) if aligned_length > 0 else 0
    coverage_ref = (ref_aligned / len(seq1) * 100) if len(seq1) > 0 else 0
    coverage_query = (query_aligned / len(seq2) * 100) if len(seq2) > 0 else 0

    return {
        'pct_identity': round(pct_identity, 2),
        'coverage_ref': round(coverage_ref, 2),
        'coverage_query': round(coverage_query, 2),
        'aligned_length': aligned_length,
        'identical_positions': matches,
        'ref_aligned_aa': ref_aligned,
        'query_aligned_aa': query_aligned,
        'ref_length': len(seq1),
        'query_length': len(seq2),
    }


# Run all validations
print(f"Running {len(validation_tasks)} pairwise alignments...\n")

results = []
failed = []

for idx, task in enumerate(validation_tasks):
    if idx % 50 == 0:
        print(f"  Progress: {idx}/{len(validation_tasks)}...")

    # For all tasks (with or without .0): look up by base gene name.
    # When .0 exists, gene_names maps to the right gene in the ref.
    # When no .0 exists, this finds all copies of the base gene in
    # both NCIB_3610 and the group strain — same pairwise approach.
    ref_proteins = get_protein_for_gene(
        task['ref_strain'], task['gene_names'],
        strain_data
    )

    grp_proteins = get_protein_for_gene(
        task['group_strain'], task['gene_names'],
        strain_data
    )

    if not ref_proteins:
        failed.append({**task, 'reason': f"Gene not found in ref strain {task['ref_strain']}"})
        continue
    if not grp_proteins:
        failed.append({**task, 'reason': f"Gene not found in group strain {task['group_strain']}"})
        continue

    # Compute ALL pairwise combinations
    all_pairs = []

    for ref_p in ref_proteins:
        for grp_p in grp_proteins:
            aln_result = calculate_aa_identity(
                ref_p['sequence'], grp_p['sequence']
            )

            if aln_result is not None:
                all_pairs.append({
                    'group': task['group'],
                    'variant': task['variant'],
                    'ref_variant': task['ref_variant'],
                    'function': task['function'],
                    'ref_strain': task['ref_strain'],
                    'ref_locus_tag': ref_p['locus_tag'],
                    'ref_length_aa': ref_p['length_aa'],
                    'ref_location': f"{ref_p['start']}-{ref_p['end']}({ref_p['strand']})",
                    'group_strain': task['group_strain'],
                    'group_locus_tag': grp_p['locus_tag'],
                    'group_length_aa': grp_p['length_aa'],
                    'group_location': f"{grp_p['start']}-{grp_p['end']}({grp_p['strand']})",
                    'pct_identity': aln_result['pct_identity'],
                    'coverage_ref': aln_result['coverage_ref'],
                    'coverage_query': aln_result['coverage_query'],
                    'aligned_length': aln_result['aligned_length'],
                    'identical_positions': aln_result['identical_positions'],
                    'length_ratio': round(min(ref_p['length_aa'], grp_p['length_aa']) /
                                         max(ref_p['length_aa'], grp_p['length_aa']) * 100, 1)
                                         if max(ref_p['length_aa'], grp_p['length_aa']) > 0 else 0,
                    'n_ref_copies': len(ref_proteins),
                    'n_group_copies': len(grp_proteins),
                })

    if not all_pairs:
        continue

    # Sort by identity: lowest first (most divergent = likely the actual variant)
    all_pairs.sort(key=lambda x: x['pct_identity'])

    if len(all_pairs) == 1:
        results.append(all_pairs[0])
    else:
        # MULTI-COPY: pick most divergent as primary
        primary = all_pairs[0]
        primary['all_pairs_identities'] = ', '.join(
            f"{p['pct_identity']}%" for p in all_pairs
        )
        primary['all_pairs_detail'] = ' | '.join(
            f"ref:{p['ref_locus_tag']} vs grp:{p['group_locus_tag']} = "
            f"{p['pct_identity']}% id, {p['coverage_ref']}% cov_ref, {p['coverage_query']}% cov_qry "
            f"({p['ref_length_aa']}AA vs {p['group_length_aa']}AA)"
            for p in all_pairs
        )
        primary['best_pair_identity'] = all_pairs[-1]['pct_identity']
        primary['worst_pair_identity'] = all_pairs[0]['pct_identity']
        primary['identity_spread'] = round(all_pairs[-1]['pct_identity'] - all_pairs[0]['pct_identity'], 2)
        results.append(primary)

print(f"\nCompleted: {len(results)} alignments")
print(f"Failed: {len(failed)} (gene not found in GBK)")

# Show multi-copy cases
multi = [r for r in results if r.get('n_ref_copies', 1) > 1 or r.get('n_group_copies', 1) > 1]
if multi:
    print(f"\n Multi-copy genes detected: {len(multi)}")
    for r in multi[:10]:
        spread = r.get('identity_spread', 0)
        print(f"  {r['group']:8s} {r['variant']:35s} "
              f"copies: ref={r['n_ref_copies']}, grp={r['n_group_copies']} | "
              f"identities: {r.get('all_pairs_identities', '?')} "
              f"{' SPREAD=' + str(spread) + '%' if spread > 10 else ''}")
    if len(multi) > 10:
        print(f"  ... and {len(multi) - 10} more")

"""
## 8. Results Summary and Flagging
"""

df_results = pd.DataFrame(results)

# Average across group strains for the same variant
if len(df_results) > 0:
    # First, compute per-variant averages across group strains
    avg_cols = ['group', 'variant', 'ref_variant', 'function', 'ref_strain']
    df_avg = df_results.groupby(['group', 'variant']).agg(
        ref_variant=('ref_variant', 'first'),
        function=('function', 'first'),
        ref_strain=('ref_strain', 'first'),
        ref_length_aa=('ref_length_aa', 'first'),
        mean_pct_identity=('pct_identity', 'mean'),
        min_pct_identity=('pct_identity', 'min'),
        max_pct_identity=('pct_identity', 'max'),
        mean_coverage_ref=('coverage_ref', 'mean'),
        mean_coverage_query=('coverage_query', 'mean'),
        n_strains_compared=('pct_identity', 'count'),
        group_length_aa_avg=('group_length_aa', 'mean'),
        length_ratio_avg=('length_ratio', 'mean'),
    ).reset_index()

    df_avg['pct_identity'] = df_avg['mean_pct_identity'].round(2)
    df_avg['coverage_ref'] = df_avg['mean_coverage_ref'].round(2)
    df_avg['coverage_query'] = df_avg['mean_coverage_query'].round(2)

    print(f"Raw alignments: {len(df_results)}")
    print(f"Unique variants (averaged across strains): {len(df_avg)}")

    # Switch to averaged results for downstream analysis
    df_results_raw = df_results.copy()  # keep raw
    df_results = df_avg

if len(df_results) > 0:
    # ============================================================
    # CLASSIFICATION THRESHOLDS
    # ============================================================
    # Three simple categories:
    #   1. Likely artefact:  ≥95% AA identity AND ≥80% length ratio
    #      → Roary should have merged these (its threshold was 80%)
    #   2. Likely truncated: <80% length ratio
    #      → assembly fragment, not a real full-length variant
    #   3. True variant:     ≥80% length ratio AND <95% AA identity
    #      → genuinely divergent protein, Roary correctly split it
    ARTEFACT_THRESHOLD = 95    # ≥95% AA identity → likely clustering artefact
    LENGTH_RATIO_THRESHOLD = 80  # <80% length ratio → likely truncated

    # ============================================================
    # CLASSIFY EACH VARIANT
    # ============================================================
    def classify_variant(row):
        pct_id = row['pct_identity']
        len_ratio = row['length_ratio_avg']

        # TRUNCATED: length ratio <80% → partial/fragmented gene
        if len_ratio < LENGTH_RATIO_THRESHOLD:
            return f'LIKELY TRUNCATED ({len_ratio:.0f}% length ratio)'

        # Full-length comparisons
        if pct_id >= ARTEFACT_THRESHOLD:
            return 'LIKELY ARTEFACT (≥95% AA identity)'
        else:
            return f'VARIANT (<95% AA identity)'

    df_results['classification'] = df_results.apply(classify_variant, axis=1)

    # Multi-strain flag
    df_results.loc[df_results['n_strains_compared'] > 1, 'multi_strain'] = True
    df_results['multi_strain'] = df_results['multi_strain'].fillna(False)

    # ============================================================
    # SUMMARY STATISTICS
    # ============================================================
    print("=" * 80)
    print("VARIANT VALIDATION SUMMARY")
    print("=" * 80)
    print(f"\nTotal unique variants validated: {len(df_results)}")
    print(f"Mean % AA identity: {df_results['pct_identity'].mean():.1f}%")
    print(f"Median % AA identity: {df_results['pct_identity'].median():.1f}%")

    print(f"\nClassification breakdown:")
    for cls, count in df_results['classification'].value_counts().items():
        print(f"  {count:4d}  {cls}")

    # ============================================================
    # PER-GROUP SUMMARY
    # ============================================================
    print(f"\n{'Group':10s} {'n':>4s} {'Mean%':>7s} {'Min%':>7s} {'Max%':>7s} {'Artefact':>9s} {'Truncated':>10s} {'Variant':>8s}")
    print("-" * 70)
    for grp in sorted(df_results['group'].unique()):
        sub = df_results[df_results['group'] == grp]
        n_art = sub['classification'].str.contains('ARTEFACT').sum()
        n_trc = sub['classification'].str.contains('TRUNCATED').sum()
        n_var = sub['classification'].str.startswith('VARIANT').sum()
        print(f"{grp:10s} {len(sub):4d} {sub['pct_identity'].mean():7.1f} "
              f"{sub['pct_identity'].min():7.1f} {sub['pct_identity'].max():7.1f} "
              f"{n_art:9d} {n_trc:10d} {n_var:8d}")

    # ============================================================
    # DETAILED SECTIONS
    # ============================================================
    artefacts = df_results[df_results['classification'].str.contains('ARTEFACT')].sort_values('pct_identity', ascending=False)
    truncated = df_results[df_results['classification'].str.contains('TRUNCATED')].sort_values('length_ratio_avg')
    variants = df_results[df_results['classification'].str.startswith('VARIANT')].sort_values('pct_identity')

    if len(artefacts) > 0:
        print(f"\n{'='*80}")
        print(f"LIKELY ARTEFACTS — ≥95% AA identity, ≥80% length ({len(artefacts)} variants):")
        print(f"{'='*80}")
        for _, row in artefacts.iterrows():
            print(f"  {row['group']:8s} {row['variant']:35s} {row['pct_identity']:6.1f}% AA, "
                  f"length ratio: {row['length_ratio_avg']:.0f}% "
                  f"(ref: {row['ref_length_aa']:.0f}AA, grp: {row['group_length_aa_avg']:.0f}AA)")

    if len(truncated) > 0:
        print(f"\n{'='*80}")
        print(f"LIKELY TRUNCATED — <80% length ratio ({len(truncated)} variants):")
        print(f"{'='*80}")
        for _, row in truncated.iterrows():
            print(f"  {row['group']:8s} {row['variant']:35s} {row['pct_identity']:6.1f}% AA, "
                  f"length ratio: {row['length_ratio_avg']:.0f}% "
                  f"(ref: {row['ref_length_aa']:.0f}AA, grp: {row['group_length_aa_avg']:.0f}AA)")

    if len(variants) > 0:
        print(f"\n{'='*80}")
        print(f"VARIANTS — <95% AA identity, ≥80% length ({len(variants)} variants):")
        print(f"{'='*80}")
        for _, row in variants.iterrows():
            print(f"  {row['group']:8s} {row['variant']:35s} {row['pct_identity']:6.1f}% AA, "
                  f"length ratio: {row['length_ratio_avg']:.0f}%")

"""
## 9. Show Failed Lookups

These are variants where the gene couldn't be found in the GFF or FAA file.
May need manual curation of gene names.
"""

if failed:
    df_failed = pd.DataFrame(failed)
    print(f"Failed lookups: {len(df_failed)}")
    print(f"\n{'Group':8s} {'Variant':35s} {'Reason'}")
    print('-' * 90)
    for _, row in df_failed.iterrows():
        print(f"{row['group']:8s} {row['variant']:35s} {row['reason']}")

    # Save failed list
    df_failed.to_csv(os.path.join(OUTPUT_DIR, 'failed_lookups.csv'), index=False)
    print(f"\nSaved to {OUTPUT_DIR}/failed_lookups.csv")
else:
    print("No failed lookups!")

"""
## 10. Save Results
"""

if len(df_results) > 0:
    df_out = df_results.sort_values(['group', 'pct_identity']).reset_index(drop=True)

    out_xlsx = os.path.join(OUTPUT_DIR, 'variant_validation_results.xlsx')
    with pd.ExcelWriter(out_xlsx, engine='openpyxl') as writer:
        # All results (averaged)
        df_out.to_excel(writer, sheet_name='All Results', index=False)

        # Raw per-strain results
        df_results_raw.sort_values(['group', 'variant', 'group_strain']).to_excel(
            writer, sheet_name='Raw Per-Strain', index=False)

        # Separate sheets by classification
        for cls_name, sheet_name in [
            ('ARTEFACT', 'Likely Artefacts'),
            ('TRUNCATED', 'Likely Truncated'),
            ('VARIANT', 'Variants'),
        ]:
            if cls_name == 'VARIANT':
                sub = df_out[df_out['classification'].str.startswith(cls_name)]
            else:
                sub = df_out[df_out['classification'].str.contains(cls_name)]
            if len(sub) > 0:
                sub.to_excel(writer, sheet_name=sheet_name, index=False)

        # Per-group summary
        summary = df_out.groupby('group').agg(
            n_total=('variant', 'count'),
            mean_identity=('pct_identity', 'mean'),
            min_identity=('pct_identity', 'min'),
            max_identity=('pct_identity', 'max'),
            n_artefacts=('classification', lambda x: x.str.contains('ARTEFACT').sum()),
            n_truncated=('classification', lambda x: x.str.contains('TRUNCATED').sum()),
            n_variants=('classification', lambda x: x.str.startswith('VARIANT').sum()),
        ).round(1)
        summary.to_excel(writer, sheet_name='Group Summary')

    # Also save CSV
    df_out.to_csv(os.path.join(OUTPUT_DIR, 'variant_validation_results.csv'), index=False)

    print(f"Results saved to: {out_xlsx}")
    print(f"\nSheets: All Results, Raw Per-Strain, Likely Artefacts, Likely Truncated, Variants, Group Summary")

"""
## 11. Visualization (Optional)
"""

import matplotlib.pyplot as plt
import numpy as np

if len(df_results) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Histogram of all identities
    ax1 = axes[0]
    ax1.hist(df_results['pct_identity'], bins=30, edgecolor='black', alpha=0.7, color='steelblue')
    ax1.axvline(x=95, color='red', linestyle='--', label='95% threshold (artefact?)')
    ax1.axvline(x=80, color='orange', linestyle='--', label='80% threshold')
    ax1.set_xlabel('% Amino Acid Identity', fontsize=12)
    ax1.set_ylabel('Number of Variants', fontsize=12)
    ax1.set_title('Distribution of Variant Divergence', fontsize=14)
    ax1.legend()

    # Per-group boxplot
    ax2 = axes[1]
    groups_sorted = df_results.groupby('group')['pct_identity'].median().sort_values().index
    data_by_group = [df_results[df_results['group'] == g]['pct_identity'].values for g in groups_sorted]
    bp = ax2.boxplot(data_by_group, labels=groups_sorted, vert=True, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
    ax2.axhline(y=95, color='red', linestyle='--', alpha=0.5)
    ax2.axhline(y=80, color='orange', linestyle='--', alpha=0.5)
    ax2.set_ylabel('% AA Identity', fontsize=12)
    ax2.set_title('Variant Divergence by KD Group', fontsize=14)
    ax2.tick_params(axis='x', rotation=90)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'variant_validation_plot.png'), dpi=150, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'variant_validation_plot.svg'), bbox_inches='tight')
    plt.show()
    print("Plot saved!")

"""
# 12. Edit Results
"""

# ============================================================
# 12a. Build lookup: variant name -> validation result
# ============================================================

# df_results should exist from step 8 (averaged across strains)
# df_results_raw should exist from step 8 (per-strain)

validation_lookup = {}  # key: (group, variant_name) -> row from df_results

if len(df_results) > 0:
    for _, row in df_results.iterrows():
        key = (row['group'], row['variant'])
        validation_lookup[key] = row

print(f"Validation lookup built: {len(validation_lookup)} entries")

# Quick check
sample_keys = list(validation_lookup.keys())[:5]
for k in sample_keys:
    r = validation_lookup[k]
    print(f"  {k[0]:8s} {k[1]:35s} -> {r['pct_identity']:.1f}% {r['classification']}")

# ============================================================
# 12b. Process each group sheet and produce corrected output
# ============================================================

xls_feat = pd.ExcelFile(FEATURES_PATH)

# Storage for corrected sheets and summary
corrected_sheets = {}  # group_name -> DataFrame
summary_rows = []

for sheet in xls_feat.sheet_names:
    if sheet.lower() == 'removed':
        continue

    df_feat = pd.read_excel(xls_feat, sheet_name=sheet, header=None)

    # ---- Process PRESENT features (columns A-D) ----
    present_rows = []
    n_orig_present = 0
    n_artefact_present = 0
    n_validated_present = 0

    for i in range(2, len(df_feat)):
        gene = df_feat.iloc[i, 0]
        if not pd.notna(gene) or str(gene).strip() in ['/', '']:
            continue

        gene_str = str(gene).strip()
        func = str(df_feat.iloc[i, 1]).strip() if df_feat.shape[1] > 1 and pd.notna(df_feat.iloc[i, 1]) else ''
        broad = str(df_feat.iloc[i, 2]).strip() if df_feat.shape[1] > 2 and pd.notna(df_feat.iloc[i, 2]) else ''
        loc = str(df_feat.iloc[i, 3]).strip() if df_feat.shape[1] > 3 and pd.notna(df_feat.iloc[i, 3]) else ''

        n_orig_present += 1
        is_variant = bool(re.search(r'\.\d+$', gene_str))

        # Look up validation result
        key = (sheet, gene_str)
        val_result = validation_lookup.get(key, None)

        if val_result is not None:
            pct_id = val_result['pct_identity']
            classification = val_result['classification']
            ref_var = val_result['ref_variant']
            n_strains = val_result.get('n_strains_compared', 1)

            # Check if artefact
            if 'ARTEFACT' in classification:
                # Check multi-copy: if there's a big spread, it's NOT an artefact
                identity_spread = val_result.get('identity_spread', 0)
                if pd.notna(identity_spread) and identity_spread > 10:
                    aa_note = f'{pct_id:.1f}% AA, cov {val_result["coverage_ref"]:.0f}%ref (multi-copy, spread {identity_spread:.0f}%)'
                    status = ' KEEP (multi-copy)'
                    n_validated_present += 1
                else:
                    aa_note = f'{pct_id:.1f}% AA, cov {val_result["coverage_ref"]:.0f}%ref to {ref_var}'
                    status = ' ARTEFACT — REMOVE'
                    n_artefact_present += 1
            else:
                aa_note = f'{pct_id:.1f}% AA, cov {val_result["coverage_ref"]:.0f}%ref/{val_result["coverage_query"]:.0f}%qry to {ref_var}'
                status = classification
                n_validated_present += 1
        elif is_variant:
            aa_note = 'NOT VALIDATED (gene not found in GBK)'
            status = ' UNVALIDATED'
            n_validated_present += 1
        else:
            aa_note = '—'
            status = '— (not a variant)'
            n_validated_present += 1

        present_rows.append({
            'gene': gene_str,
            'function': func,
            'broad_function': broad,
            'location': loc,
            'AA_identity': aa_note,
            'status': status,
        })

    # ---- Process ABSENT features (columns I-L) ----
    absent_rows = []
    n_orig_absent = 0
    n_artefact_absent = 0
    n_validated_absent = 0

    for i in range(2, len(df_feat)):
        if df_feat.shape[1] <= 8:
            break
        gene = df_feat.iloc[i, 8]
        if not pd.notna(gene) or str(gene).strip() in ['/', '']:
            continue

        gene_str = str(gene).strip()
        func = str(df_feat.iloc[i, 9]).strip() if df_feat.shape[1] > 9 and pd.notna(df_feat.iloc[i, 9]) else ''
        broad = str(df_feat.iloc[i, 10]).strip() if df_feat.shape[1] > 10 and pd.notna(df_feat.iloc[i, 10]) else ''
        loc = str(df_feat.iloc[i, 11]).strip() if df_feat.shape[1] > 11 and pd.notna(df_feat.iloc[i, 11]) else ''

        n_orig_absent += 1
        is_variant = bool(re.search(r'\.\d+$', gene_str))

        # If the corresponding present variant was flagged as artefact,
        # this absent .0 should also be removed
        corresponding_artefact = False
        if is_variant:
            base = re.sub(r'\.\d+$', '', gene_str)
            for prow in present_rows:
                if prow['gene'].startswith(base + '.') and 'ARTEFACT' in prow['status']:
                    corresponding_artefact = True
                    break

        if corresponding_artefact:
            status = ' ARTEFACT — REMOVE (matches removed present variant)'
            gene_note = ''
            n_artefact_absent += 1
        else:
            # Check if any copy of the base gene exists in group strain genomes
            # This distinguishes true gene loss from allelic replacement
            base_check = re.sub(r'\.\d+$', '', gene_str) if is_variant else gene_str
            gene_names_check = [n.strip() for n in base_check.split(',')]

            grp_strains_list = GROUP_STRAINS.get(sheet, [])
            found_in_any = False
            for gs in grp_strains_list:
                hits = get_protein_for_gene(gs, gene_names_check, strain_data)
                if hits:
                    found_in_any = True
                    break

            if found_in_any:
                gene_note = 'other variant in genome'
                status = '— (absent, but other gene variant is present in genome)'
            else:
                gene_note = 'gene missing from genome'
                status = '— (absent, gene missing from genome)'
            n_validated_absent += 1

        absent_rows.append({
            'gene': gene_str,
            'function': func,
            'broad_function': broad,
            'location': loc,
            'gene_check': gene_note,
            'status': status,
        })

    # ---- Store ----
    df_present = pd.DataFrame(present_rows)
    df_absent = pd.DataFrame(absent_rows)

    corrected_sheets[sheet] = {
        'present': df_present,
        'absent': df_absent,
    }

    summary_rows.append({
        'group': sheet,
        'orig_present': n_orig_present,
        'artefact_present': n_artefact_present,
        'corrected_present': n_validated_present,
        'orig_absent': n_orig_absent,
        'artefact_absent': n_artefact_absent,
        'corrected_absent': n_validated_absent,
        'orig_total': n_orig_present + n_orig_absent,
        'corrected_total': n_validated_present + n_validated_absent,
        'removed_total': n_artefact_present + n_artefact_absent,
    })

df_summary = pd.DataFrame(summary_rows)

# ---- Print summary ----
print('=' * 90)
print('CORRECTED FEATURES SUMMARY')
print('=' * 90)
print(f"\n{'Group':10s} {'Orig':>6s} {'Removed':>8s} {'Corrected':>10s} {'Orig_P':>7s} {'Art_P':>6s} {'Corr_P':>7s} {'Orig_A':>7s} {'Art_A':>6s} {'Corr_A':>7s}")
print('-' * 90)

total_removed = 0
for _, row in df_summary.iterrows():
    total_removed += row['removed_total']
    marker = ' ' if row['removed_total'] > 0 else ''
    print(f"{row['group']:10s} {row['orig_total']:6d} {row['removed_total']:8d} {row['corrected_total']:10d} "
          f"{row['orig_present']:7d} {row['artefact_present']:6d} {row['corrected_present']:7d} "
          f"{row['orig_absent']:7d} {row['artefact_absent']:6d} {row['corrected_absent']:7d}{marker}")

print(f"\nTotal features removed as artefacts: {total_removed}")
print(f"Groups affected: {(df_summary['removed_total'] > 0).sum()} / {len(df_summary)}")

zeroed = df_summary[df_summary['corrected_total'] == 0]
if len(zeroed) > 0:
    print(f"\n  Groups with ZERO features after correction:")
    for _, row in zeroed.iterrows():
        print(f"    {row['group']} (had {row['orig_total']} features, all removed)")

# ============================================================
# 12c. Write corrected features Excel
# ============================================================

out_corrected = os.path.join(OUTPUT_DIR, 'corrected_features_analysis.xlsx')

with pd.ExcelWriter(out_corrected, engine='openpyxl') as writer:

    df_summary.to_excel(writer, sheet_name='SUMMARY', index=False)

    for group_name, data in corrected_sheets.items():
        df_p = data['present']
        df_a = data['absent']

        max_rows = max(len(df_p), len(df_a)) if len(df_p) > 0 or len(df_a) > 0 else 0

        combined_data = []

        combined_data.append({
            'P_gene': 'PRESENT FEATURES', 'P_function': '', 'P_broad': '',
            'P_location': '', 'P_AA_identity': '', 'P_status': '',
            'gap': '',
            'A_gene': 'ABSENT FEATURES', 'A_function': '', 'A_broad': '',
            'A_location': '', 'A_gene_check': '', 'A_status': '',
        })
        combined_data.append({
            'P_gene': 'gene', 'P_function': 'function', 'P_broad': 'broad function',
            'P_location': 'location', 'P_AA_identity': 'AA% vs ref', 'P_status': 'status',
            'gap': '',
            'A_gene': 'gene', 'A_function': 'function', 'A_broad': 'broad function',
            'A_location': 'location', 'A_gene_check': 'gene in genome?', 'A_status': 'status',
        })

        for j in range(max_rows):
            row = {}
            if j < len(df_p):
                row['P_gene'] = df_p.iloc[j]['gene']
                row['P_function'] = df_p.iloc[j]['function']
                row['P_broad'] = df_p.iloc[j]['broad_function']
                row['P_location'] = df_p.iloc[j]['location']
                row['P_AA_identity'] = df_p.iloc[j]['AA_identity']
                row['P_status'] = df_p.iloc[j]['status']
            else:
                row['P_gene'] = ''
                row['P_function'] = ''
                row['P_broad'] = ''
                row['P_location'] = ''
                row['P_AA_identity'] = ''
                row['P_status'] = ''

            row['gap'] = ''

            if j < len(df_a):
                row['A_gene'] = df_a.iloc[j]['gene']
                row['A_function'] = df_a.iloc[j]['function']
                row['A_broad'] = df_a.iloc[j]['broad_function']
                row['A_location'] = df_a.iloc[j]['location']
                row['A_gene_check'] = df_a.iloc[j].get('gene_check', '')
                row['A_status'] = df_a.iloc[j]['status']
            else:
                row['A_gene'] = ''
                row['A_function'] = ''
                row['A_broad'] = ''
                row['A_location'] = ''
                row['A_gene_check'] = ''
                row['A_status'] = ''

            combined_data.append(row)

        df_combined = pd.DataFrame(combined_data)
        df_combined.to_excel(writer, sheet_name=group_name, index=False, header=False)

print(f"\n Corrected features analysis saved to:")
print(f"   {out_corrected}")
print(f"\nSheets: SUMMARY + {len(corrected_sheets)} group sheets")

# ============================================================
# 13. Quick Overview: What Changed Per Group
# ============================================================

print('=' * 90)
print('DETAILED CHANGES PER GROUP')
print('=' * 90)

for group_name, data in sorted(corrected_sheets.items()):
    df_p = data['present']
    df_a = data['absent']

    if len(df_p) == 0 and len(df_a) == 0:
        continue

    removed_p = df_p[df_p['status'].str.contains('ARTEFACT', na=False)] if len(df_p) > 0 else pd.DataFrame()
    removed_a = df_a[df_a['status'].str.contains('ARTEFACT', na=False)] if len(df_a) > 0 else pd.DataFrame()

    n_removed = len(removed_p) + len(removed_a)
    n_total = len(df_p) + len(df_a)

    if n_removed == 0:
        continue

    print(f"\n--- {group_name} ({n_total} features, {n_removed} removed) ---")

    if len(removed_p) > 0:
        print(f"  Present features REMOVED:")
        for _, row in removed_p.iterrows():
            print(f"     {row['gene']:35s} {row['AA_identity']}")

    if len(removed_a) > 0:
        print(f"  Absent features REMOVED (paired with removed present variant):")
        for _, row in removed_a.iterrows():
            print(f"     {row['gene']}")

print(f"\n\n{'='*90}")
print('VALIDATED VARIANT IDENTITIES BY GROUP')
print('=' * 90)

for group_name, data in sorted(corrected_sheets.items()):
    df_p = data['present']
    if len(df_p) == 0:
        continue

    variants = df_p[
        (df_p['AA_identity'] != '—') &
        (~df_p['status'].str.contains('ARTEFACT', na=False)) &
        (~df_p['status'].str.contains('UNVALIDATED', na=False))
    ]

    if len(variants) == 0:
        continue

    print(f"\n  {group_name} ({len(variants)} validated variants):")
    for _, row in variants.iterrows():
        icon = ''
        if 'TRUNCATED' in row['status']:
            icon = ''
        elif row['status'].startswith('VARIANT'):
            icon = ''
        elif 'multi-copy' in row['status'].lower():
            icon = ''

        print(f"    {icon} {row['gene']:35s} {row['AA_identity']}")
