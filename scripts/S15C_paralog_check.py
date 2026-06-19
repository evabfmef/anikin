"""
# Paralog Check: Multi-Copy Gene Detection

**Author:** Eva Stare

For each group-specific gene identified in the features analysis, this notebook scans
all 67 genomes to determine whether the gene exists in multiple copies (paralogs or
split genes).

"""

import pandas as pd
import re, os, glob
from collections import defaultdict
from Bio import SeqIO
from urllib.parse import unquote
import warnings
warnings.filterwarnings('ignore')
print('Ready!')

BASE_DIR = 'input/'  # UPDATE THIS PATH
FEATURES_PATH = os.path.join(BASE_DIR, 'Features_analysis_67_output.xlsx')
GBK_DIR = os.path.join(BASE_DIR, '67_strains_gbks')
OUTPUT_DIR = os.path.join(BASE_DIR, '03_paralog_check_output_v5')
os.makedirs(OUTPUT_DIR, exist_ok=True)

gbk_files = sorted(
    glob.glob(os.path.join(GBK_DIR, '*.gbk')) +
    glob.glob(os.path.join(GBK_DIR, '*.gb')) +
    glob.glob(os.path.join(GBK_DIR, '*.gbff'))
)
print(f'GBK files: {len(gbk_files)}')

"""
## 1. Parse all GBK files — index gene name → genomic positions per strain
"""

# ============================================================
# STRAIN NAME <-> FILENAME MAPPING
# Needed because non-PS strains have GCF accession filenames
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

# Build reverse mapping: filename_base -> strain_name
FILENAME_TO_STRAIN = {v: k for k, v in STRAIN_TO_FILENAME.items()}


def parse_all_gbk(gbk_files):
    """
    Returns: {strain: {gene_name: [list of (start, end, strand, length_aa, translation)]}}
    Gene names extracted from [gene=XXX] in /product or /gene qualifier.

    Uses FILENAME_TO_STRAIN to convert GCF filenames back to strain names.
    """
    data = {}
    for path in gbk_files:
        bn = os.path.basename(path)
        file_base = re.sub(r'(_edit2?|_edited?)?\.gb[kf]?$', '', bn)

        # Map filename back to strain name (e.g., GCF-009662275.1 -> 73)
        strain = FILENAME_TO_STRAIN.get(file_base, file_base)

        genes = defaultdict(list)
        try:
            for rec in SeqIO.parse(path, 'genbank'):
                for feat in rec.features:
                    if feat.type != 'CDS':
                        continue
                    # Extract gene name
                    names = set()
                    for g in feat.qualifiers.get('gene', []):
                        g = g.strip()
                        if g and not g.startswith(('CUSTOM_', 'BSU')):
                            names.add(g)
                    if not names:
                        for prod in feat.qualifiers.get('product', []):
                            for m in re.finditer(r'\[gene(?:%3D|=)([^\]]+)\]', prod, re.IGNORECASE):
                                gn = unquote(m.group(1)).strip()
                                if gn and not gn.startswith(('BSU', 'CUSTOM_')):
                                    names.add(gn)
                    if not names:
                        continue
                    # Get translation
                    tr = feat.qualifiers.get('translation', [None])[0]
                    if not tr:
                        try:
                            tr = str(feat.location.extract(rec.seq).translate(to_stop=True))
                        except:
                            tr = ''
                    start = int(feat.location.start)
                    end = int(feat.location.end)
                    strand = feat.location.strand
                    length_aa = len(tr) if tr else 0
                    entry = (start, end, strand, length_aa, tr)
                    for gn in names:
                        genes[gn].append(entry)
        except Exception as e:
            print(f'  ERROR {bn}: {e}')
        # Dedup by position within each gene name
        for gn in genes:
            seen = set()
            unique = []
            for ent in genes[gn]:
                key = (ent[0], ent[1], ent[2])  # start, end, strand
                if key not in seen:
                    seen.add(key)
                    unique.append(ent)
            genes[gn] = unique
        data[strain] = dict(genes)
    return data

print('Parsing GBK files...')
strain_genes = parse_all_gbk(gbk_files)
print(f'Parsed {len(strain_genes)} strains')
avg = sum(len(v) for v in strain_genes.values()) / len(strain_genes)
print(f'Average genes per strain: {avg:.0f}')

# Verify that key strains are present with correct names
check_strains = ['NCIB_3610', 'PS-216', 'NRS6085', '73', 'BS16045', 'PS-108']
for s in check_strains:
    if s in strain_genes:
        print(f'   {s}: {len(strain_genes[s])} genes')
    else:
        print(f'    {s} NOT FOUND — check mapping!')

print(f'\nAll strain names: {sorted(strain_genes.keys())}')

"""
## 2. Collect all gene names from Features analysis
"""

xls = pd.ExcelFile(FEATURES_PATH)
varpat = re.compile(r'^(.+)\.(\d+)$')

# Collect all genes mentioned in the Features file
# For each: group, full feature name, base gene names to search, variant suffix (if any), section
all_features = []

for sheet in xls.sheet_names:
    if sheet.lower() == 'removed':
        continue
    df = pd.read_excel(xls, sheet_name=sheet, header=None)
    section = 'unique'
    for idx in range(2, len(df)):
        val = str(df.iloc[idx, 0]).strip() if pd.notna(df.iloc[idx, 0]) else ''
        if val.lower().startswith('absent'):
            section = 'absent'
            continue
        if not val or val in ('/', 'nan') or val.lower() == sheet.lower():
            continue
        # Parse gene names — strip variant suffixes for searching
        parts = [p.strip() for p in val.split(',')]
        search_names = []
        suffix = None
        for p in parts:
            m = varpat.match(p)
            if m:
                search_names.append(m.group(1))
                suffix = int(m.group(2))
            else:
                search_names.append(p)

        func = str(df.iloc[idx, 1])[:80] if pd.notna(df.iloc[idx, 1]) else ''
        all_features.append(dict(
            group=sheet, full_name=val, search_names=search_names,
            variant_suffix=suffix, section=section, function=func
        ))

print(f'Total features in Features analysis: {len(all_features)}')
n_with_var = sum(1 for f in all_features if f['variant_suffix'] is not None)
n_no_var = len(all_features) - n_with_var
print(f'  With variant suffix: {n_with_var}')
print(f'  Without variant suffix (plain genes): {n_no_var}')

"""
## 3. For each feature, check every strain for multiple copies
"""

# KD group assignments for 67 strains (26 groups)
# From the 67-strain variant validation (KD_variant_validation_gbk_v6)
group_strains = {
    'g3':    ['RO-F-3'],
    'g4':    ['RS-D-2'],
    'g5':    ['RO-DD-2'],
    'g7':    ['RO-A-4'],
    'g19':   ['NRS6085'],
    'g25_A': ['PS-13', 'PS-14', 'PS-18', 'PS-30', 'PS-31', 'PS-51', 'PS-65',
              'PS-68', 'PS-96', 'PS-168', 'PS-210', 'PS-216', 'PS-233', 'PS-237',
              'MB8_B1', 'MB8_B7', 'MB9_B1', 'MB9_B6', 'P8_B1', 'P9_B1', 'NCIB_3610'],
    'g29':   ['BS16045'],
    'g30':   ['NRS6105', 'NRS6145'],
    'g31':   ['PS-52', 'PS-53'],
    'g34':   ['NRS6128'],
    'g39':   ['PS-209'],
    'g44':   ['PS-196'],
    'g50':   ['PS-108', 'PS-109', 'PS-119', 'PS-130', 'PS-131'],
    'g52':   ['NRS6160'],
    'g53':   ['PS-20', 'PS-24', 'PS-25'],
    'g54':   ['PS-160'],
    'g55':   ['PS-93', 'PS-95'],
    'g56':   ['PS-217', 'PS-218'],
    'g57':   ['PS-15'],
    'g60_A': ['PS-64', 'PS-194', 'NRS6186', 'NRS6181'],
    'g60_B': ['NRS6132', 'NRS6099', 'NRS6187'],
    'g63_A': ['PS-54', 'PS-55', 'PS-149', 'NRS6103', 'NRS6118', 'NRS6127', 'NRS6153'],
    'g82':   ['RO-FF-1'],
    'g92':   ['KF24'],
    'g95':   ['PS-263'],
    'g96':   ['73'],
}

# Verify all group strains exist in parsed data
all_strain_names = sorted(strain_genes.keys())
print(f'Total parsed strains: {len(all_strain_names)}')

missing_from_gbk = []
total_group_strains = 0
for grp, strains in sorted(group_strains.items()):
    total_group_strains += len(strains)
    for s in strains:
        if s not in strain_genes:
            missing_from_gbk.append((grp, s))

print(f'Total strains across all groups: {total_group_strains}')
if missing_from_gbk:
    print(f'\n  {len(missing_from_gbk)} group strains NOT FOUND in parsed GBK data:')
    for grp, s in missing_from_gbk:
        print(f'    {grp}: {s}')
else:
    print(' All group strains found in parsed GBK data!')

# Strains parsed but not in any group (the 3 extras + any others)
grouped = set()
for strains in group_strains.values():
    grouped.update(strains)
extras = sorted(set(all_strain_names) - grouped)
if extras:
    print(f'\nStrains in GBK folder but not in any group ({len(extras)}): {extras}')
    print('  (These are fine — they will be used for ALL-STRAINS context checks)')

SPLIT_GENE_DISTANCE_BP = 5000

def analyze_copies(all_copies_sorted):
    """Analyze a list of copies (sorted by start) and return distances, strands, flags."""
    lengths = [c[3] for c in all_copies_sorted]
    positions = [(c[0], c[1]) for c in all_copies_sorted]
    strands = [c[2] for c in all_copies_sorted]

    distances = []
    for k in range(len(all_copies_sorted) - 1):
        distances.append(all_copies_sorted[k+1][0] - all_copies_sorted[k][1])

    same_strand = len(set(strands)) == 1
    min_distance = min(distances) if distances else 999999
    has_adjacent = min_distance < SPLIT_GENE_DISTANCE_BP

    # Classification — treat ADJACENT (close + same strand) as potential split too
    if has_adjacent and same_strand:
        if min(lengths)/max(lengths) < 0.9:
            proximity_flag = f'LIKELY SPLIT GENE ({min_distance}bp apart, same strand, unequal lengths)'
        else:
            proximity_flag = f'LIKELY SPLIT GENE ({min_distance}bp apart, same strand, similar lengths)'
    elif has_adjacent and not same_strand:
        proximity_flag = f'Close copies ({min_distance}bp apart, different strands)'
    else:
        proximity_flag = f'Distant paralogs (min {min_distance}bp apart)'

    is_likely_split = has_adjacent and same_strand

    return dict(
        lengths=lengths, positions=positions, strands=strands,
        distances=distances, min_distance=min_distance,
        same_strand=same_strand, proximity_flag=proximity_flag,
        is_likely_split=is_likely_split,
    )


def check_paralogs(all_features, strain_genes, scope='group'):
    """Check for paralogs. scope='group' or 'all'."""
    paralog_rows = []
    features_with_paralogs = []

    for feat in all_features:
        search_names = feat['search_names']
        group = feat['group']

        if scope == 'group':
            strains_to_check = group_strains.get(group, [])
        else:
            strains_to_check = all_strain_names

        if not strains_to_check:
            continue

        n_group = len(group_strains.get(group, []))
        strains_with_multicopy = []

        for strain in sorted(strains_to_check):
            all_copies = []
            seen_pos = set()
            for sn in search_names:
                for entry in strain_genes.get(strain, {}).get(sn, []):
                    pos_key = (entry[0], entry[1], entry[2])
                    if pos_key not in seen_pos:
                        seen_pos.add(pos_key)
                        all_copies.append(entry)

            if len(all_copies) > 1:
                all_copies_sorted = sorted(all_copies, key=lambda c: c[0])
                info = analyze_copies(all_copies_sorted)
                in_group = strain in group_strains.get(group, [])

                strains_with_multicopy.append((strain, len(all_copies_sorted),
                                                info['lengths'], info['positions']))
                paralog_rows.append(dict(
                    group=group, feature=feat['full_name'], section=feat['section'],
                    variant_suffix=f'.{feat["variant_suffix"]}' if feat['variant_suffix'] is not None else 'none',
                    function=feat['function'], strain=strain,
                    in_group='YES' if in_group else 'no',
                    n_copies=len(all_copies_sorted),
                    copy_lengths_aa=str(info['lengths']),
                    copy_positions=str(info['positions']),
                    copy_strands=str(info['strands']),
                    max_length_aa=max(info['lengths']),
                    min_length_aa=min(info['lengths']),
                    length_ratio=round(min(info['lengths'])/max(info['lengths']), 3) if max(info['lengths']) > 0 else 0,
                    likely_truncation='YES' if min(info['lengths'])/max(info['lengths']) < 0.5 else 'no',
                    distances_bp=str(info['distances']),
                    min_distance_bp=info['min_distance'],
                    same_strand='yes' if info['same_strand'] else 'no',
                    proximity_flag=info['proximity_flag'],
                    is_likely_split='YES' if info['is_likely_split'] else 'no',
                ))

        if strains_with_multicopy:
            feature_flags = [r for r in paralog_rows
                             if r['feature'] == feat['full_name'] and r['group'] == group]
            n_split = sum(1 for r in feature_flags if r['is_likely_split'] == 'YES')
            n_total_par = len(strains_with_multicopy)
            in_grp = sum(1 for s,_,_,_ in strains_with_multicopy if s in group_strains.get(group, []))
            out_grp = n_total_par - in_grp

            features_with_paralogs.append(dict(
                group=group, feature=feat['full_name'], section=feat['section'],
                variant_suffix=f'.{feat["variant_suffix"]}' if feat['variant_suffix'] is not None else 'none',
                function=feat['function'],
                group_size=n_group,
                n_in_group=in_grp, n_outside_group=out_grp,
                n_total=n_total_par,
                n_likely_split=n_split,
                strains=', '.join(f'{s}({n}){"*" if s in group_strains.get(group,[]) else ""}'
                                  for s,n,_,_ in strains_with_multicopy)
            ))

    return paralog_rows, features_with_paralogs


# ============================================================
# Run BOTH scopes
# ============================================================
print("\nRunning GROUP-ONLY check...")
group_rows, group_features = check_paralogs(all_features, strain_genes, scope='group')
print(f"  Features with paralogs in group: {len(group_features)}")
print(f"  Detail rows: {len(group_rows)}")

print("\nRunning ALL-STRAINS check...")
all_rows, all_features_par = check_paralogs(all_features, strain_genes, scope='all')
print(f"  Features with paralogs anywhere: {len(all_features_par)}")
print(f"  Detail rows: {len(all_rows)}")

# ============================================================
# Cross-reference: for each split in GROUP, is it also split outside group?
# ============================================================
print(f'\n{"="*70}')
print(f'SPLIT GENE UNIQUENESS CHECK')
print(f'{"="*70}')
print(f'For each gene with likely splits in GROUP strains, check if the same')
print(f'split pattern occurs in strains OUTSIDE the group.')
print()

# Get group-level split features
group_split_features = set()
for r in group_rows:
    if r['is_likely_split'] == 'YES':
        group_split_features.add((r['group'], r['feature']))

for (grp, feat) in sorted(group_split_features):
    all_entries = [r for r in all_rows if r['group'] == grp and r['feature'] == feat]
    in_group_splits = [r for r in all_entries if r['in_group'] == 'YES' and r['is_likely_split'] == 'YES']
    outside_splits = [r for r in all_entries if r['in_group'] == 'no' and r['is_likely_split'] == 'YES']
    outside_nonsplit = [r for r in all_entries if r['in_group'] == 'no' and r['is_likely_split'] == 'no']

    n_grp_split = len(in_group_splits)
    n_out_split = len(outside_splits)
    n_out_nonsplit = len(outside_nonsplit)
    n_group_total = len(group_strains.get(grp, []))

    if n_out_split == 0 and n_out_nonsplit == 0:
        uniqueness = "UNIQUE to group (no copies outside)"
    elif n_out_split == 0 and n_out_nonsplit > 0:
        uniqueness = f"SPLIT UNIQUE to group ({n_out_nonsplit} strains outside have intact distant paralogs)"
    elif n_out_split > 0:
        uniqueness = f" SPLIT ALSO OUTSIDE group ({n_out_split} outside strains also split) — CHECK MANUALLY"
        out_strains = [r['strain'] for r in outside_splits]
        uniqueness += f" [{', '.join(out_strains[:5])}{'...' if len(out_strains) > 5 else ''}]"
    else:
        uniqueness = "?"

    print(f"  {grp:8s} {feat:35s} grp_split={n_grp_split}/{n_group_total}  "
          f"outside_split={n_out_split}  outside_other={n_out_nonsplit}  → {uniqueness}")

# ============================================================
# Print GROUP summary
# ============================================================
print(f'\n{"="*70}')
print(f'GROUP RESULTS: Paralogs within group members')
print(f'{"="*70}')

if group_features:
    print(f'{"Group":8s} {"Feature":40s} {"Sec":8s} {"Suf":7s} {"Grp":>4s} {"#Par":>5s} {"#Spl":>5s}  Strains')
    print('-' * 130)
    for f in sorted(group_features, key=lambda x: (-x['n_likely_split'], -x['n_in_group'], x['group'])):
        spl = ' SPLIT!' if f['n_likely_split'] > 0 else ''
        allm = ' [ALL]' if f['n_in_group'] == f['group_size'] else ''
        print(f'{f["group"]:8s} {f["feature"]:40s} {f["section"]:8s} {f["variant_suffix"]:7s} '
              f'{f["group_size"]:>4d} {f["n_in_group"]:>5d} {f["n_likely_split"]:>5d}  '
              f'{f["strains"][:50]}{spl}{allm}')

# ============================================================
# Print ALL summary (context)
# ============================================================
print(f'\n{"="*70}')
print(f'ALL STRAINS CONTEXT')
print(f'{"="*70}')
only_outside = [f for f in all_features_par if f['n_in_group'] == 0]
both = [f for f in all_features_par if f['n_in_group'] > 0 and f['n_outside_group'] > 0]
only_inside = [f for f in all_features_par if f['n_outside_group'] == 0]
print(f'  Paralogs only in group strains: {len(only_inside)}')
print(f'  Paralogs in group AND other strains: {len(both)}')
print(f'  Paralogs only outside group: {len(only_outside)}')

"""
## 4. Save results
"""

import pandas as pd

out_file = os.path.join(OUTPUT_DIR, '67_paralog_check_results.xlsx')

with pd.ExcelWriter(out_file, engine='openpyxl') as w:
    # ---- GROUP SCOPE ----
    if group_features:
        df = pd.DataFrame(group_features)
        df = df.sort_values(['n_likely_split', 'n_in_group'], ascending=[False, False])
        df.to_excel(w, sheet_name='GROUP - features summary', index=False)

    if group_rows:
        df = pd.DataFrame(group_rows)
        df.sort_values(['group', 'feature', 'strain']).to_excel(
            w, sheet_name='GROUP - all details', index=False)

        df_split = df[df['is_likely_split'] == 'YES']
        if len(df_split):
            df_split.to_excel(w, sheet_name='GROUP - likely splits', index=False)

    # ---- ALL STRAINS SCOPE ----
    if all_features_par:
        df = pd.DataFrame(all_features_par)
        df = df.sort_values(['n_likely_split', 'n_total'], ascending=[False, False])
        df.to_excel(w, sheet_name='ALL - features summary', index=False)

    if all_rows:
        df = pd.DataFrame(all_rows)
        df.sort_values(['group', 'feature', 'strain']).to_excel(
            w, sheet_name='ALL - all details', index=False)

        df_split = df[df['is_likely_split'] == 'YES']
        if len(df_split):
            df_split.to_excel(w, sheet_name='ALL - likely splits', index=False)

    # ---- SPLIT UNIQUENESS REPORT ----
    uniqueness_rows = []
    group_split_features_set = set()
    for r in group_rows:
        if r['is_likely_split'] == 'YES':
            group_split_features_set.add((r['group'], r['feature']))

    for (grp, feat) in sorted(group_split_features_set):
        all_entries = [r for r in all_rows if r['group'] == grp and r['feature'] == feat]
        in_grp_split = [r for r in all_entries if r['in_group'] == 'YES' and r['is_likely_split'] == 'YES']
        out_split = [r for r in all_entries if r['in_group'] == 'no' and r['is_likely_split'] == 'YES']
        out_nonsplit = [r for r in all_entries if r['in_group'] == 'no' and r['is_likely_split'] == 'no']
        n_group = len(group_strains.get(grp, []))

        if len(out_split) == 0 and len(out_nonsplit) == 0:
            status = "UNIQUE to group (no copies outside)"
        elif len(out_split) == 0:
            status = f"Split unique to group ({len(out_nonsplit)} outside have distant paralogs only)"
        else:
            out_strains = ', '.join(r['strain'] for r in out_split[:10])
            status = f"SPLIT ALSO OUTSIDE ({len(out_split)} strains) — CHECK: {out_strains}"

        uniqueness_rows.append(dict(
            group=grp, feature=feat,
            n_group_strains_split=len(in_grp_split), group_size=n_group,
            n_outside_split=len(out_split),
            n_outside_distant=len(out_nonsplit),
            uniqueness=status,
            outside_split_strains=', '.join(r['strain'] for r in out_split),
        ))

    if uniqueness_rows:
        df_uniq = pd.DataFrame(uniqueness_rows)
        df_uniq.to_excel(w, sheet_name='SPLIT uniqueness check', index=False)

    # ---- Clean features ----
    par_keys = set((f['group'], f['feature']) for f in all_features_par)
    clean = [f for f in all_features if (f['group'], f['full_name']) not in par_keys]
    if clean:
        pd.DataFrame(clean).to_excel(w, sheet_name='Clean (no paralogs)', index=False)

print(f'\nSaved: {out_file}')
print(f'Sheets:')
print(f'  GROUP - features summary: {len(group_features)}')
print(f'  GROUP - all details: {len(group_rows)}')
n_gs = sum(1 for r in group_rows if r.get("is_likely_split") == "YES")
print(f'  GROUP - likely splits: {n_gs}')
print(f'  ALL - features summary: {len(all_features_par)}')
print(f'  ALL - all details: {len(all_rows)}')
n_as = sum(1 for r in all_rows if r.get("is_likely_split") == "YES")
print(f'  ALL - likely splits: {n_as}')
print(f'  SPLIT uniqueness check: {len(uniqueness_rows) if uniqueness_rows else 0}')
par_keys = set((f["group"], f["feature"]) for f in all_features_par)
n_clean = len([f for f in all_features if (f["group"], f["full_name"]) not in par_keys])
print(f'  Clean (no paralogs): {n_clean}')
