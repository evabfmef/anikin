#!/usr/bin/env python3
"""
S6_genomovar_concordance.py - Concordance of ANI genomovars with the core-genome phylogeny

Tests whether each ANI-defined genomovar (>=99.5% ANI; Supplementary Table S8)
forms an exclusive clade in the core-genome maximum-likelihood phylogeny of the
313-strain panel. Corresponds to Supplementary Methods 4.

Method
------
For every edge of the tree, the two sets of strains on either side of the
corresponding split are recorded. A multi-strain genomovar is concordant if its
membership exactly matches one of these sets, i.e. it contains all members of
that genomovar and no strains from any other genomovar. Because both sides of
every split are considered, the criterion is independent of tree rooting.
Concordance is the percentage of multi-strain genomovars meeting it. Singleton
genomovars are reported but excluded from the percentage, because every single
tip satisfies the criterion trivially. Bootstrap support is read from the input
tree for the matching edges.

Each split is identified in two independent ways (descendant tip sets, and
connected components after removing each edge from an undirected graph); the
script stops if the two disagree. No tree topology, group membership or support
value is inferred or modified.

For discordant genomovars, the closest split is reported together with the
strains missing from it and the strains from other genomovars it contains.

Inputs
------
  --tree       core-genome Newick tree with bootstrap support (313 strains)
  --workbook   Supplementary Tables workbook; sheet "8" holds one genomovar per
               column, with group labels in row 3 and members below
  --aliases    optional JSON dictionary mapping table labels to tree labels

Outputs (--outdir)
------------------
  genomovar_results.csv   per-genomovar status, support values and, for
                          discordant groups, the nearest split
  strain_mapping.csv      table label -> tree label for every strain
  summary.json            counts, concordance percentage, software versions
                          and SHA-256 checksums of all inputs and this script

Usage
-----
  python S6_genomovar_concordance.py --tree core313.nwk \
      --workbook "Supplementary Tables isme.xlsx" --outdir output/S6

Environment
-----------
Python >= 3.10 with openpyxl.
"""

import argparse
import csv
import hashlib
import json
import platform
import re
from collections import Counter
from dataclasses import dataclass, field
from importlib.metadata import version
from pathlib import Path


# ---------------------------------------------------------------------------
# 1. Newick parsing
# ---------------------------------------------------------------------------
@dataclass(eq=False)
class Node:
    """Tree node; `tips` holds the names of all descendant tips."""
    children: list = field(default_factory=list)
    name: str | None = None
    support: float | None = None
    length: float | None = None
    parent: object = None
    tips: frozenset = field(default_factory=frozenset)


def parse_newick(text):
    """Read one Newick tree, including numeric [support] comments after lengths.

    Supports numeric internal-node labels (bootstrap support) and quoted tip
    names. Non-numeric comments are ignored; extended Newick annotations are
    not decoded.
    """
    pattern = r"\[[^\]]*\]|'(?:[^']|'')*'|[(),:;]|[^\s(),:;\[\]']+"
    tokens = re.findall(pattern, text)
    pos = 0
    nodes = []

    def numeric_support(node, token):
        try:
            value = float(token)
        except ValueError:
            return
        if node.support is not None and node.support != value:
            raise ValueError("Conflicting support annotations on one edge.")
        node.support = value

    def subtree():
        nonlocal pos
        node = Node()
        nodes.append(node)
        if tokens[pos] == "(":
            pos += 1
            while True:
                child = subtree()
                child.parent = node
                node.children.append(child)
                if tokens[pos] == ",":
                    pos += 1
                    continue
                if tokens[pos] != ")":
                    raise ValueError("Expected closing parenthesis in Newick tree.")
                pos += 1
                break
            # Internal-node label (bootstrap support)
            if tokens[pos] not in (":", ",", ")", ";") and not tokens[pos].startswith("["):
                numeric_support(node, tokens[pos])
                pos += 1
        else:
            token = tokens[pos]
            if token in (")", ",", ":", ";") or token.startswith("["):
                raise ValueError("Missing terminal label in Newick tree.")
            node.name = token[1:-1].replace("''", "'") if token.startswith("'") else token
            pos += 1
        # Branch length and bracketed support, in either order
        while pos < len(tokens):
            token = tokens[pos]
            if token == ":":
                pos += 1
                node.length = float(tokens[pos])
                pos += 1
            elif token.startswith("["):
                numeric_support(node, token[1:-1])
                pos += 1
            else:
                break
        return node

    # Skip leading tree annotations such as [&U]
    while pos < len(tokens) and tokens[pos].startswith("["):
        pos += 1
    try:
        root = subtree()
    except IndexError as exc:
        raise ValueError("Incomplete Newick tree.") from exc
    if tokens[pos:] != [";"]:
        raise ValueError("Expected exactly one semicolon-terminated Newick tree.")

    labels = [n.name for n in nodes if not n.children]
    if len(labels) != len(set(labels)):
        raise ValueError("Duplicate tip labels in tree.")

    def descendants(node):
        node.tips = (frozenset([node.name]) if not node.children else
                     frozenset().union(*(descendants(c) for c in node.children)))
        return node.tips

    descendants(root)
    return root, nodes


# ---------------------------------------------------------------------------
# 2. Genomovar assignments (Supplementary Table S8)
# ---------------------------------------------------------------------------
def read_groups(workbook_path, sheet, header_row, aliases):
    """Read one genomovar per column, with members listed below the header row."""
    from openpyxl import load_workbook

    workbook = load_workbook(workbook_path, read_only=True, data_only=True)
    try:
        rows = list(workbook[sheet].values)
    finally:
        workbook.close()

    headers = rows[header_row - 1]
    groups = {}
    mappings = []
    all_labels = []
    for column, heading in enumerate(headers):
        if heading is None:
            continue
        group = str(heading).strip()
        if group in groups:
            raise ValueError(f"Duplicate group heading: {group}")
        members = []
        for row_number, row in enumerate(rows[header_row:], start=header_row + 1):
            if column >= len(row) or row[column] is None:
                continue
            original = str(row[column]).strip()
            if not original:
                continue
            label = aliases.get(original, original)
            members.append(label)
            all_labels.append(label)
            mappings.append({"genomovar": group, "table_label": original,
                             "tree_label": label, "row": row_number,
                             "column": column + 1})
        if not members:
            raise ValueError(f"Empty group: {group}")
        groups[group] = frozenset(members)

    duplicates = [s for s, count in Counter(all_labels).items() if count > 1]
    if duplicates:
        raise ValueError(f"Duplicate strain assignments after aliases: {duplicates}")
    if not groups:
        raise ValueError("No genomovar assignments found.")
    return groups, mappings


# ---------------------------------------------------------------------------
# 3. Split-based concordance
# ---------------------------------------------------------------------------
def analyse(root, nodes, groups):
    """Classify each genomovar as singleton, concordant or discordant."""
    universe = root.tips
    assigned = frozenset().union(*groups.values())
    if assigned != universe:
        raise ValueError(f"Strain labels differ. Tree-only: {sorted(universe - assigned)}; "
                         f"table-only: {sorted(assigned - universe)}")

    # Both sides of every split (root-independent)
    edges = [(i, n) for i, n in enumerate(nodes) if n.parent is not None]
    sides = [(i, n, side) for i, n in edges for side in (n.tips, universe - n.tips)]

    # Independent check: remove each edge from an undirected graph and collect
    # the tips of the resulting component.
    ids = {n: i for i, n in enumerate(nodes)}
    adjacency = {i: set() for i in range(len(nodes))}
    for i, n in edges:
        j = ids[n.parent]
        adjacency[i].add(j)
        adjacency[j].add(i)
    graph_sides = set()
    for i, n in edges:
        parent_id = ids[n.parent]
        seen, stack = set(), [i]
        while stack:
            current = stack.pop()
            if current in seen:
                continue
            seen.add(current)
            stack.extend(j for j in adjacency[current] if j not in seen and
                         {current, j} != {i, parent_id})
        side = frozenset(nodes[j].name for j in seen if not nodes[j].children)
        graph_sides.update((side, universe - side))

    results = []
    for group, members in groups.items():
        matches = [(i, n) for i, n in edges if members in (n.tips, universe - n.tips)]
        if bool(matches) != (members in graph_sides):
            raise RuntimeError("Independent split verification failed.")
        row = {"genomovar": group, "size": len(members),
               "members": ";".join(sorted(members)),
               "status": "singleton" if len(members) == 1 else
                         ("concordant" if matches else "discordant"),
               "support_values": ";".join(str(n.support) for _, n in matches
                                          if n.support is not None)
                                 if len(members) > 1 else "",
               "matching_edges": ";".join(str(i) for i, _ in matches)
                                 if len(members) > 1 else "",
               "nearest_split_difference": "", "missing_from_nearest": "",
               "outsiders_in_nearest": ""}
        # For discordant groups, report the closest split (smallest symmetric difference)
        if row["status"] == "discordant":
            i, n, side = min(sides, key=lambda x: (len(members ^ x[2]), len(x[2]), x[0]))
            row.update(nearest_split_difference=len(members ^ side),
                       missing_from_nearest=";".join(sorted(members - side)),
                       outsiders_in_nearest=";".join(sorted(side - members)))
        results.append(row)

    multi = sum(r["size"] > 1 for r in results)
    concordant = sum(r["status"] == "concordant" for r in results)
    summary = {"tree_tips": len(universe), "genomovars": len(groups),
               "singletons": sum(r["status"] == "singleton" for r in results),
               "multi_strain_genomovars": multi, "concordant": concordant,
               "discordant": [r["genomovar"] for r in results if r["status"] == "discordant"],
               "concordance_percent": 100 * concordant / multi if multi else None,
               "independent_graph_check_passed": True}
    return summary, results


# ---------------------------------------------------------------------------
# 4. Output
# ---------------------------------------------------------------------------
def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Concordance of ANI genomovars with the core-genome phylogeny (S6).")
    parser.add_argument("--tree", type=Path, default=Path("core313.nwk"),
                        help="Core-genome Newick tree with bootstrap support.")
    parser.add_argument("--workbook", type=Path, default=Path("Supplementary Tables isme.xlsx"),
                        help="Supplementary Tables workbook.")
    parser.add_argument("--sheet", default="8", help="Sheet with genomovar assignments.")
    parser.add_argument("--header-row", type=int, default=3,
                        help="Excel row holding the genomovar labels.")
    parser.add_argument("--aliases", type=Path, default=None,
                        help="JSON dictionary mapping table labels to tree labels.")
    parser.add_argument("--outdir", type=Path, default=Path("output/S6"),
                        help="Output directory.")
    args = parser.parse_args()

    aliases = json.loads(args.aliases.read_text(encoding="utf-8-sig")) if args.aliases else {}
    root, nodes = parse_newick(args.tree.read_text(encoding="utf-8-sig"))
    groups, mappings = read_groups(args.workbook, args.sheet, args.header_row, aliases)
    summary, results = analyse(root, nodes, groups)

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "genomovar_results.csv", results)
    write_csv(args.outdir / "strain_mapping.csv", mappings)

    # Provenance: software versions and checksums of every input and this script
    used_paths = [args.tree, args.workbook, Path(__file__)]
    if args.aliases:
        used_paths.append(args.aliases)
    summary["provenance"] = {
        "python_version": platform.python_version(),
        "openpyxl_version": version("openpyxl"),
        "sheet": args.sheet, "header_row": args.header_row, "aliases": aliases,
        "input_and_script_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in used_paths}}
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"\nSaved to: {args.outdir}")


if __name__ == "__main__":
    main()