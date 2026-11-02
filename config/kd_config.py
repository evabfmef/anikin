"""
kd_config.py - Shared configuration loader for S14 and the S15 series
---------------------------------------------------------------------
Author: Eva Stare

Loads config_S14_S17.yml and the metadata TSVs it points to, so that
thresholds, reference-strain choices and strain-level data are defined once
and are reviewable without reading script internals.

Usage in a script:

    from kd_config import load_config
    cfg = load_config()                    # reads config_S14_S17.yml
    cfg = load_config(analysis="39")       # override panel
    cfg = load_config("my_config.yml")     # alternative config file

Paths containing "{analysis}" are substituted automatically.

Command-line overrides supported by the scripts that use this module:
    --config=FILE      configuration file  (default: config_S14_S17.yml)
    --analysis=39|67   strain panel        (overrides config)
"""

import os
import sys

import pandas as pd
import yaml

DEFAULT_CONFIG = "config_S14_S17.yml"


def parse_cli(argv=None):
    """Parse --key=value arguments. Returns a dict; ignores unknown keys."""
    argv = sys.argv[1:] if argv is None else argv
    out = {}
    for a in argv:
        if a in ("--help", "-h"):
            out["help"] = True
            continue
        if a.startswith("--") and "=" in a:
            k, v = a[2:].split("=", 1)
            out[k.replace("-", "_")] = v
    return out


def _subst(obj, analysis):
    """Recursively substitute {analysis} in string values."""
    if isinstance(obj, str):
        return obj.replace("{analysis}", str(analysis))
    if isinstance(obj, dict):
        return {k: _subst(v, analysis) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_subst(v, analysis) for v in obj]
    return obj


def load_config(config_file=None, analysis=None, argv=None):
    """
    Load the configuration file, apply command-line overrides, and substitute
    {analysis} placeholders throughout.

    Precedence: explicit argument > command line > config file default.
    """
    cli = parse_cli(argv)
    config_file = config_file or cli.get("config") or DEFAULT_CONFIG

    if not os.path.exists(config_file):
        raise FileNotFoundError(
            f"Configuration file not found: {config_file}\n"
            f"  Expected {DEFAULT_CONFIG} in the working directory, "
            f"or pass --config=FILE."
        )

    with open(config_file) as fh:
        cfg = yaml.safe_load(fh)

    analysis = str(analysis or cli.get("analysis") or cfg.get("analysis"))
    if analysis not in ("39", "67"):
        raise ValueError(
            f"analysis must be '39' or '67', got {analysis!r}. "
            f"Set `analysis:` in {config_file} or pass --analysis=67."
        )

    cfg = _subst(cfg, analysis)
    cfg["analysis"] = analysis
    cfg["_config_file"] = config_file
    return cfg


# ---------------------------------------------------------------------------
# Metadata loaders
# ---------------------------------------------------------------------------

def _require(path, what):
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{what} not found: {path}\n"
            f"  Check the `metadata:` block of the configuration file."
        )
    return path


def load_kd_groups(cfg, representatives=False):
    """
    Return {kd_group: [strain, ...]} from the metadata TSV.

    Set representatives=True to load the representatives table instead of the
    full membership table.
    """
    key = "kd_groups_representatives" if representatives else "kd_groups"
    path = _require(cfg["metadata"][key], "KD group table")
    df = pd.read_csv(path, sep="\t")
    for col in ("kd_group", "strain"):
        if col not in df.columns:
            raise ValueError(
                f"{path} is missing required column '{col}'. "
                f"Expected columns: kd_group, strain."
            )
    return df.groupby("kd_group")["strain"].apply(list).to_dict()


def load_strain_aliases(cfg, alias_type="filename_or_assembly"):
    """
    Return {canonical_strain: alias_value} for the requested alias type.

    Strains without an alias of that type map to themselves, so the result is
    safe to use as a total lookup.
    """
    path = _require(cfg["metadata"]["strain_aliases"], "Strain alias table")
    df = pd.read_csv(path, sep="\t")
    required = {"canonical_strain", "alias_type", "alias_value"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{path} is missing required column(s): {sorted(missing)}. "
            f"Expected: canonical_strain, alias_type, alias_value."
        )
    sub = df[df["alias_type"] == alias_type]
    mapping = sub.set_index("canonical_strain")["alias_value"].to_dict()
    return {s: mapping.get(s, s) for s in df["canonical_strain"].unique()}


def load_pinned_exclusions(cfg):
    """
    Return {kd_group: [strain, ...]} of a-priori-fixed vertex-cover solutions.

    Returns an empty dict if pinning is disabled in the config, or if the file
    is absent (pinning is optional: without it, S14 takes the first minimum
    cover in sorted order, which is still deterministic).
    """
    if not cfg.get("s14", {}).get("use_pinned_exclusions", True):
        return {}

    path = cfg["metadata"].get("pinned_exclusions")
    if not path or not os.path.exists(path):
        return {}

    df = pd.read_csv(path, sep="\t")
    for col in ("kd_group", "strain"):
        if col not in df.columns:
            raise ValueError(
                f"{path} is missing required column '{col}'. "
                f"Expected columns: kd_group, strain[, rationale]."
            )
    return df.groupby("kd_group")["strain"].apply(list).to_dict()


def load_manual_curation(cfg):
    """Return the manual curation table as a DataFrame, or None if absent."""
    path = cfg["metadata"].get("manual_curation")
    if not path or not os.path.exists(path):
        return None
    return pd.read_csv(path, sep="\t")


def describe(cfg):
    """One-line-per-item summary printed at script start, for the run log."""
    lines = [
        f"Configuration:      {cfg['_config_file']}",
        f"Analysis panel:     {cfg['analysis']}",
    ]
    return "\n".join(lines)
