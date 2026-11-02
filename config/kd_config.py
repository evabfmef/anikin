"""
kd_config.py - Shared configuration loader for the S17 series (S17A-S17D)
---------------------------------------------------------------------
Author: Eva Stare

Loads config_S16_S19.yml and the metadata TSVs it points to (KD group
memberships, representative strains, genome accessions), so that thresholds,
reference-strain choices and strain-level data are defined once and are
reviewable without reading script internals.

Usage in a script:

    from kd_config import load_config
    cfg = load_config()                    # reads config_S16_S19.yml
    cfg = load_config(analysis="39")       # override panel
    cfg = load_config("my_config.yml")     # alternative config file

Paths containing "{analysis}" are substituted automatically.

Command-line overrides supported by the scripts that use this module:
    --config=FILE      configuration file  (default: config_S16_S19.yml, looked up
                                           in the working directory, then config/)
    --analysis=39|67   strain panel        (overrides config)
"""

import os
import sys

import pandas as pd
import yaml

DEFAULT_CONFIG = "config_S16_S19.yml"

# Searched in order when no --config is given, so scripts work whether the
# config sits in the working directory or in the repository's config/ folder.
DEFAULT_CONFIG_LOCATIONS = [DEFAULT_CONFIG, os.path.join("config", DEFAULT_CONFIG)]


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
    config_file = config_file or cli.get("config")
    if not config_file:
        config_file = next((p for p in DEFAULT_CONFIG_LOCATIONS
                            if os.path.exists(p)), DEFAULT_CONFIG)

    if not os.path.exists(config_file):
        raise FileNotFoundError(
            f"Configuration file not found: {config_file}\n"
            f"  Expected {DEFAULT_CONFIG} in the working directory or in "
            f"config/, or pass --config=FILE."
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


def load_strain_aliases(cfg):
    """
    Return {strain_name: genome_file_base} from metadata/accessions.tsv.

    Public genomes are stored under their assembly accession (e.g.
    GCF-009662275.1 for strain 73); in-house genomes under their own name, so
    the result is safe to use as a total lookup.
    """
    path = _require(cfg["metadata"]["accessions"], "Accession table")
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = {"strain_name", "source_filename"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{path} is missing required column(s): {sorted(missing)}. "
            f"Expected at least: strain_name, source_filename."
        )
    names = df["source_filename"].fillna(df["strain_name"])
    return dict(zip(df["strain_name"], names))


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
