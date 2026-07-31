# Environment

This directory documents the software environment used for analysis and figure
generation in the anikin project (Stare et al., 2026).

## Contents

- `requirements.txt` — Python packages for the main environment (S1–S9, S10A, S11–S18)
- `requirements_GEM.txt` — Python packages for the GEM environment (S10B, S10C)
- `environment.yml` — conda specification for the main environment
- `R_packages.txt` — R packages and versions used in the analysis and visualization workflows

## Two Python environments

Two environments are required. They cannot be combined.

| File | Scripts | Python | cobra |
|---|---|---|---|
| `requirements.txt` / `environment.yml` | S1–S9, S10A, S11–S18 | 3.12 | 0.29.1 |
| `requirements_GEM.txt` | S10B, S10C | 3.11 | 0.26.0 |

### Why two

S10B (draft strain-specific model generation) requires cobra 0.26.0. Newer
releases were tested and did not reproduce this step correctly, so the pin is
required rather than conservative.

cobra 0.26.0 in turn pins `pandas~=1.0` and `numpy~=1.13`. These cannot be
satisfied alongside the pandas used by every other script: installing both
sets together fails, with pip attempting to build pandas 1.5.3 from source under
Python 3.12.

S10A (bidirectional BLASTp and homology matrix generation) runs under the main
environment and is unaffected.

### Main environment

Conda (recommended):

```bash
conda env create -f environment/environment.yml
conda activate anikin
```

Or pip:

```bash
python3.12 -m venv anikin
source anikin/bin/activate
pip install -r environment/requirements.txt
```

### GEM environment (S10B, S10C)

```bash
conda create -n anikin_gem python=3.11
conda activate anikin_gem
pip install -r environment/requirements_GEM.txt
```

## Reproducibility

Scripts that use random number generation set a fixed seed (42). Seeded results
are reproducible within a given scikit-learn version but are not guaranteed to
be identical across versions, so the scikit-learn pin (1.6.1) is required rather
than advisory.

Beyond the two pinned cases above (cobra 0.26.0 for S10B–S10C, scikit-learn
1.6.1 for S16), the package lists were compiled from the environment used during
the study and are provided for transparency.

## Non-Python dependencies

Bioinformatics tools (fastANI, Roary, PHASTER, antiSMASH, Prokka, MUSCLE,
BLAST+, Canu, NCBI Datasets CLI, Dorado) are not installed by these files. See
`docs/software_versions.md` for the versions used and `R_packages.txt` for the R
environment.