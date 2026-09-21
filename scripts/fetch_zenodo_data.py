#!/usr/bin/env python
"""
Download the PHENOMS validation datasets (EGFR WT/del747-749, PLCg2) from Zenodo.

DOI: 10.5281/zenodo.22879722

Usage:
    python scripts/fetch_zenodo_data.py             # both datasets
    python scripts/fetch_zenodo_data.py --egfr-only
    python scripts/fetch_zenodo_data.py --plcg2-only
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
_RECORD_ID = "22879722"
_BASE_URL = f"https://zenodo.org/records/{_RECORD_ID}/files"

_EGFR_DIR = _REPO_ROOT / "protein_pdb_exports"
_EGFR_FILES = [
    "egfr_3_19_26_wt_rep_1_protein_centered_fit_nopbc.pdb",
    "egfr_3_19_26_wt_rep_2_protein_centered_fit_nopbc.pdb",
    "egfr_3_19_26_wt_rep_3_protein_centered_fit_nopbc.pdb",
    "egfr_3_19_26_del747-749_rep1_protein_centered_fit_nopbc.pdb",
    "egfr_3_19_26_del747-749_rep2_protein_centered_fit_nopbc.pdb",
    "egfr_3_19_26_del747-749_rep3_protein_centered_fit_nopbc.pdb",
]
# Local filenames drop the egfr_/plcg2_ prefix Zenodo needed (its bucket API has no subfolders).
_EGFR_LOCAL_NAMES = [name.removeprefix("egfr_") for name in _EGFR_FILES]

_PLCG2_DIR = _REPO_ROOT / "plcg2_data"
_PLCG2_FILES = ["plcg2_wtPLCg2.prmtop", "plcg2_wtPLCg2-1.3500ns.nc"]
_PLCG2_LOCAL_NAMES = [name.removeprefix("plcg2_") for name in _PLCG2_FILES]


def _download(remote_name: str, local_path: Path) -> None:
    if local_path.exists():
        print(f"  skip (already present): {local_path}")
        return
    local_path.parent.mkdir(parents=True, exist_ok=True)
    url = f"{_BASE_URL}/{remote_name}?download=1"
    print(f"  downloading {remote_name} -> {local_path}")
    urllib.request.urlretrieve(url, local_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--egfr-only", action="store_true")
    group.add_argument("--plcg2-only", action="store_true")
    args = parser.parse_args()

    if not args.plcg2_only:
        print("EGFR (WT / del747-749), ~1GB total:")
        for remote, local in zip(_EGFR_FILES, _EGFR_LOCAL_NAMES):
            _download(remote, _EGFR_DIR / local)

    if not args.egfr_only:
        print("PLCg2, ~7.5GB total:")
        for remote, local in zip(_PLCG2_FILES, _PLCG2_LOCAL_NAMES):
            _download(remote, _PLCG2_DIR / local)

    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
