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
_PLCG2_PRMTOP = ("plcg2_wtPLCg2.prmtop", "wtPLCg2.prmtop")
# The trajectory exceeded what a single upload to Zenodo's bucket API could
# reliably complete (large uploads were failing with a proxy timeout), so
# it's stored as 15 sequential 500MB parts and reassembled here.
_PLCG2_NC_PARTS = [f"plcg2_wtPLCg2-1.3500ns.nc.part{i:02d}" for i in range(15)]
_PLCG2_NC_LOCAL_NAME = "wtPLCg2-1.3500ns.nc"


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
        prmtop_remote, prmtop_local = _PLCG2_PRMTOP
        _download(prmtop_remote, _PLCG2_DIR / prmtop_local)

        nc_path = _PLCG2_DIR / _PLCG2_NC_LOCAL_NAME
        if nc_path.exists():
            print(f"  skip (already present): {nc_path}")
        else:
            part_paths = []
            for i, remote in enumerate(_PLCG2_NC_PARTS):
                part_path = _PLCG2_DIR / f".{_PLCG2_NC_LOCAL_NAME}.part{i:02d}"
                _download(remote, part_path)
                part_paths.append(part_path)
            print(f"  reassembling {len(part_paths)} parts -> {nc_path}")
            with open(nc_path, "wb") as out:
                for part_path in part_paths:
                    with open(part_path, "rb") as part:
                        while chunk := part.read(1024 * 1024):
                            out.write(chunk)
            for part_path in part_paths:
                part_path.unlink()

    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
