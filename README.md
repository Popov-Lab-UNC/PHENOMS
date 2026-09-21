# PHENOMS

Python-based Hydrogen-Deuterium Exchange of Molecular Dynamics Simulations

![PHENOMS workflow](./assets/phenoms.png)

**Documentation:** [https://popov-lab-unc.github.io/PHENOMS/](https://popov-lab-unc.github.io/PHENOMS/)

PHENOMS analyzes backbone H-bond networks from MD trajectories for HDX-MS–style interpretation: Rust-accelerated Baker–Hubbard detection, replicate/`ComparisonSet` workflows, occupancy heatmaps, differential protection, and connectivity exports. Backbone N–O is the default; all-bond mode and native traj+topology inputs are optional.

## Install

Install [Rust](https://rustup.rs/) first, then:

```bash
conda env create -f environment.yml
conda activate phenoms
pip install -e .
```

Or with [uv](https://docs.astral.sh/uv/):

```bash
uv venv
source .venv/bin/activate
uv pip install -e .
```

`pip install -e .` (or `uv pip install -e .`) builds the Rust kernel
(`phenoms.phenoms_hbond_rs`) into the same package via
[maturin](https://www.maturin.rs/). Without a Rust toolchain at install time,
PHENOMS falls back to MDTraj for detection.

GROMACS `.tpr` topologies require the optional `MDAnalysis` dependency:
`pip install -e ".[gromacs]"` when building from source (or
`pip install "phenoms[gromacs]"` if you installed from PyPI). Every other
topology format (`.gro`, `.pdb`, `.prmtop`, …) needs no extra — see
[Inputs](https://popov-lab-unc.github.io/PHENOMS/usage.html#inputs) for the
full file-type breakdown by engine.

## Quick start

```python
from phenoms import SimulationSet, ComparisonSet

sim = SimulationSet(
    pdb_files=["rep1.pdb", "rep2.pdb"],
    resid_range=(50, 70),
    sub_frames=100,
)
sim.run()
# Writes CSVs, heatmap plots, and a structure-colored PDB to a fresh timestamped
# dir under default_output_root() by default. Pass output_dir=<path> to control
# where, or output_dir=False for in-memory-only results.

a = SimulationSet(["a1.pdb", "a2.pdb"], sub_frames=100)
b = SimulationSet(["b1.pdb", "b2.pdb"], sub_frames=100)
cmp = ComparisonSet(a, b, label_a="apo", label_b="holo")
a.run(); b.run()
cmp.compare()
# Writes comparison.csv, a difference plot, aligned heatmaps, and a differential
# structure-colored PDB the same way — same output_dir=/False controls apply.
```

CLI (optional): `phenoms prep`, `phenoms run`, `phenoms compare` — see the [docs](https://popov-lab-unc.github.io/PHENOMS/).

Outputs default to `./phenom_outputs/` (or `$PHENOMS_OUTPUT_DIR`).

## Benchmarks

Kernel timing (Docker, 4 CPUs, 3×250 frames, 2,722 atoms): Rust **0.047 s** vs MDTraj **0.803 s** (~17×) vs MDAnalysis **1.126 s** (~24×). Reproduce via [docker/README.md](./docker/README.md).

## Validation data

EGFR (WT vs. del747-749) and PLCg2 MD trajectories used to validate PHENOMS are on Zenodo: [10.5281/zenodo.22879722](https://doi.org/10.5281/zenodo.22879722) (CC0). Get them with `python scripts/fetch_zenodo_data.py` — see [docs/data.rst](./docs/data.rst) for details.

## License

[MIT](./LICENSE)
