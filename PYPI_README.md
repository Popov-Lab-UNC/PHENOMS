# PHENOMS

Python-based Hydrogen-Deuterium Exchange of Molecular Dynamics Simulations

PHENOMS analyzes backbone H-bond networks from MD trajectories for HDX-MS–style
interpretation: Rust-accelerated Baker–Hubbard detection, replicate/`ComparisonSet`
workflows, occupancy heatmaps, differential protection, and connectivity exports.

## Install

```bash
pip install phenoms
```

Reading GROMACS `.tpr` topologies needs one optional extra:

```bash
pip install "phenoms[gromacs]"
```

Every other supported format (`.gro`/`.pdb` topologies, OpenMM `.dcd`, AMBER
`.nc`/`.prmtop`, multi-frame PDBs) works with the plain install above — see
the [Inputs docs](https://popov-lab-unc.github.io/PHENOMS/usage.html#inputs)
for the full file-type breakdown.

## Quick start

```python
from phenoms import SimulationSet, default_output_root

sim = SimulationSet(
    pdb_files=["rep1.pdb", "rep2.pdb"],
    sub_frames=100,
    output_dir=default_output_root() / "my_run",
)
sim.run()
```

A CLI is also available: `phenoms prep`, `phenoms run`, `phenoms compare`.

**Docs:** https://popov-lab-unc.github.io/PHENOMS/
**Source:** https://github.com/Brandon-Cole/PHENOMS

## License

MIT
