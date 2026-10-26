import os
from pathlib import Path

import mdtraj as md
import numpy as np
import pytest

from phenoms.io import load_trajectory
from phenoms.preprocess import (
    detect_replicate_input,
    discover_replicate_dirs,
    normalize_replicate_to_pdb,
)

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
MINIMAL_PDB = os.path.join(FIXTURES_DIR, "minimal.pdb")


def test_detect_gromacs_requires_topology(tmp_path):
    rep = tmp_path / "rep1"
    rep.mkdir()
    (rep / "traj.xtc").write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="GROMACS replicate"):
        detect_replicate_input(rep)


def test_detect_openmm_pair(tmp_path):
    rep = tmp_path / "rep_openmm"
    rep.mkdir()
    (rep / "traj.dcd").write_text("", encoding="utf-8")
    (rep / "top.pdb").write_text("MODEL        1\nENDMDL\nEND\n", encoding="utf-8")
    out = detect_replicate_input(rep)
    assert out.engine == "openmm"
    assert out.trajectory_path.name == "traj.dcd"
    assert out.topology_path.name == "top.pdb"


def test_detect_amber_pair(tmp_path):
    rep = tmp_path / "rep_amber"
    rep.mkdir()
    (rep / "traj.nc").write_text("", encoding="utf-8")
    (rep / "sys.prmtop").write_text("", encoding="utf-8")
    out = detect_replicate_input(rep)
    assert out.engine == "amber"
    assert out.trajectory_path.name == "traj.nc"
    assert out.topology_path.name == "sys.prmtop"


def test_discover_single_replicate_dir(tmp_path):
    rep = tmp_path / "rep_single"
    rep.mkdir()
    (rep / "traj.nc").write_text("", encoding="utf-8")
    (rep / "sys.parm7").write_text("", encoding="utf-8")
    reps = discover_replicate_dirs(rep)
    assert reps == [rep.resolve()]


def test_detect_gromacs_tpr_preferred(tmp_path):
    rep = tmp_path / "rep_gmx"
    rep.mkdir()
    (rep / "traj.xtc").write_text("", encoding="utf-8")
    (rep / "topol.tpr").write_text("", encoding="utf-8")
    out = detect_replicate_input(rep)
    assert out.engine == "gromacs"
    assert out.topology_path.name == "topol.tpr"


def test_load_trajectory_with_tpr_topology():
    """
    .tpr topology is read via MDAnalysis (no GROMACS binary required) and handed
    to MDTraj as a converted PDB topology.
    """
    pytest.importorskip("MDAnalysis")
    datafiles = pytest.importorskip("MDAnalysisTests.datafiles")

    traj = load_trajectory(datafiles.XTC, top=datafiles.TPR)
    assert traj.n_frames > 0
    assert traj.n_atoms > 0
    assert traj.topology.n_bonds > 0


def test_load_trajectory_tpr_without_mdanalysis_raises(monkeypatch, tmp_path):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "MDAnalysis":
            raise ImportError("no MDAnalysis")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    tpr = tmp_path / "topol.tpr"
    tpr.write_text("", encoding="utf-8")
    xtc = tmp_path / "traj.xtc"
    xtc.write_text("", encoding="utf-8")

    with pytest.raises(ValueError, match="phenoms\\[gromacs\\]"):
        load_trajectory(str(xtc), top=str(tpr))


def _make_synthetic_gromacs_replicate(tmp_path, n_frames=20, dt_ps=100.0):
    """A small multi-frame .xtc + .pdb topology replicate with distinct,
    deterministic per-frame coordinates and explicit frame times."""
    base = md.load(MINIMAL_PDB)
    xyz = np.tile(base.xyz, (n_frames, 1, 1)).astype(np.float32)
    rng = np.random.default_rng(0)
    xyz = xyz + rng.normal(scale=0.05, size=xyz.shape).astype(np.float32)
    time = (np.arange(n_frames, dtype=np.float32) * dt_ps)
    unitcell_lengths = np.tile([2.0, 2.0, 2.0], (n_frames, 1)).astype(np.float32)
    unitcell_angles = np.tile([90.0, 90.0, 90.0], (n_frames, 1)).astype(np.float32)
    traj = md.Trajectory(
        xyz, base.topology, time=time,
        unitcell_lengths=unitcell_lengths, unitcell_angles=unitcell_angles,
    )
    rep_dir = tmp_path / "rep1"
    rep_dir.mkdir()
    traj.save_xtc(str(rep_dir / "traj.xtc"))
    base.save_pdb(str(rep_dir / "top.pdb"))
    return rep_dir


class TestNormalizeReplicateStreaming:
    """
    normalize_replicate_to_pdb streams the raw trajectory via mdtraj.iterload
    instead of loading it whole, to avoid OOM on large raw trajectories (the
    biggest file a user has, since it's read before frame_dt_ps thinning).
    Chunk size must not change the result: fit/thinning are per-frame
    operations, so chunking is purely a memory-use detail.
    """

    def test_chunk_size_does_not_change_output(self, tmp_path):
        rep_dir = _make_synthetic_gromacs_replicate(tmp_path)
        rep = detect_replicate_input(rep_dir)

        out_small_chunks = tmp_path / "small_chunks.pdb"
        out_one_chunk = tmp_path / "one_chunk.pdb"
        normalize_replicate_to_pdb(
            rep, out_small_chunks, frame_dt_ps=300.0, apply_imaging=False, chunk_size=3
        )
        normalize_replicate_to_pdb(
            rep, out_one_chunk, frame_dt_ps=300.0, apply_imaging=False, chunk_size=1000
        )

        t_small = md.load(str(out_small_chunks))
        t_one = md.load(str(out_one_chunk))
        assert t_small.n_frames == t_one.n_frames == 7  # times 0,300,...,1800
        assert np.allclose(t_small.xyz, t_one.xyz, atol=1e-6)

    def test_frame_dt_thinning_keeps_expected_frame_count(self, tmp_path):
        rep_dir = _make_synthetic_gromacs_replicate(tmp_path)
        rep = detect_replicate_input(rep_dir)
        out = tmp_path / "thinned.pdb"
        normalize_replicate_to_pdb(
            rep, out,
            frame_dt_ps=300.0, start_ps=200.0, end_ps=1600.0,
            apply_imaging=True, fit=True, chunk_size=3,
        )
        traj = md.load(str(out))
        # start_ps=200 drops t=0,100; frame_dt=300 then keeps 200,500,800,1100,1400; end_ps=1600 is a no-op here.
        assert traj.n_frames == 5


def test_discover_set_with_subdirs(tmp_path):
    root = tmp_path / "wt_set"
    rep_a = root / "rep_a"
    rep_b = root / "rep_b"
    rep_a.mkdir(parents=True)
    rep_b.mkdir(parents=True)
    (rep_a / "traj.dcd").write_text("", encoding="utf-8")
    (rep_a / "top.pdb").write_text("MODEL        1\nENDMDL\nEND\n", encoding="utf-8")
    (rep_b / "traj.nc").write_text("", encoding="utf-8")
    (rep_b / "sys.prmtop").write_text("", encoding="utf-8")
    reps = discover_replicate_dirs(root)
    assert reps == [rep_a.resolve(), rep_b.resolve()]
