"""
Tests for SimulationSet (run, getters, optional bond statistics).
"""

import os
import mdtraj as md
import numpy as np
import pytest
import pandas as pd

from phenoms import SimulationSet, default_output_root
from phenoms.io import load_and_select_residues
from phenoms.hbond import process_frames

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
MINIMAL_PDB = os.path.join(FIXTURES_DIR, "minimal.pdb")
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KRAS_PDB = os.path.join(_REPO_ROOT, "H_bond_work", "rep1_md_Kras_only_500ns_none_nojump_center.pdb")


@pytest.fixture
def minimal_pdb_path():
    return MINIMAL_PDB


@pytest.fixture
def kras_pdb_path():
    """Path to Kras trajectory (501 frames). Skip if file not present."""
    if not os.path.isfile(KRAS_PDB):
        pytest.skip(f"Kras trajectory not found: {KRAS_PDB}")
    return KRAS_PDB


class TestSimulationSet:
    def test_init_single_file(self):
        sim = SimulationSet(["a.pdb"], resid_range=(1, 10), sub_frames=5)
        assert len(sim.pdb_files) == 1
        assert sim.resid_range == (1, 10)
        assert sim.sub_frames == 5

    def test_init_replicates(self):
        sim = SimulationSet(["a.pdb", "b.pdb"], resid_range=None)
        assert len(sim.pdb_files) == 2
        assert sim.resid_range is None

    def test_run_whole_protein_single_replicate(self, minimal_pdb_path):
        sim = SimulationSet(
            [minimal_pdb_path],
            resid_range=None,
            sub_frames=1,
        )
        try:
            sim.run(n_jobs=1)
        except ValueError as e:
            if "No bonds" in str(e):
                pytest.skip("Trajectory has no bonds (minimal fixture has no H atoms)")
            raise
        assert len(sim.get_pivot_tables()) == 1
        assert len(sim.get_hbond_dfs()) == 1
        assert sim.get_bond_labels_sorted() is not None
        assert sim.bond_statistics is None

    def test_run_with_bond_statistics_threshold(self, minimal_pdb_path):
        sim = SimulationSet(
            [minimal_pdb_path],
            resid_range=None,
            sub_frames=1,
            bond_statistics_threshold=0.5,
        )
        try:
            sim.run(n_jobs=1)
        except ValueError as e:
            if "No bonds" in str(e):
                pytest.skip("Trajectory has no bonds (minimal fixture has no H atoms)")
            raise
        assert sim.bond_statistics is not None
        assert "mean_lifetimes" in sim.bond_statistics
        assert "std_lifetimes" in sim.bond_statistics
        assert "mean_break_frequencies" in sim.bond_statistics
        assert "std_break_frequencies" in sim.bond_statistics

    def test_get_pivot_tables_before_run_raises(self):
        sim = SimulationSet(["a.pdb"])
        with pytest.raises(ValueError, match="No data"):
            sim.plot_heatmaps()

    def test_plot_aggregated_heatmap_after_run(self, minimal_pdb_path):
        sim = SimulationSet([minimal_pdb_path], sub_frames=1)
        try:
            sim.run(n_jobs=1)
        except ValueError as e:
            if "No bonds" in str(e):
                pytest.skip("Trajectory has no bonds (minimal fixture has no H atoms)")
            raise
        # Should not raise; we don't capture plt here
        sim.plot_aggregated_heatmap(save_path=os.path.join(FIXTURES_DIR, "agg.png"))
        if os.path.exists(os.path.join(FIXTURES_DIR, "agg.png")):
            os.remove(os.path.join(FIXTURES_DIR, "agg.png"))

    def test_run_kras_first_10ns(self, kras_pdb_path):
        """Integration test: run on Kras rep1, first 100 frames (~10 ns)."""
        sim = SimulationSet(
            [kras_pdb_path],
            resid_range=None,
            sub_frames=100,
        )
        sim.run(n_jobs=1)
        assert len(sim.get_pivot_tables()) == 1
        assert len(sim.get_hbond_dfs()) == 1
        pivot = sim.get_pivot_tables()[0]
        assert pivot.shape[0] > 0 and pivot.shape[1] == 100

    def test_qc_sees_full_trajectory_even_when_sub_frames_caps_detection(self, tmp_path):
        """
        run(sub_frames=...) streams only the first sub_frames frames off disk for
        detection (memory optimization), but qc=True must still see every frame —
        RMSD convergence is meaningless over an artificially truncated series.
        """
        base = md.load(MINIMAL_PDB)
        n_frames = 15
        rng = np.random.default_rng(0)
        xyz = np.repeat(base.xyz, n_frames, axis=0).astype(np.float32)
        xyz += rng.normal(scale=0.001, size=xyz.shape).astype(np.float32)
        full = md.Trajectory(xyz, base.topology)
        pdb_path = tmp_path / "rep.pdb"
        full.save_pdb(str(pdb_path))

        sim = SimulationSet([str(pdb_path)], sub_frames=3)
        try:
            sim.run(n_jobs=1, qc=True, qc_fail_on_nonconverged=False)
        except ValueError as e:
            if "No bonds" in str(e):
                pytest.skip("Trajectory has no bonds (minimal fixture has no H atoms)")
            raise
        rmsd_entry = sim.get_qc_report()["rmsd_convergence"][0]
        assert rmsd_entry["n_points"] == n_frames
        assert sim.get_pivot_tables()[0].shape[1] == 3

    def test_run_kras_region(self, kras_pdb_path):
        """Integration test: run on Kras rep1 with residue range, first 50 frames."""
        sim = SimulationSet(
            [kras_pdb_path],
            resid_range=(1, 50),
            sub_frames=50,
        )
        sim.run(n_jobs=1)
        assert len(sim.get_pivot_tables()) == 1
        assert len(sim.get_hbond_dfs()) == 1


class TestDefaultOutputs:
    """
    SimulationSet.run() writes a standard artifact bundle (tables, plots, a
    structure-colored reference PDB) by default, to a fresh timestamped
    directory when output_dir isn't given. output_dir=False opts back out to
    pure in-memory results (the pre-existing behavior).
    """

    def _bundle_paths(self, output_dir):
        return {
            "raw_data": output_dir / "raw_data",
            "manifest": output_dir / "manifest.json",
            "aggregated_heatmap": output_dir / "plots" / "aggregated_heatmap.png",
            "structure_bfactors": output_dir / "structure_bfactors.pdb",
            "structure_bfactors_viewer_html": output_dir / "structure_bfactors_viewer.html",
            "connectivity_graph_html": output_dir / "connectivity_graph.html",
        }

    def test_default_output_dir_is_timestamped_under_default_root(self, kras_pdb_path):
        sim = SimulationSet([kras_pdb_path], sub_frames=5)
        assert sim.output_dir.parent == default_output_root()
        assert sim.output_dir.name.startswith("simulationset_")

    def test_run_writes_default_bundle(self, kras_pdb_path):
        sim = SimulationSet([kras_pdb_path], sub_frames=5)
        sim.run(n_jobs=1)
        for label, path in self._bundle_paths(sim.output_dir).items():
            assert path.exists(), f"missing default output: {label} ({path})"
        viewer_html = (sim.output_dir / "structure_bfactors_viewer.html").read_text(encoding="utf-8")
        assert "NGL" in viewer_html
        graph_html = (sim.output_dir / "connectivity_graph.html").read_text(encoding="utf-8")
        assert "vis-network" in graph_html

    def test_output_dir_false_disables_writing(self, kras_pdb_path):
        sim = SimulationSet([kras_pdb_path], sub_frames=5, output_dir=False)
        sim.run(n_jobs=1)
        assert sim.output_dir is None

    def test_explicit_output_dir_is_used_and_resolved(self, kras_pdb_path, tmp_path):
        target = tmp_path / "custom_run"
        sim = SimulationSet([kras_pdb_path], sub_frames=5, output_dir=target)
        assert sim.output_dir == target.resolve()
        sim.run(n_jobs=1)
        for path in self._bundle_paths(sim.output_dir).values():
            assert path.exists()

    def test_structure_bfactors_reference_is_frame_0_of_input_pdb(self, kras_pdb_path, tmp_path):
        """
        For pdb_files= input, the structure-coloring reference is the input PDB
        itself (frame 0) -- no separate reference file needs to be written.
        """
        sim = SimulationSet([kras_pdb_path], sub_frames=5, output_dir=tmp_path / "run")
        sim.run(n_jobs=1)
        assert not (sim.output_dir / "reference_frame0.pdb").exists()
        assert (sim.output_dir / "structure_bfactors.pdb").exists()

    def test_structure_bfactors_reference_is_frame_0_for_trajectory_input(self, tmp_path):
        """
        For trajectories=/topology= input there is no ready reference PDB, so
        frame 0 of the trajectory is written out as reference_frame0.pdb and
        used as the coloring target.
        """
        base = md.load(MINIMAL_PDB)
        xtc_path = tmp_path / "traj.xtc"
        pdb_path = tmp_path / "top.pdb"
        base.save_xtc(str(xtc_path))
        base.save_pdb(str(pdb_path))

        sim = SimulationSet.from_trajectories(
            [str(xtc_path)], topology=str(pdb_path), output_dir=tmp_path / "run"
        )
        # Skip real detection (minimal fixture has no bonds) and drive the
        # output-writing step directly with a fabricated pivot table.
        pivot = pd.DataFrame(
            {0: [1, 0], 1: [1, 1], 2: [0, 1]},
            index=["ALA1 -- ALA1", "GLY2 -- ALA1"],
        )
        sim._pivot_tables = [pivot]

        out_path = sim.output_dir / "structure_bfactors.pdb"
        sim.output_dir.mkdir(parents=True, exist_ok=True)
        sim._write_default_structure_bfactors(out_path)

        assert (sim.output_dir / "reference_frame0.pdb").exists()
        assert out_path.exists()
