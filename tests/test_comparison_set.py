import os

import pandas as pd
import pytest

from phenoms import SimulationSet, default_output_root
from phenoms.comparison_set import ComparisonSet

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KRAS_PDB = os.path.join(_REPO_ROOT, "H_bond_work", "rep1_md_Kras_only_500ns_none_nojump_center.pdb")


@pytest.fixture
def kras_pdb_path():
    if not os.path.isfile(KRAS_PDB):
        pytest.skip(f"Kras trajectory not found: {KRAS_PDB}")
    return KRAS_PDB


class _DummySet:
    def __init__(self, pivot_tables, resid_range=None):
        self._pivot_tables = pivot_tables
        self.resid_range = resid_range
        self.pdb_files = ["dummy.pdb"]

    def get_pivot_tables(self):
        return self._pivot_tables

    def get_plot_pivot_tables(self):
        return self._pivot_tables

    def get_plot_bond_labels_sorted(self):
        labels = set()
        for pt in self._pivot_tables:
            labels |= set(pt.index.tolist())
        return sorted(labels)

    def _plot_region_str(self):
        return "Entire protein"


def test_export_connectivity_community_graph_html_writes_outputs(tmp_path):
    labels = ["ALA1 -- GLY2", "GLY2 -- SER3", "THR10 -- VAL11"]
    frames = [0, 1, 2]
    pt_a = pd.DataFrame(
        [[1, 1, 0], [1, 0, 0], [0, 0, 1]],
        index=labels,
        columns=frames,
    )
    pt_b = pd.DataFrame(
        [[0, 1, 1], [0, 1, 1], [1, 1, 0]],
        index=labels,
        columns=frames,
    )
    set_a = _DummySet([pt_a])
    set_b = _DummySet([pt_b])
    comp = ComparisonSet(set_a, set_b, label_a="a", label_b="b", output_dir=False)
    comp.compare(flip_difference=True)

    html_out = tmp_path / "community_diff.html"
    nodes_csv = tmp_path / "community_nodes_diff.csv"
    summary_csv = tmp_path / "community_summary_diff.csv"
    comp.export_connectivity_community_graph_html(
        str(html_out),
        graph_mode="diff",
        community_nodes_csv_path=str(nodes_csv),
        community_summary_csv_path=str(summary_csv),
    )

    assert html_out.exists()
    assert nodes_csv.exists()
    assert summary_csv.exists()
    assert not pd.read_csv(nodes_csv).empty
    assert not pd.read_csv(summary_csv).empty


class TestDefaultOutputs:
    """
    ComparisonSet.compare() writes a standard bundle by default, mirroring
    SimulationSet.run(): comparison CSV/manifest, a difference plot, aligned
    heatmaps, and a differential structure-colored PDB.
    """

    def _bundle_paths(self, output_dir):
        return {
            "raw_data": output_dir / "raw_data" / "comparison.csv",
            "manifest": output_dir / "manifest.json",
            "difference_plot": output_dir / "plots" / "difference.png",
            "structure_bfactors_diff": output_dir / "structure_bfactors_diff.pdb",
            "structure_bfactors_diff_viewer_html": output_dir / "structure_bfactors_diff_viewer.html",
            "connectivity_graph_html": output_dir / "connectivity_graph.html",
        }

    def _built_sets(self, kras_pdb_path):
        a = SimulationSet([kras_pdb_path], resid_range=(1, 50), sub_frames=5, output_dir=False)
        b = SimulationSet([kras_pdb_path], resid_range=(1, 50), sub_frames=5, output_dir=False)
        a.run(n_jobs=1)
        b.run(n_jobs=1)
        return a, b

    def test_default_output_dir_is_timestamped_under_default_root(self, kras_pdb_path):
        a, b = self._built_sets(kras_pdb_path)
        comp = ComparisonSet(a, b, label_a="a", label_b="b")
        assert comp.output_dir.parent == default_output_root()
        assert comp.output_dir.name.startswith("comparisonset_")

    def test_compare_writes_default_bundle(self, kras_pdb_path):
        a, b = self._built_sets(kras_pdb_path)
        comp = ComparisonSet(a, b, label_a="a", label_b="b")
        comp.compare()
        for label, path in self._bundle_paths(comp.output_dir).items():
            assert path.exists(), f"missing default comparison output: {label} ({path})"
        graph_html = (comp.output_dir / "connectivity_graph.html").read_text(encoding="utf-8")
        assert "Edge" in graph_html  # diff-mode legend, unlike SimulationSet's occupancy-only graph
        heatmap_pngs = list((comp.output_dir / "plots" / "heatmaps").glob("*.png"))
        assert heatmap_pngs, "expected at least one aligned heatmap PNG"

    def test_output_dir_false_disables_writing(self, kras_pdb_path):
        a, b = self._built_sets(kras_pdb_path)
        comp = ComparisonSet(a, b, label_a="a", label_b="b", output_dir=False)
        comp.compare()
        assert comp.output_dir is None

    def test_explicit_output_dir_is_used_and_resolved(self, kras_pdb_path, tmp_path):
        a, b = self._built_sets(kras_pdb_path)
        target = tmp_path / "custom_comparison"
        comp = ComparisonSet(a, b, label_a="a", label_b="b", output_dir=target)
        assert comp.output_dir == target.resolve()
        comp.compare()
        for path in self._bundle_paths(comp.output_dir).values():
            assert path.exists()
