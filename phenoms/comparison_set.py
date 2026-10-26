"""
Two simulation sets: differential protection, comparison CSV, difference plot, PDB B-factors, PCA/t-SNE/Isomap.

Use plot_heatmaps_both() to get heatmaps for both sets with the same H-bonds in the same order
(missing bonds in one set shown as 0). Single-set heatmaps remain via SimulationSet.plot_heatmaps().
"""

import json
from pathlib import Path

from phenoms.analysis import extract_residue_numbers
from phenoms.outputs import timestamped_run_dir
from phenoms.plotting import plot_difference, plot_heatmap, plot_heatmap_with_legend


class ComparisonSet:
    """Two simulation sets (e.g. no ligand vs ligand): compare, plot difference, write PDB, run PCA/t-SNE/Isomap."""

    def __init__(self, set_a, set_b, label_a="set_a", label_b="set_b", *, output_dir=None):
        """
        Parameters
        ----------
        set_a, set_b : SimulationSet
            Both must already be ``.run()``.
        label_a, label_b : str
            Column names / plot labels for the two sets.
        output_dir : str, pathlib.Path, False, or None
            Where the standard comparison bundle (comparison CSV, manifest,
            difference plot, aligned heatmaps, a differential
            structure-colored PDB plus a standalone HTML viewer for it, and
            an interactive HTML residue connectivity Δ graph) is written
            after a successful :meth:`compare`. ``None`` (default) creates a fresh timestamped
            directory under :func:`phenoms.default_output_root`, matching
            :class:`~phenoms.SimulationSet`. Pass an explicit path to control
            the location, or ``False`` to disable all default output writing.
        """
        self.set_a = set_a
        self.set_b = set_b
        self.label_a = label_a
        self.label_b = label_b
        self._comparison_df = None
        self._comparison_params = None
        if output_dir is False:
            self.output_dir = None
        elif output_dir is None:
            self.output_dir = timestamped_run_dir("comparisonset")
        else:
            self.output_dir = Path(output_dir).expanduser().resolve()

    def _require_both_run(self):
        if not self.set_a.get_pivot_tables() or not self.set_b.get_pivot_tables():
            raise ValueError("Both SimulationSets must be run first (.run()).")

    def compare(
        self,
        impute_threshold=0.4,
        flip_difference=False,
        donor_aggregation="sum",
    ):
        self._require_both_run()
        from phenoms.comparison import compare_two_sets
        self._comparison_df = compare_two_sets(
            self.set_a.get_pivot_tables(),
            self.set_b.get_pivot_tables(),
            label_a=self.label_a,
            label_b=self.label_b,
            impute_threshold=impute_threshold,
            flip_difference=flip_difference,
            donor_aggregation=donor_aggregation,
        )
        self._comparison_params = {
            "impute_threshold": impute_threshold,
            "flip_difference": flip_difference,
            "donor_aggregation": donor_aggregation,
        }
        if self.output_dir is not None:
            self._write_default_outputs()
        return self._comparison_df

    def _write_default_outputs(self):
        """
        Write the standard comparison bundle to ``self.output_dir``:
        ``comparison.csv`` + manifest (:meth:`export_comparison_artifacts`), a
        difference plot, aligned heatmaps for both sets, a differential
        structure-colored PDB plus a standalone HTML viewer for it (see
        :meth:`_write_default_structure_bfactors`), and an interactive HTML
        residue connectivity Δ graph (:meth:`export_connectivity_graph_html`).
        """
        self.export_comparison_artifacts(self.output_dir)
        plots_dir = self.output_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        self.plot_difference(save_path=str(plots_dir / "difference.png"))
        heatmaps_dir = plots_dir / "heatmaps"
        heatmaps_dir.mkdir(parents=True, exist_ok=True)
        self.plot_heatmaps_both(save_dir=str(heatmaps_dir))

        bfactor_pdb = self.output_dir / "structure_bfactors_diff.pdb"
        self._write_default_structure_bfactors(bfactor_pdb)

        from phenoms.viewer_bfactor_html import export_bfactor_colored_pdb_viewer_html

        export_bfactor_colored_pdb_viewer_html(
            bfactor_pdb,
            self.output_dir / "structure_bfactors_diff_viewer.html",
            title=f"Δ H-bond occupancy ({self.label_b} vs {self.label_a}, colored by B-factor)",
        )
        self.export_connectivity_graph_html(
            self.output_dir / "connectivity_graph.html", graph_mode="diff"
        )

    def _write_default_structure_bfactors(self, output_path):
        """
        Color a reference structure by per-residue occupancy difference
        (``label_b`` − ``label_a``, per :meth:`compare`'s ``flip_difference``).
        The reference structure comes from ``set_a``, via the same frame-0
        convention as :class:`~phenoms.SimulationSet`'s own default output
        (see :meth:`~phenoms.SimulationSet._resolve_reference_pdb`).
        """
        ref_pdb = self.set_a._resolve_reference_pdb(self.output_dir)
        self.write_pdb_bfactors(str(ref_pdb), str(output_path))

    def export_connectivity_graph_html(
        self,
        output_html_path,
        *,
        graph_mode: str = "diff",
        top_k_edges: int = 30,
        min_abs_edge_delta: float | None = None,
        directed: bool = False,
        title: str | None = None,
    ) -> None:
        """
        Export an interactive HTML residue connectivity graph derived from bond-level deltas.

        graph_mode:
          - "diff": edges sized by abs(Δ occupancy) and colored by sign (sign consistent with flip_difference in compare()).
          - "set_a" / "set_b": edges sized by occupancy (no sign).
        """
        self._require_both_run()
        from phenoms.connectivity import (
            bond_delta_from_pivot_tables,
            bond_occupancy_from_pivot_tables,
            build_residue_graph_from_bond_deltas,
            export_residue_graph_html,
        )

        unified = self.get_unified_bond_labels()
        # Use plot-filtered pivot tables (region-aware), then align to the unified bond list.
        pivots_a = [pt.reindex(unified, fill_value=0) for pt in self.set_a.get_plot_pivot_tables()]
        pivots_b = [pt.reindex(unified, fill_value=0) for pt in self.set_b.get_plot_pivot_tables()]

        flip_difference = False
        impute_threshold = None
        if self._comparison_params is not None:
            flip_difference = bool(self._comparison_params.get("flip_difference", False))
            impute_threshold = self._comparison_params.get("impute_threshold", None)

        if graph_mode == "diff":
            if title is None:
                title = f"Residue H-bond connectivity Δ ({self.label_b} vs {self.label_a})"
            bond_values = bond_delta_from_pivot_tables(
                pivots_a,
                pivots_b,
                unified,
                flip_difference=flip_difference,
                impute_threshold=impute_threshold,
            )
        elif graph_mode == "set_a":
            if title is None:
                title = f"Residue H-bond connectivity occupancy ({self.label_a})"
            bond_values = bond_occupancy_from_pivot_tables(pivots_a, unified)
        elif graph_mode == "set_b":
            if title is None:
                title = f"Residue H-bond connectivity occupancy ({self.label_b})"
            bond_values = bond_occupancy_from_pivot_tables(pivots_b, unified)
        else:
            raise ValueError("graph_mode must be one of: 'diff', 'set_a', 'set_b'")

        graph = build_residue_graph_from_bond_deltas(
            bond_values,
            top_k_edges=top_k_edges,
            min_abs_edge_delta=min_abs_edge_delta,
            directed=directed,
        )
        export_residue_graph_html(
            graph,
            output_html_path,
            title=title,
            directed=directed,
            show_delta_legend=(graph_mode == "diff"),
        )

    def export_connectivity_community_graph_html(
        self,
        output_html_path,
        *,
        graph_mode: str = "diff",
        top_k_edges: int = 40,
        min_abs_edge_delta: float | None = None,
        directed: bool = False,
        resolution: float = 1.0,
        title: str | None = None,
        community_nodes_csv_path: str | None = None,
        community_summary_csv_path: str | None = None,
    ) -> None:
        """
        Export a community-aware connectivity graph and optional community CSV summaries.

        graph_mode:
          - "diff": edges use delta occupancy (sign preserved).
          - "set_a"/"set_b": edges use occupancy in that set.
        """
        self._require_both_run()
        from phenoms.connectivity import (
            bond_delta_from_pivot_tables,
            bond_occupancy_from_pivot_tables,
            build_residue_graph_with_communities,
            community_tables,
            export_residue_graph_html,
        )

        unified = self.get_unified_bond_labels()
        pivots_a = [pt.reindex(unified, fill_value=0) for pt in self.set_a.get_plot_pivot_tables()]
        pivots_b = [pt.reindex(unified, fill_value=0) for pt in self.set_b.get_plot_pivot_tables()]

        flip_difference = False
        impute_threshold = None
        if self._comparison_params is not None:
            flip_difference = bool(self._comparison_params.get("flip_difference", False))
            impute_threshold = self._comparison_params.get("impute_threshold", None)

        if graph_mode == "diff":
            if title is None:
                title = f"Community H-bond connectivity Δ ({self.label_b} vs {self.label_a})"
            bond_values = bond_delta_from_pivot_tables(
                pivots_a,
                pivots_b,
                unified,
                flip_difference=flip_difference,
                impute_threshold=impute_threshold,
            )
        elif graph_mode == "set_a":
            if title is None:
                title = f"Community H-bond connectivity occupancy ({self.label_a})"
            bond_values = bond_occupancy_from_pivot_tables(pivots_a, unified)
        elif graph_mode == "set_b":
            if title is None:
                title = f"Community H-bond connectivity occupancy ({self.label_b})"
            bond_values = bond_occupancy_from_pivot_tables(pivots_b, unified)
        else:
            raise ValueError("graph_mode must be one of: 'diff', 'set_a', 'set_b'")

        result = build_residue_graph_with_communities(
            bond_values,
            top_k_edges=top_k_edges,
            min_abs_edge_delta=min_abs_edge_delta,
            directed=directed,
            resolution=resolution,
        )
        node_df, summary_df = community_tables(result)
        node_df["graph_mode"] = graph_mode
        node_df["label_a"] = self.label_a
        node_df["label_b"] = self.label_b
        summary_df["graph_mode"] = graph_mode
        summary_df["label_a"] = self.label_a
        summary_df["label_b"] = self.label_b

        if community_nodes_csv_path:
            node_df.to_csv(community_nodes_csv_path, index=False)
        if community_summary_csv_path:
            summary_df.to_csv(community_summary_csv_path, index=False)

        summary_stats = {
            "n_nodes": len(result.graph.nodes),
            "n_edges": len(result.graph.edges),
            "n_communities": len(set(result.community_by_node.values())),
            "modularity": float(result.modularity),
        }
        export_residue_graph_html(
            result.graph,
            output_html_path,
            title=title,
            directed=directed,
            summary_stats=summary_stats,
            show_delta_legend=(graph_mode == "diff"),
        )

    def get_comparison_df(self):
        if self._comparison_df is None:
            self.compare()
        return self._comparison_df

    def get_unified_bond_labels(self):
        """Union of all bond labels from both sets, sorted (for aligned heatmaps)."""
        self._require_both_run()
        if self.set_a.resid_range != self.set_b.resid_range:
            raise ValueError(
                "plotting resid_range must match between set_a and set_b for aligned heatmaps/PCA."
            )
        bonds_a = set(self.set_a.get_plot_bond_labels_sorted())
        bonds_b = set(self.set_b.get_plot_bond_labels_sorted())
        return sorted(bonds_a | bonds_b, key=extract_residue_numbers)

    def plot_heatmaps_both(self, save_dir=None, use_legend=False):
        """
        Plot heatmaps for all replicates in both sets using a unified bond list,
        so the same H-bonds appear in the same order on every heatmap (missing bonds = 0).
        """
        self._require_both_run()
        unified = self.get_unified_bond_labels()
        region_str = self.set_a._plot_region_str()
        plot_fn = plot_heatmap_with_legend if use_legend else plot_heatmap
        for i, pt in enumerate(self.set_a.get_plot_pivot_tables()):
            pt_aligned = pt.reindex(unified, fill_value=0)
            stem = Path(self.set_a.pdb_files[i]).stem if i < len(self.set_a.pdb_files) else f"{self.label_a}_rep{i+1}"
            title = f"H-Bond Occupation - {stem} ({region_str})" if use_legend else f"{stem} ({region_str})"
            path = f"{save_dir}/{stem}.png" if save_dir else None
            plot_fn(pt_aligned, title, save_path=path)
        for i, pt in enumerate(self.set_b.get_plot_pivot_tables()):
            pt_aligned = pt.reindex(unified, fill_value=0)
            stem = Path(self.set_b.pdb_files[i]).stem if i < len(self.set_b.pdb_files) else f"{self.label_b}_rep{i+1}"
            title = f"H-Bond Occupation - {stem} ({region_str})" if use_legend else f"{stem} ({region_str})"
            path = f"{save_dir}/{stem}.png" if save_dir else None
            plot_fn(pt_aligned, title, save_path=path)

    def save_comparison_csv(self, path):
        self.get_comparison_df().to_csv(path, index=False)

    def export_comparison_artifacts(self, output_dir):
        """
        Write ``raw_data/comparison.csv`` (bond-level Δ table) and ``manifest.json``
        under ``output_dir``. Call after :meth:`compare` (or rely on :meth:`get_comparison_df`
        which triggers ``compare()`` with defaults).
        """
        output_dir = Path(output_dir).expanduser().resolve()
        raw = output_dir / "raw_data"
        raw.mkdir(parents=True, exist_ok=True)
        self.get_comparison_df().to_csv(raw / "comparison.csv", index=False)
        unified = self.get_unified_bond_labels()
        manifest = {
            "label_a": self.label_a,
            "label_b": self.label_b,
            "n_unified_bonds": len(unified),
            "comparison_params": self._comparison_params,
        }
        (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")

    def plot_difference(
        self,
        diff_threshold=0.2,
        title=None,
        save_path=None,
        diff_threshold_mode="manual",
        impute_small_differences=True,
        show_labels=True,
        interpolate=True,
        interpolation_points=900,
        autocorr_kappa=1.96,
        return_meta=False,
    ):
        df = self.get_comparison_df()
        if title is None:
            title = f"Change in Protection ({self.label_b} − {self.label_a}) (Entire protein)"
        return plot_difference(
            df,
            title=title,
            diff_threshold=diff_threshold,
            save_path=save_path,
            diff_threshold_mode=diff_threshold_mode,
            impute_small_differences=impute_small_differences,
            show_labels=show_labels,
            interpolate=interpolate,
            interpolation_points=interpolation_points,
            autocorr_kappa=autocorr_kappa,
            return_meta=return_meta,
        )

    def write_pdb_bfactors(self, pdb_path, output_path, value_column="Difference_clipped"):
        from phenoms.structure import write_pdb_bfactors
        write_pdb_bfactors(pdb_path, self.get_comparison_df(), output_path, value_column=value_column)

    def run_pca(self, group_labels=None, plot=True, title="PCA"):
        from phenoms.dimensionality_reduction import aggregate_replicate_data, run_pca as _run_pca, plot_manifold
        pivots_a, pivots_b = self.set_a.get_plot_pivot_tables(), self.set_b.get_plot_pivot_tables()
        all_pivots = pivots_a + pivots_b
        group_labels = group_labels or [self.label_a] * len(pivots_a) + [self.label_b] * len(pivots_b)
        data = aggregate_replicate_data(all_pivots)
        scores, var_ratio = _run_pca(data, n_components=2)
        if plot:
            plot_manifold(scores, group_labels, title=f"{title} ({self.set_a._plot_region_str()})")
        return scores, var_ratio

    def run_tsne(self, group_labels=None, perplexity=2, plot=True, title="t-SNE"):
        from phenoms.dimensionality_reduction import aggregate_replicate_data, run_tsne, plot_manifold
        pivots_a, pivots_b = self.set_a.get_plot_pivot_tables(), self.set_b.get_plot_pivot_tables()
        all_pivots = pivots_a + pivots_b
        group_labels = group_labels or [self.label_a] * len(pivots_a) + [self.label_b] * len(pivots_b)
        data = aggregate_replicate_data(all_pivots)
        scores = run_tsne(data, perplexity=perplexity)
        if plot:
            plot_manifold(scores, group_labels, title=f"{title} ({self.set_a._plot_region_str()})")
        return scores

    def run_isomap(self, group_labels=None, n_neighbors=5, plot=True, title="Isomap"):
        from phenoms.dimensionality_reduction import aggregate_replicate_data, run_isomap, plot_manifold
        pivots_a, pivots_b = self.set_a.get_plot_pivot_tables(), self.set_b.get_plot_pivot_tables()
        all_pivots = pivots_a + pivots_b
        group_labels = group_labels or [self.label_a] * len(pivots_a) + [self.label_b] * len(pivots_b)
        data = aggregate_replicate_data(all_pivots)
        scores = run_isomap(data, n_neighbors=n_neighbors)
        if plot:
            plot_manifold(scores, group_labels, title=f"{title} ({self.set_a._plot_region_str()})")
        return scores

    def run_manifold_suite(self, group_labels=None, perplexity=2, n_neighbors=5, random_state=42, plot=True):
        """PCA, t-SNE, Isomap on both sets’ replicates (two-group coloring by default)."""
        from phenoms.dimensionality_reduction import run_manifold_suite as _suite
        pivots_a, pivots_b = self.set_a.get_plot_pivot_tables(), self.set_b.get_plot_pivot_tables()
        all_pivots = pivots_a + pivots_b
        group_labels = group_labels or [self.label_a] * len(pivots_a) + [self.label_b] * len(pivots_b)
        return _suite(
            all_pivots,
            group_labels,
            perplexity=perplexity,
            n_neighbors=n_neighbors,
            random_state=random_state,
            plot=plot,
        )
