"""
SimulationSet: single or replicate MD runs, region or whole-protein H-bond analysis.
"""

import json
import os
from pathlib import Path

import pandas as pd
from tqdm import tqdm
from phenoms.qc import check_mdp_key_consistency, rmsd_convergence_report


def default_n_jobs():
    """Parallel workers for MDTraj frame loop: all CPUs minus two (minimum 1)."""
    return max(1, (os.cpu_count() or 4) - 2)

from phenoms.io import load_and_select_residues, load_trajectory, normalize_topology_list
from phenoms.outputs import timestamped_run_dir
from phenoms.hbond import hbond_occupancy_table, process_frames, process_frames_all
from phenoms.analysis import extract_residue_numbers, create_pivot_table, calculate_bond_statistics, fluctuating_bonds
from phenoms.plotting import (
    plot_heatmap,
    plot_heatmap_with_legend,
    plot_top_hbonds,
    plot_aggregated_heatmap,
    plot_bond_lifetimes_with_error_bars,
    plot_break_frequencies_with_error_bars,
    plot_difference,
)


class SimulationSet:
    """
    Single or replicate MD simulations with H-bond analysis.

    Simple form (unchanged): pass multi-frame PDB paths via ``pdb_files``.

    Native trajectories: pass ``trajectories=`` plus ``topology=`` / ``topologies=``,
    or use :meth:`from_trajectories`.

    - ``backbone_only=True`` (default): Baker–Hubbard backbone N–O bonds (HDX-style).
    - ``backbone_only=False``: all detected donor/acceptor H-bonds.
    - ``resid_range``: filters heatmaps / manifold plots only; detection still runs on protein.
    """

    def __init__(
        self,
        pdb_files=None,
        resid_range=None,
        sub_frames=None,
        bond_statistics_threshold=None,
        output_dir=None,
        *,
        trajectories=None,
        topology=None,
        topologies=None,
        backbone_only=True,
    ):
        """
        Parameters
        ----------
        pdb_files : list of str or str or None
            Multi-frame PDB paths (simple / recommended form). Mutually exclusive with
            ``trajectories``.
        resid_range : tuple (int, int) or None
            Residue range for plot filtering; None = whole protein.
        sub_frames : int or None
            Number of frames to analyze per replicate. Omit or ``None`` = entire trajectory.
        bond_statistics_threshold : float or None
            If set, compute bond statistics (lifetime, break frequency) and
            mean ± std across replicates; used for bar plots and aggregated heatmap.
        output_dir : str, pathlib.Path, False, or None
            Where the standard artifact bundle (per-replicate CSVs, manifest,
            heatmap/statistics plots, a structure-colored reference PDB plus a
            standalone HTML viewer for it, and an interactive HTML residue
            connectivity graph) is written after a successful :meth:`run`.
            ``None`` (default) creates a
            fresh timestamped directory under :func:`phenoms.default_output_root`
            so a run always leaves output behind without naming a path. Pass an
            explicit path to control the location, or ``False`` to disable all
            default output writing and keep results in memory only (via the
            ``get_*`` accessors).
        trajectories : list of str or str or None
            Native trajectory files (``.xtc``, ``.trr``, ``.dcd``, ``.nc``, …).
            Requires ``topology`` or ``topologies``.
        topology : str or path-like or None
            Shared topology for all ``trajectories``.
        topologies : list of str or None
            Per-replicate topologies (same length as ``trajectories``).
        backbone_only : bool
            If True (default), detect backbone N–O H-bonds only. If False, detect all
            Baker–Hubbard H-bonds.
        """
        if pdb_files is not None and trajectories is not None:
            raise ValueError("Pass either pdb_files= or trajectories=, not both.")
        if pdb_files is None and trajectories is None:
            raise ValueError("Provide pdb_files= (simple form) or trajectories= (native MD).")

        if trajectories is not None:
            if isinstance(trajectories, (str, Path)):
                trajectories = [trajectories]
            input_files = [str(Path(p).expanduser()) for p in trajectories]
            tops = normalize_topology_list(
                len(input_files), topology=topology, topologies=topologies
            )
            if tops is None or any(t is None for t in tops):
                raise ValueError(
                    "Native trajectories require topology= (shared) or topologies= "
                    "(one per replicate)."
                )
            self.topologies = tops
            self.input_kind = "trajectory"
        else:
            if isinstance(pdb_files, (str, Path)):
                pdb_files = [pdb_files]
            input_files = [str(Path(p).expanduser()) for p in pdb_files]
            if topology is not None or topologies is not None:
                raise ValueError("topology=/topologies= only apply with trajectories=.")
            self.topologies = None
            self.input_kind = "pdb"

        # Kept for backward compatibility (artifact stems / replicate labels).
        self.pdb_files = input_files
        self.resid_range = resid_range
        self.sub_frames = sub_frames
        self.bond_statistics_threshold = bond_statistics_threshold
        self.backbone_only = bool(backbone_only)
        if output_dir is False:
            self.output_dir = None
        elif output_dir is None:
            self.output_dir = timestamped_run_dir("simulationset")
        else:
            self.output_dir = Path(output_dir).expanduser().resolve()

        self._hbond_dfs = []
        self._pivot_tables = []
        self._bond_labels_sorted = []
        self._bond_statistics = None
        self._qc_report = None

    @classmethod
    def from_trajectories(
        cls,
        trajectories,
        topology=None,
        topologies=None,
        *,
        resid_range=None,
        sub_frames=None,
        bond_statistics_threshold=None,
        output_dir=None,
        backbone_only=True,
    ):
        """Build a SimulationSet from native MD trajectory + topology file(s)."""
        return cls(
            trajectories=trajectories,
            topology=topology,
            topologies=topologies,
            resid_range=resid_range,
            sub_frames=sub_frames,
            bond_statistics_threshold=bond_statistics_threshold,
            output_dir=output_dir,
            backbone_only=backbone_only,
        )

    def run(
        self,
        n_jobs=None,
        use_rust=True,
        *,
        qc=False,
        mdp_files=None,
        skip_mdp_consistency=False,
        qc_fail_on_nonconverged=True,
        qc_last_fraction=0.2,
        qc_mean_tolerance=0.05,
        qc_slope_tolerance=1e-4,
    ):
        """
        Load all trajectories, run Baker–Hubbard detection, build pivot tables.
        If bond_statistics_threshold was set, compute per-replicate stats and mean ± std.

        Parameters
        ----------
        n_jobs : int or None
            MDTraj parallel workers per replicate (Rust path ignores this). ``None`` =
            :func:`default_n_jobs` (all CPUs minus two).
        use_rust : bool
            If True, use Rust extension when available. Set False if Rust returns no bonds
            for your topology (MDTraj fallback).
        qc : bool
            If True, run a simple replicate QC pass (RMSD convergence and optional MDP checks).
        mdp_files : list[str] or None
            Optional list of .mdp paths to compare for parameter consistency.
        skip_mdp_consistency : bool
            If True, skip MDP checks even if mdp_files are supplied.
        qc_fail_on_nonconverged : bool
            If True, raise RuntimeError when any replicate fails RMSD convergence.
        """
        if n_jobs is None:
            n_jobs = default_n_jobs()
        process_fn = process_frames if self.backbone_only else process_frames_all
        all_hbonds_dfs = []
        all_bond_labels = set()
        qc_report = {
            "enabled": bool(qc),
            "rmsd_convergence": [],
            "all_replicates_converged": True,
            "mdp_consistency": None,
            "backbone_only": self.backbone_only,
        }

        replicate_bar = tqdm(self.pdb_files, desc="Replicates", unit="replicate")
        for i, input_path in enumerate(replicate_bar):
            replicate_bar.set_postfix_str(Path(input_path).name)
            top = None if self.topologies is None else self.topologies[i]
            # QC (e.g. RMSD convergence) needs every frame regardless of sub_frames, so
            # only cap the read when QC is off — otherwise behavior for qc=True is unchanged.
            max_frames = self.sub_frames if not qc else None
            trajectory = load_and_select_residues(
                input_path, resid_range=None, top=top, max_frames=max_frames
            )
            if qc:
                rep = rmsd_convergence_report(
                    trajectory,
                    last_fraction=qc_last_fraction,
                    mean_tolerance=qc_mean_tolerance,
                    slope_tolerance=qc_slope_tolerance,
                )
                rep["input_file"] = str(input_path)
                rep["pdb_file"] = str(input_path)  # backward-compatible key
                qc_report["rmsd_convergence"].append(rep)
                if not rep["converged"]:
                    qc_report["all_replicates_converged"] = False
            hbonds_df = process_fn(
                trajectory,
                sub_frames=self.sub_frames,
                n_jobs=n_jobs,
                use_rust=use_rust,
            )
            all_bond_labels.update(hbonds_df["Bond Label"].unique())
            all_hbonds_dfs.append(hbonds_df)

        bond_labels_sorted = sorted(
            all_bond_labels,
            key=lambda x: extract_residue_numbers(x),
        )
        self._bond_labels_sorted = bond_labels_sorted

        all_pivot_tables = []
        for hbonds_df in all_hbonds_dfs:
            pivot = create_pivot_table(hbonds_df, bond_labels_sorted)
            all_pivot_tables.append(pivot)

        self._hbond_dfs = all_hbonds_dfs
        self._pivot_tables = all_pivot_tables
        self._qc_report = qc_report

        if qc and mdp_files and not skip_mdp_consistency:
            mdp_rep = check_mdp_key_consistency(mdp_files)
            self._qc_report["mdp_consistency"] = mdp_rep
            if not mdp_rep["ok"]:
                details = json.dumps(mdp_rep["mismatches"], indent=2)
                raise RuntimeError(
                    "MDP consistency check failed for key parameters.\n"
                    f"Mismatched keys and values:\n{details}"
                )

        if qc and qc_fail_on_nonconverged and not self._qc_report["all_replicates_converged"]:
            failed = [r for r in self._qc_report["rmsd_convergence"] if not r["converged"]]
            details = "\n".join(
                f"- {x['pdb_file']}: reason={x['reason']}, rel_diff={x.get('relative_mean_diff')}, "
                f"slope={x.get('last_window_slope')}"
                for x in failed
            )
            raise RuntimeError(
                "QC failed: not all replicates passed RMSD convergence.\n"
                f"Failed replicates:\n{details}"
            )

        if self.bond_statistics_threshold is not None:
            self._compute_bond_statistics()

        if self.output_dir is not None:
            self._write_default_outputs()

        return self

    def _write_default_outputs(self):
        """
        Write the standard artifact bundle to ``self.output_dir``: raw/occupancy/
        pivot CSVs and a manifest (:meth:`export_run_artifacts`), per-replicate
        and aggregated heatmaps, bond-statistics bar plots (if
        ``bond_statistics_threshold`` was set), a structure-colored reference
        PDB plus a standalone HTML viewer for it, and an interactive HTML
        residue connectivity graph (:meth:`export_connectivity_graph_html`).
        """
        self.export_run_artifacts(self.output_dir)
        plots_dir = self.output_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        self.plot_heatmaps(save_dir=str(plots_dir))
        self.plot_aggregated_heatmap(save_path=str(plots_dir / "aggregated_heatmap.png"))
        if self.bond_statistics_threshold is not None:
            self.plot_bond_statistics(
                save_lifetime_path=str(plots_dir / "bond_lifetimes.png"),
                save_break_path=str(plots_dir / "break_frequencies.png"),
            )
        bfactor_pdb = self.output_dir / "structure_bfactors.pdb"
        self._write_default_structure_bfactors(bfactor_pdb)

        from phenoms.viewer_bfactor_html import export_bfactor_colored_pdb_viewer_html

        export_bfactor_colored_pdb_viewer_html(
            bfactor_pdb,
            self.output_dir / "structure_bfactors_viewer.html",
            title="H-bond occupancy variance (colored by B-factor)",
        )
        self.export_connectivity_graph_html(self.output_dir / "connectivity_graph.html")

    def _resolve_reference_pdb(self, extra_output_dir):
        """
        Path to a single-frame reference structure for B-factor coloring: the
        first replicate's PDB for ``pdb_files=`` input (frame 0 embedded), or
        frame 0 of the first replicate's trajectory for ``trajectories=`` input
        (no ready PDB exists there, so it's written out as
        ``reference_frame0.pdb`` under ``extra_output_dir``).

        PHENOMS does not auto-trim equilibration on either path, so if your
        simulation includes a burn-in period, either pre-trim it first (e.g.
        ``phenoms prep --start-ps``) or treat this reference structure as a
        coloring target only — it plays no part in H-bond detection, which
        already ran on the full requested frame range. Shared with
        :class:`~phenoms.ComparisonSet`'s own default structure coloring.
        """
        if self.input_kind == "pdb":
            return Path(self.pdb_files[0])
        top = None if self.topologies is None else self.topologies[0]
        frame0 = load_trajectory(self.pdb_files[0], top=top, max_frames=1)
        ref_pdb = Path(extra_output_dir) / "reference_frame0.pdb"
        frame0.save_pdb(str(ref_pdb))
        return ref_pdb

    def _write_default_structure_bfactors(self, output_path):
        """
        Color a reference structure (see :meth:`_resolve_reference_pdb`) by
        per-residue H-bond variance across replicates (see
        :meth:`write_structure_bfactors`).
        """
        ref_pdb = self._resolve_reference_pdb(self.output_dir)
        self.write_structure_bfactors(str(ref_pdb), str(output_path), metric="variance")

    def export_run_artifacts(self, output_dir):
        """
        Write per-replicate H-bond tables, Polars-backed occupancy summaries, pivots,
        and a small ``manifest.json`` under ``output_dir/raw_data/``.

        Safe to call again after :meth:`run` with a different path.
        """
        if not self._hbond_dfs:
            raise ValueError("No data. Run .run() first.")
        output_dir = Path(output_dir).expanduser().resolve()
        raw = output_dir / "raw_data"
        raw.mkdir(parents=True, exist_ok=True)
        for i, df in enumerate(self._hbond_dfs):
            stem = Path(self.pdb_files[i]).stem
            df.to_csv(raw / f"{stem}_hbonds.csv", index=False)
            hbond_occupancy_table(df).to_csv(raw / f"{stem}_occupancy.csv", index=False)
        for i, pivot in enumerate(self._pivot_tables):
            stem = Path(self.pdb_files[i]).stem
            pivot.to_csv(raw / f"{stem}_pivot.csv")
        manifest = {
            "input_kind": self.input_kind,
            "pdb_files": [str(p) for p in self.pdb_files],
            "topologies": self.topologies,
            "resid_range": self.resid_range,
            "sub_frames": self.sub_frames,
            "bond_statistics_threshold": self.bond_statistics_threshold,
            "backbone_only": self.backbone_only,
        }
        (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        if self._qc_report is not None:
            (output_dir / "qc_report.json").write_text(
                json.dumps(self._qc_report, indent=2, default=str),
                encoding="utf-8",
            )

    def _compute_bond_statistics(self):
        """Fill self._bond_statistics with mean/std lifetime and break frequency."""
        threshold = self.bond_statistics_threshold
        all_lifetimes = []
        all_breaks = []
        for pivot in self._pivot_tables:
            life, breaks = calculate_bond_statistics(pivot, threshold)
            all_lifetimes.append(life)
            all_breaks.append(breaks)
        lifetimes_df = pd.DataFrame(all_lifetimes)
        breaks_df = pd.DataFrame(all_breaks)
        self._bond_statistics = {
            "mean_lifetimes": lifetimes_df.mean(axis=0),
            "std_lifetimes": lifetimes_df.std(axis=0),
            "mean_break_frequencies": breaks_df.mean(axis=0),
            "std_break_frequencies": breaks_df.std(axis=0),
        }

    def get_hbond_dfs(self):
        """Return list of per-replicate H-bond DataFrames (Bond Label, Frame, ...)."""
        return self._hbond_dfs

    def get_occupancy_tables(self):
        """Return list of per-replicate occupancy summary DataFrames."""
        if not self._hbond_dfs:
            raise ValueError("No data. Run .run() first.")
        return [hbond_occupancy_table(df) for df in self._hbond_dfs]

    def get_qc_report(self):
        """Return QC report dict from latest run (or None if QC was not enabled)."""
        return getattr(self, "_qc_report", None)

    def get_pivot_tables(self):
        """Return list of per-replicate pivot tables (bond label x frame)."""
        return self._pivot_tables

    def get_bond_labels_sorted(self):
        """Return sorted list of all bond labels (union across replicates)."""
        return self._bond_labels_sorted

    def _plot_region_str(self) -> str:
        if self.resid_range is None:
            return "Entire protein"
        start, end = self.resid_range
        return f"Residues {start}-{end}"

    def _filter_bond_labels_for_plot(self, bond_labels):
        """
        Filter bond labels to those fully contained within `self.resid_range`.

        A bond label is included iff BOTH donor and acceptor residue numbers are
        within the range. Donor residue number is the first integer in the bond label.
        """
        if self.resid_range is None:
            return list(bond_labels)
        start, end = self.resid_range
        out = []
        for b in bond_labels:
            d, a = extract_residue_numbers(b)
            if start <= d <= end and start <= a <= end:
                out.append(b)
        return out

    def get_plot_bond_labels_sorted(self):
        """Bond labels (sorted) to use for heatmaps and manifold plots."""
        if not self._bond_labels_sorted:
            return []
        return self._filter_bond_labels_for_plot(self._bond_labels_sorted)

    def get_plot_pivot_tables(self):
        """Reindex each replicate pivot table to the plot bond list."""
        plot_labels = self.get_plot_bond_labels_sorted()
        return [pt.reindex(plot_labels, fill_value=0) for pt in self._pivot_tables]

    @property
    def bond_statistics(self):
        """None or dict with mean_lifetimes, std_lifetimes, mean_break_frequencies, std_break_frequencies."""
        return self._bond_statistics

    def plot_heatmaps(self, save_dir=None, use_legend=False):
        """
        Plot one heatmap per replicate (whole protein or region, depending on resid_range).
        Bond rows = union of bonds found in this set's replicates. For comparison of two sets
        with aligned bond lists, use ComparisonSet.plot_heatmaps_both() instead.
        """
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        plot_pivots = self.get_plot_pivot_tables()
        region_str = self._plot_region_str()
        for i, pivot in enumerate(plot_pivots):
            pdb_path = self.pdb_files[i] if i < len(self.pdb_files) else None
            stem = Path(pdb_path).stem if pdb_path else f"Replicate_{i + 1}"
            title_name = f"{stem} ({region_str})"
            path = f"{save_dir}/{stem}.png" if save_dir else None
            if use_legend:
                plot_heatmap_with_legend(
                    pivot,
                    f"H-Bond Occupation over Time - {title_name}",
                    save_path=path,
                )
            else:
                plot_heatmap(pivot, title_name, save_path=path)

    def plot_top_hbonds(self, threshold=0.5, save_path=None):
        """Bar plot of top H-bonds by average lifetime (presence fraction > threshold)."""
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        plot_top_hbonds(
            self._pivot_tables,
            threshold=threshold,
            sub_frames=self.sub_frames,
            save_path=save_path,
        )

    def plot_bond_statistics(
        self,
        lifetime_title="Average Lifetime of Hydrogen Bonds with Standard Deviation (Threshold = 50%)",
        break_title="Frequency of Breaks of Hydrogen Bonds with Standard Deviation",
        save_lifetime_path=None,
        save_break_path=None,
    ):
        """Bar plots for mean lifetime and break frequency with error bars (requires bond_statistics_threshold)."""
        if self._bond_statistics is None:
            raise ValueError("Bond statistics not computed. Set bond_statistics_threshold and run .run().")
        plot_bond_lifetimes_with_error_bars(
            self._bond_statistics["mean_lifetimes"],
            self._bond_statistics["std_lifetimes"],
            title=lifetime_title,
            save_path=save_lifetime_path,
        )
        plot_break_frequencies_with_error_bars(
            self._bond_statistics["mean_break_frequencies"],
            self._bond_statistics["std_break_frequencies"],
            title=break_title,
            save_path=save_break_path,
        )

    def plot_aggregated_heatmap(self, save_path=None):
        """Average presence of each bond across all replicates."""
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        plot_aggregated_heatmap(
            self.get_plot_pivot_tables(),
            save_path=save_path,
            region_str=self._plot_region_str(),
        )

    def get_fluctuating_bonds(self, quantile=0.9):
        """Bonds with highest variance (fluctuation) across frames. For single-set highlighting."""
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        return fluctuating_bonds(self._pivot_tables, quantile=quantile)

    def _default_replicate_labels(self):
        return [
            self.pdb_files[i] if i < len(self.pdb_files) else f"rep{i + 1}"
            for i in range(len(self._pivot_tables))
        ]

    def run_pca(self, group_labels=None, plot=True, title="PCA"):
        """PCA on this set’s replicates (one point per replicate). Same bond space as two-set case."""
        from phenoms.dimensionality_reduction import aggregate_replicate_data, run_pca as _run_pca, plot_manifold
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        group_labels = group_labels or self._default_replicate_labels()
        data = aggregate_replicate_data(self.get_plot_pivot_tables())
        scores, var_ratio = _run_pca(data, n_components=2)
        if plot:
            plot_manifold(scores, group_labels, title=f"{title} ({self._plot_region_str()})")
        return scores, var_ratio

    def run_tsne(self, group_labels=None, perplexity=2, plot=True, title="t-SNE", random_state=42):
        """t-SNE on this set’s replicates."""
        from phenoms.dimensionality_reduction import aggregate_replicate_data, run_tsne, plot_manifold
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        group_labels = group_labels or self._default_replicate_labels()
        data = aggregate_replicate_data(self.get_plot_pivot_tables())
        scores = run_tsne(data, perplexity=perplexity, random_state=random_state)
        if plot:
            plot_manifold(scores, group_labels, title=f"{title} ({self._plot_region_str()})")
        return scores

    def run_isomap(self, group_labels=None, n_neighbors=5, plot=True, title="Isomap"):
        """Isomap on this set’s replicates."""
        from phenoms.dimensionality_reduction import aggregate_replicate_data, run_isomap, plot_manifold
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        group_labels = group_labels or self._default_replicate_labels()
        data = aggregate_replicate_data(self.get_plot_pivot_tables())
        scores = run_isomap(data, n_neighbors=n_neighbors)
        if plot:
            plot_manifold(scores, group_labels, title=f"{title} ({self._plot_region_str()})")
        return scores

    def run_manifold_suite(self, group_labels=None, perplexity=2, n_neighbors=5, random_state=42, plot=True):
        """
        PCA, t-SNE, and Isomap in one call (single-set replicates as points).
        group_labels: e.g. one label per replicate; default = PDB filenames.
        """
        from phenoms.dimensionality_reduction import run_manifold_suite as _suite
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        group_labels = group_labels or self._default_replicate_labels()
        return _suite(
            self.get_plot_pivot_tables(),
            group_labels,
            perplexity=perplexity,
            n_neighbors=n_neighbors,
            random_state=random_state,
            plot=plot,
        )

    def write_structure_bfactors(self, pdb_path, output_path, metric="variance"):
        """
        Write PDB with B-factors = per-residue metric (for PyMOL/Chimera).
        metric 'variance': highlight fluctuating residues; 'mean': mean occupancy.
        Requires biopython. Single-set only (no comparison).
        """
        from phenoms.comparison import per_donor_metric_from_pivots
        from phenoms.structure import write_pdb_bfactors
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        donor_df = per_donor_metric_from_pivots(
            self._pivot_tables, metric=metric, donor_aggregation="sum"
        )
        write_pdb_bfactors(pdb_path, donor_df, output_path, value_column="Average")

    def export_connectivity_graph_html(
        self,
        output_html_path,
        *,
        top_k_edges: int = 30,
        min_abs_edge_delta: float | None = None,
        directed: bool = False,
        title: str | None = None,
    ) -> None:
        """
        Export an interactive HTML residue connectivity graph from this set's
        own H-bond occupancy (edges sized by occupancy, no second condition
        needed). For a two-condition delta graph, use
        :meth:`~phenoms.ComparisonSet.export_connectivity_graph_html` instead.
        """
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        from phenoms.connectivity import (
            bond_occupancy_from_pivot_tables,
            build_residue_graph_from_bond_deltas,
            export_residue_graph_html,
        )

        unified = self.get_plot_bond_labels_sorted()
        pivots = self.get_plot_pivot_tables()
        if title is None:
            title = "Residue H-bond connectivity (occupancy)"
        bond_values = bond_occupancy_from_pivot_tables(pivots, unified)
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
            show_delta_legend=False,
        )
