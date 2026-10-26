Outputs
=======

Core analysis outputs
---------------------

``SimulationSet.run()`` writes this bundle by default — see :ref:`default-outputs`
below for exactly when and where. Under ``output_dir/raw_data/``:

* **``*_hbonds.csv``** — frame-level donor/hydrogen/acceptor indices, frame, bond label
* **``*_occupancy.csv``** — present-frame counts and occupancy fractions
* **``*_pivot.csv``** — bond-label × frame matrices for heatmaps and comparisons
* **``manifest.json``** — run metadata (inputs, ``sub_frames``, ``resid_range``,
  ``backbone_only``, …)
* **``qc_report.json``** — present when QC was enabled

And directly under ``output_dir``:

* **``plots/``** — per-replicate heatmaps, an aggregated heatmap, and (when
  ``bond_statistics_threshold`` is set) lifetime/break-frequency bar plots
* **``structure_bfactors.pdb``** — frame-0 reference structure colored by
  per-residue H-bond variance across replicates (see :ref:`default-outputs` for
  where "frame 0" comes from on each input path)
* **``structure_bfactors_viewer.html``** — standalone, no-install 3D viewer for
  the same structure (opens directly in a browser; built with NGL)
* **``connectivity_graph.html``** — interactive residue H-bond network, edges
  sized by occupancy (:meth:`~phenoms.SimulationSet.export_connectivity_graph_html`)

Differential / comparison outputs
---------------------------------

``ComparisonSet.compare()`` writes this bundle by default too (see
:ref:`default-outputs`) — no need to call ``export_comparison_artifacts`` or
``plot_difference`` yourself:

* **``raw_data/comparison.csv``**, **``manifest.json``** — per-donor occupancy
  in each condition, raw and clipped difference (clipped to ``[-1, 1]`` for
  visualization / B-factors); see :meth:`~phenoms.ComparisonSet.export_comparison_artifacts`
* **``plots/difference.png``** — :meth:`~phenoms.ComparisonSet.plot_difference`
  with manual or autocorrelation-based thresholds
* **``plots/heatmaps/``** — aligned per-replicate heatmaps for both sets, via
  :meth:`~phenoms.ComparisonSet.plot_heatmaps_both`
* **``structure_bfactors_diff.pdb``** — differential B-factors (``label_b`` −
  ``label_a``) on ``set_a``'s frame-0 reference structure
* **``structure_bfactors_diff_viewer.html``** — standalone, no-install 3D
  viewer for the same differential structure (built with NGL)
* **``connectivity_graph.html``** — interactive residue Δ-occupancy network,
  ``graph_mode="diff"`` (:meth:`~phenoms.ComparisonSet.export_connectivity_graph_html`)

Opt-in only (not part of the default bundle, since they need an explicit
choice): single-condition connectivity graphs (``graph_mode="set_a"``/``"set_b"``)
and the community-aware variant
(:meth:`~phenoms.ComparisonSet.export_connectivity_community_graph_html`).

Default location
----------------

A fresh, timestamped directory under ``./phenom_outputs/`` (or
``$PHENOMS_OUTPUT_DIR``) per run — see :ref:`default-outputs` for the exact
naming and how to override it with ``output_dir=`` / ``--output-dir``, or
disable default writing entirely with ``output_dir=False``.
