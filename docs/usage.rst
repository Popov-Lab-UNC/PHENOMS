Usage
=====

Inputs
------

PHENOMS accepts three input shapes; pick whichever matches what you already have:

* **Multi-frame PDBs** — one PDB per replicate, already protein-centered/fitted.
  A multi-model PDB carries atoms/bonds *and* per-frame coordinates in one
  file, so no separate topology is ever needed here. Simplest and recommended
  when available (see "Simple form" below).
* **Native trajectory + topology** — engine-native trajectories (``.xtc``, ``.trr``,
  ``.dcd``, ``.nc``, ``.mdcrd``, …) loaded directly with a shared or per-replicate
  topology. These trajectory formats carry coordinates only — no atom names,
  bonds, or residues — so a separate topology file is **always required** for
  this path. No intermediate files are written.
* **Engine folders (GROMACS / OpenMM / AMBER)** — point at raw simulation output
  directories; PHENOMS discovers the trajectory/topology pair per replicate
  (same requirement as above), applies PBC imaging/centering/fitting, and
  produces normalized PDBs for you.

All three end up feeding :class:`~phenoms.SimulationSet`; see the matching
sections below for code.

File types by engine
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 15 25 30 30

   * - Engine
     - Trajectory
     - Topology
     - Extra install needed?
   * - GROMACS
     - ``.xtc``, ``.trr``
     - ``.tpr`` (preferred — carries bonds and residue naming exactly as
       simulated)
     - **Yes** — ``pip install "phenoms[gromacs]"``. Without it, ``.tpr``
       raises a clear error telling you to install the extra or convert with
       ``gmx trjconv``.
   * - GROMACS
     - ``.xtc``, ``.trr``
     - ``.gro`` or ``.pdb``
     - No — read natively by MDTraj, core install only.
   * - OpenMM
     - ``.dcd``, ``.xtc``
     - ``.pdb`` or ``.prmtop``
     - No — read natively by MDTraj, core install only.
   * - AMBER
     - ``.nc``, ``.mdcrd``
     - ``.prmtop``, ``.parm7``, ``.prm7``
     - No — read natively by MDTraj, core install only.
   * - Any engine
     - multi-model ``.pdb``
     - *(embedded in the same file)*
     - No.

The ``.tpr`` extra is the only optional-dependency gate on the read path:
every other combination above works with a plain ``pip install phenoms``.
If you have both a ``.tpr`` and a ``.gro``/``.pdb`` for the same replicate and
don't want the extra, point PHENOMS at the ``.gro``/``.pdb`` instead — you
only lose the exact bond/residue naming GROMACS itself used.

.. _default-outputs:

Output location and default artifacts
--------------------------------------

``SimulationSet.run()`` and :class:`~phenoms.ComparisonSet`'s ``.compare()``
both always leave a standard artifact bundle on disk unless you opt out — you
don't need to pass ``output_dir=`` or call any ``get_*``/``plot_*`` method
yourself to get results. ``SimulationSet.run()`` writes:

* ``raw_data/*_hbonds.csv``, ``*_occupancy.csv``, ``*_pivot.csv``, ``manifest.json``
  (and ``qc_report.json`` when ``qc=True``) — see :meth:`~phenoms.SimulationSet.export_run_artifacts`.
* ``plots/`` — a heatmap per replicate, an aggregated heatmap across
  replicates, and (when ``bond_statistics_threshold`` is set) lifetime/break-frequency
  bar plots.
* ``structure_bfactors.pdb`` — a reference structure with B-factors set to
  per-residue H-bond variance across replicates, for coloring in
  PyMOL/Chimera. The reference structure is **frame 0**: the input PDB itself
  for ``pdb_files=`` input, or (for ``trajectories=``/``topology=`` input,
  where no ready PDB exists) frame 0 of the first replicate's trajectory,
  written alongside the other artifacts as ``reference_frame0.pdb``. PHENOMS
  does not auto-trim equilibration on either path — if your simulation
  includes a burn-in period, either pre-trim it first (``phenoms prep
  --start-ps ...``) or treat frame 0 as a coloring reference only; it has no
  effect on H-bond detection, which already ran over the full requested frame
  range.
* ``structure_bfactors_viewer.html`` — a standalone HTML viewer for
  ``structure_bfactors.pdb`` that renders directly in a browser (built with
  NGL, no local software needed) — meant to make the structure result
  immediately visible to someone who doesn't already have PyMOL/Chimera set up.
* ``connectivity_graph.html`` — an interactive residue H-bond network, edges
  sized by occupancy (:meth:`~phenoms.SimulationSet.export_connectivity_graph_html`).

:class:`~phenoms.ComparisonSet`'s ``.compare()`` writes the equivalent bundle
for a comparison:

* ``raw_data/comparison.csv``, ``manifest.json`` — see
  :meth:`~phenoms.ComparisonSet.export_comparison_artifacts`.
* ``plots/difference.png`` and ``plots/heatmaps/`` — see
  :meth:`~phenoms.ComparisonSet.plot_difference` and
  :meth:`~phenoms.ComparisonSet.plot_heatmaps_both`.
* ``structure_bfactors_diff.pdb`` — colored by occupancy difference
  (``label_b`` − ``label_a``), using ``set_a``'s frame 0 as the reference
  structure (same convention as above).
* ``structure_bfactors_diff_viewer.html`` — standalone HTML viewer for the
  differential structure, same as above.
* ``connectivity_graph.html`` — interactive residue Δ-occupancy network
  (``graph_mode="diff"``), see
  :meth:`~phenoms.ComparisonSet.export_connectivity_graph_html`.

Both HTML viewers embed their JS library (NGL / vis-network) from a CDN at
export time so the file is standalone when opened later; if no network is
available at export time, the file still gets written but falls back to
loading that library from the CDN when you open it instead.

Still opt-in (not part of either default bundle, since they need an explicit
choice): single-condition connectivity graphs on a comparison
(``graph_mode="set_a"``/``"set_b"``) and the community-aware variant
(:meth:`~phenoms.ComparisonSet.export_connectivity_community_graph_html`).

Where a bundle goes, by ``output_dir=`` (same for ``SimulationSet`` and
``ComparisonSet``):

* **Omitted (default)** — a fresh, timestamped directory under
  :func:`phenoms.default_output_root` (``./phenom_outputs/`` unless
  ``PHENOMS_OUTPUT_DIR`` is set), so every run gets its own directory and
  nothing is overwritten:

  .. code-block:: bash

     export PHENOMS_OUTPUT_DIR=/path/to/my_phenom_runs

* **A path** — written there instead, e.g. ``output_dir=default_output_root() / "my_run"``.
* **``False``** — disables all default output writing; results stay in memory
  only, via ``sim.get_hbond_dfs()`` / ``get_pivot_tables()`` / ``get_qc_report()``,
  or the individual ``plot_*``/``write_structure_bfactors`` methods called by hand.

Simple form: multi-frame PDBs
-----------------------------

This is the recommended path when you already have protein-centered,
fitted multi-frame PDBs. Backbone N–O analysis is the default
(``backbone_only=True``):

.. code-block:: python

   from phenoms import SimulationSet, default_output_root

   sim = SimulationSet(
       pdb_files=["rep1.pdb", "rep2.pdb"],
       resid_range=(50, 70),       # plot filter only; None = whole protein
       sub_frames=100,             # None = all frames
       backbone_only=True,
       output_dir=default_output_root() / "my_run" / "set_a",
   )
   sim.run()
   # raw_data/*_hbonds.csv, *_occupancy.csv, *_pivot.csv, manifest.json

One-shot helper:

.. code-block:: python

   from phenoms import run_backbone_hbond_analysis

   sim = run_backbone_hbond_analysis(
       ["r1.pdb", "r2.pdb"],
       sub_frames=100,
       plot_heatmaps=False,
       output_dir=default_output_root() / "quick",
   )

Native trajectories
-------------------

Load engine-native files without writing intermediate PDBs. As with the
simple form, ``output_dir`` is optional — see `Output location and default
artifacts`_ for what gets written when it's omitted:

.. code-block:: python

   from phenoms import SimulationSet, default_output_root

   sim = SimulationSet.from_trajectories(
       ["rep1.xtc", "rep2.xtc"],
       topology="system.pdb",          # or topologies=[...]; .gro/.tpr/.prmtop also work
       sub_frames=200,
       backbone_only=True,
       output_dir=default_output_root() / "my_run",
   )
   sim.run()

Equivalent constructor form:

.. code-block:: python

   SimulationSet(
       trajectories=["rep1.xtc", "rep2.xtc"],
       topology="system.pdb",
       output_dir=default_output_root() / "my_run",
   )

GROMACS ``.tpr`` topologies are read via the optional `MDAnalysis
<https://www.mdanalysis.org/>`_ dependency (``pip install "phenoms[gromacs]"``);
without it, supply a ``.gro``/``.pdb`` topology instead.

For PBC imaging / centering / fitting first, use
:func:`phenoms.prepare_set_from_dir` or ``phenoms prep``.

All-bond mode
-------------

Backbone N–O remains the HDX-style default. Opt into all donor/acceptor classes:

.. code-block:: python

   sim = SimulationSet(["rep1.pdb"], backbone_only=False, sub_frames=100)
   sim.run()

Comparison sets
---------------

:class:`~phenoms.ComparisonSet` may be constructed before ``.run()``; methods
that need pivots validate at call time. ``.compare()`` writes the standard
comparison bundle by default (see `Output location and default artifacts`_
above), same as ``SimulationSet.run()``:

.. code-block:: python

   from phenoms import SimulationSet, ComparisonSet

   a = SimulationSet(["a1.pdb", "a2.pdb"], resid_range=(50, 70), sub_frames=100)
   b = SimulationSet(["b1.pdb", "b2.pdb"], resid_range=(50, 70), sub_frames=100)
   cmp = ComparisonSet(a, b, label_a="apo", label_b="holo")
   a.run()
   b.run()
   cmp.compare()
   # raw_data/comparison.csv, plots/difference.png, plots/heatmaps/, structure_bfactors_diff.pdb
   # cmp.export_connectivity_graph_html("network.html")  # opt-in, not part of the default bundle

Engine folders (GROMACS / OpenMM / AMBER)
-----------------------------------------

Directory helpers discover replicates, normalize to protein-only multi-frame
PDBs, then return ready ``SimulationSet`` objects:

.. code-block:: python

   from phenoms import simulation_set_from_dir, comparison_sets_from_dirs

   sim = simulation_set_from_dir(
       input_dir="./wt_set",
       prepared_dir="./prepared/wt",
       frame_dt_ps=1000,
       start_ps=1000,
       end_ps=500000,
   ).run()

   set_a, set_b, comp = comparison_sets_from_dirs(
       dir_a="./wt_set",
       dir_b="./mut_set",
       prepared_dir_a="./prepared/wt",
       prepared_dir_b="./prepared/mut",
       label_a="wt",
       label_b="mut",
   )
   set_a.run()
   set_b.run()
   comp.compare()

Required files per replicate directory — see `File types by engine`_ above for
the full breakdown; summary:

* **GROMACS**: ``.xtc``/``.trr`` + topology. If both a ``.tpr`` and a
  ``.gro``/``.pdb`` are present, ``.tpr`` is preferred (needs
  ``phenoms[gromacs]``); otherwise ``.gro``/``.pdb`` is used with no extra.
* **OpenMM**: ``.dcd``/``.xtc`` + ``.pdb``/``.prmtop``
* **AMBER**: ``.nc``/``.mdcrd`` + ``.prmtop``/``.parm7``

Direct detection helpers
------------------------

Skip the ``SimulationSet`` wrapper when you only need bond tables:

.. code-block:: python

   from phenoms import detect_hbonds, detect_hbonds_with_occupancy

   df_bb = detect_hbonds("trajectory.pdb", backbone_only=True)
   df = detect_hbonds("rep.xtc", top="system.pdb")
   hbonds_df, occ_df = detect_hbonds_with_occupancy(
       "trajectory.pdb",
       output_csv_path="hbond_occupancy.csv",
   )

Quality control
---------------

Optional fail-fast QC inside ``SimulationSet.run``:

.. code-block:: python

   sim = SimulationSet(["rep1.pdb", "rep2.pdb"], sub_frames=250)
   sim.run(
       qc=True,
       mdp_files=["rep1.mdp", "rep2.mdp"],  # optional GROMACS consistency
       qc_fail_on_nonconverged=True,
   )
   report = sim.get_qc_report()

Checks include:

* RMSD window-shift / drift convergence per replicate
* MDP key consistency when ``mdp_files`` are provided

Renumbering & GROMACS export scripts
------------------------------------

For mismatched residue numbering across mutants/constructs:

.. code-block:: bash

   python scripts/renumber_many_to_reference.py \
     --reference ref.pdb --mobile-dir ./exports \
     --mobile-glob "*.pdb" --output-dir ./exports_renumbered

Standalone GROMACS → PDB helper:

.. code-block:: bash

   python scripts/preprocess_gromacs_to_pdb.py \
     --run-dir /path/to/run --out-dir /path/to/exports --frame-dt-ps 1000
