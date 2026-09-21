PHENOMS documentation
=====================

**PHENOMS** (Python-based Hydrogen-Deuterium Exchange of Molecular Dynamics
Simulations) quantifies local H-bond networks from MD trajectories to support
HDX-MS–style interpretation—signals that global metrics like RMSD/RMSF often miss.

.. image:: ../assets/phenoms.png
   :alt: PHENOMS workflow
   :align: center
   :width: 85%

What you get
------------

* **Baker–Hubbard H-bond detection** with an optional Rust kernel (~17–24× vs
  common Python backends in our 4-CPU Docker benchmark)
* **Backbone N–O mode by default** (HDX-relevant); opt into all donor/acceptor classes
* **Replicate-aware workflows** via :class:`~phenoms.SimulationSet` and
  :class:`~phenoms.ComparisonSet`
* **Inputs**: multi-frame PDBs, native traj+topology, or engine folders
  (GROMACS / OpenMM / AMBER)
* **Exports**: occupancy CSVs, heatmaps, differential protection, connectivity
  graphs, PDB B-factors
* **Optional CLI** (``phenoms``) wrapping the same Python API

At a glance
-----------

.. code-block:: python

   from phenoms import SimulationSet

   sim = SimulationSet(
       pdb_files=["rep1.pdb", "rep2.pdb"],
       sub_frames=100,
       backbone_only=True,  # default
   )
   sim.run()

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   usage
   cli
   outputs
   benchmarks
   data
   api/index
