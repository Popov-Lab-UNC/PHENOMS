Validation data
================

The MD trajectories used to validate PHENOMS on real systems (EGFR WT vs.
del747-749, and PLCg2) are hosted on Zenodo rather than in this repository —
too large for git, and a stable, citable DOI is the point.

**DOI:** `10.5281/zenodo.22879722 <https://doi.org/10.5281/zenodo.22879722>`_
(CC0 1.0 license)

Contents
--------

* **EGFR** — wild-type and del747-749 deletion mutant, 3 replicates each.
  CHARMM force field, GROMACS engine. Protein-centered/fitted multi-frame
  PDB trajectories.
* **PLCg2** — wild-type. AMBER force field, AMBER engine. Topology
  (``.prmtop``) plus a 3500 ns trajectory (``.nc``).

Getting the data
-----------------

.. code-block:: bash

   python scripts/fetch_zenodo_data.py             # both datasets (~8.5GB)
   python scripts/fetch_zenodo_data.py --egfr-only  # ~1GB
   python scripts/fetch_zenodo_data.py --plcg2-only # ~7.5GB

This downloads into ``protein_pdb_exports/`` (EGFR) and ``plcg2_data/``
(PLCg2) at the repository root — the same paths the benchmark scripts
already expect (see ``scripts/renumber_pdb_to_reference.py`` and
``docker/README.md``'s PLCg2 section).

Using it
--------

EGFR WT and the del747-749 mutant are numbered differently, so align them
first:

.. code-block:: bash

   python scripts/renumber_pdb_to_reference.py \\
       --reference protein_pdb_exports/3_19_26_wt_rep_1_protein_centered_fit_nopbc.pdb \\
       --mobile protein_pdb_exports/3_19_26_del747-749_rep1_protein_centered_fit_nopbc.pdb

Then load both sets with :class:`~phenoms.SimulationSet` /
:class:`~phenoms.ComparisonSet` as usual (see :doc:`usage`).

For PLCg2, see the "PLCg2 truncation benchmark" section of
``docker/README.md``.
