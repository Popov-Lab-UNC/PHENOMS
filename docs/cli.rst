Command-line interface
======================

The Python API is primary. After ``pip install -e .``, a Click CLI mirrors common
workflows:

.. code-block:: bash

   phenoms --help
   python -m phenoms --help

``phenoms prep``
----------------

Normalize GROMACS / OpenMM / AMBER replicate folders to multi-frame PDBs
(imaging / center / fit, protein-only):

.. code-block:: bash

   phenoms prep --input-dir ./wt_set --prepared-dir ./prepared/wt \
     --frame-dt-ps 1000 --start-ps 1000 --end-ps 500000

Flags: ``--no-imaging``, ``--no-center``, ``--no-fit``.

The raw trajectory (before ``--frame-dt-ps`` thinning) is read in streamed
chunks rather than loaded whole, since it's usually the largest file a user
has — this avoids OOM on long production runs regardless of how small the
final, thinned output ends up being.

``phenoms run``
---------------

Analyze PDB replicates **or** native trajectories (not both in one invocation):

.. code-block:: bash

   # PDB form
   phenoms run --pdb rep1.pdb --pdb rep2.pdb \
     --sub-frames 100 --output-dir ./out

   # Native traj + shared topology
   phenoms run --traj rep1.xtc --traj rep2.xtc --topology top.pdb

   # All-bond mode (default is backbone N–O only)
   phenoms run --pdb rep1.pdb --all-bonds

   # Optional QC
   phenoms run --pdb rep1.pdb --pdb rep2.pdb --qc

``phenoms compare``
-------------------

Prep two class directories, run analysis, and export comparison artifacts:

.. code-block:: bash

   phenoms compare \
     --dir-a ./wt --dir-b ./mut \
     --prepared-dir-a ./prepared/wt --prepared-dir-b ./prepared/mut \
     --label-a wt --label-b mut \
     --output-dir ./phenom_outputs/cli_compare
