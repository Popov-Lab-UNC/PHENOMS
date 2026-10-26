"""
Engine-aware input discovery and trajectory normalization for PHENOMS.

This module converts replicate folders from GROMACS/OpenMM/AMBER layouts into
normalized multi-frame protein-only PDBs consumable by SimulationSet.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import mdtraj as md

from phenoms.io import tpr_to_mdtraj_topology

_TRAJ_EXTS = {".xtc", ".trr", ".dcd", ".nc", ".mdcrd", ".pdb"}


@dataclass(frozen=True)
class ReplicateInput:
    replicate_dir: Path
    engine: str
    trajectory_path: Path
    topology_path: Path | None


def _files_with_exts(directory: Path, exts: set[str]) -> list[Path]:
    return sorted(
        [p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in exts],
        key=lambda p: p.name.lower(),
    )


def _pick_first(paths: Iterable[Path], suffixes: tuple[str, ...]) -> Path | None:
    suffixes = tuple(s.lower() for s in suffixes)
    for p in paths:
        if p.suffix.lower() in suffixes:
            return p
    return None


def detect_replicate_input(replicate_dir: str | Path) -> ReplicateInput:
    """
    Detect simulation engine and enforce required files for one replicate folder.

    Required files:

    - gromacs: trajectory (.xtc/.trr) + topology (.tpr preferred, else .gro/.pdb).
      Reading .tpr requires the optional MDAnalysis dependency (``pip install
      "phenoms[gromacs]"``); without it, use a .gro/.pdb topology instead.
    - openmm: trajectory (.dcd/.xtc) + topology (.pdb/.prmtop)
    - amber: trajectory (.nc/.mdcrd) + topology (.prmtop/.parm7)
    - pdb: one multi-model .pdb trajectory (already normalized input)
    """
    rep = Path(replicate_dir).expanduser().resolve()
    if not rep.is_dir():
        raise ValueError(f"Replicate path is not a directory: {rep}")

    files = [p for p in rep.iterdir() if p.is_file()]
    if not files:
        raise ValueError(f"No files found in replicate directory: {rep}")

    all_traj = _files_with_exts(rep, _TRAJ_EXTS)

    # GROMACS
    gmx_traj = _pick_first(all_traj, (".xtc", ".trr"))
    if gmx_traj is not None:
        gmx_top = _pick_first(files, (".tpr", ".gro", ".pdb"))
        if gmx_top is None:
            raise ValueError(
                f"GROMACS replicate {rep} missing topology file (.tpr/.gro/.pdb) for {gmx_traj.name}."
            )
        return ReplicateInput(rep, "gromacs", gmx_traj, gmx_top)

    # OpenMM
    openmm_traj = _pick_first(all_traj, (".dcd", ".xtc"))
    if openmm_traj is not None:
        openmm_top = _pick_first(files, (".pdb", ".prmtop"))
        if openmm_top is None:
            raise ValueError(
                f"OpenMM replicate {rep} missing topology file (.pdb/.prmtop) for {openmm_traj.name}."
            )
        return ReplicateInput(rep, "openmm", openmm_traj, openmm_top)

    # AMBER
    amber_traj = _pick_first(all_traj, (".nc", ".mdcrd"))
    if amber_traj is not None:
        amber_top = _pick_first(files, (".prmtop", ".parm7"))
        if amber_top is None:
            raise ValueError(
                f"AMBER replicate {rep} missing topology file (.prmtop/.parm7) for {amber_traj.name}."
            )
        return ReplicateInput(rep, "amber", amber_traj, amber_top)

    # PDB-only pathway (already normalized)
    pdbs = [p for p in files if p.suffix.lower() == ".pdb"]
    if len(pdbs) == 1:
        return ReplicateInput(rep, "pdb", pdbs[0], None)

    raise ValueError(
        "Could not detect replicate engine/input files in "
        f"{rep}. Expected GROMACS/OpenMM/AMBER files or a single .pdb."
    )


def discover_replicate_dirs(input_dir: str | Path) -> list[Path]:
    """
    Directory policy:
    - if input_dir itself looks like one replicate, use it
    - otherwise treat each direct child directory as a replicate
    """
    root = Path(input_dir).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f"Input path is not a directory: {root}")

    try:
        detect_replicate_input(root)
        return [root]
    except ValueError:
        pass

    reps = sorted([p for p in root.iterdir() if p.is_dir()], key=lambda p: p.name.lower())
    if not reps:
        raise ValueError(
            f"No replicate directories found in {root}. "
            "Pass either one replicate dir, or a set dir containing replicate subdirs."
        )
    return reps


def _iterload_rep_chunks(rep: ReplicateInput, *, chunk_size: int):
    """
    Yield the replicate's trajectory in chunks (``mdtraj.iterload``) instead
    of loading it in one shot, so the raw file -- usually the largest
    trajectory a user has, since it's read before any frame_dt_ps thinning --
    never has to be held in memory in full.
    """
    if rep.topology_path is None:
        yield from md.iterload(str(rep.trajectory_path), chunk=chunk_size)
        return

    if rep.topology_path.suffix.lower() == ".tpr":
        converted = tpr_to_mdtraj_topology(rep.topology_path, rep.trajectory_path)
        try:
            yield from md.iterload(str(rep.trajectory_path), top=str(converted), chunk=chunk_size)
        finally:
            converted.unlink(missing_ok=True)
        return

    yield from md.iterload(str(rep.trajectory_path), top=str(rep.topology_path), chunk=chunk_size)


def normalize_replicate_to_pdb(
    replicate: ReplicateInput,
    out_pdb_path: str | Path,
    *,
    frame_dt_ps: float | None = 1000.0,
    start_ps: float | None = None,
    end_ps: float | None = None,
    apply_imaging: bool = True,
    center: bool = True,
    fit: bool = True,
    chunk_size: int = 1000,
) -> Path:
    """
    Normalize one replicate:
    image/unwrap (best-effort) -> center -> fit -> protein-only -> save PDB.

    Streams the raw trajectory in chunks instead of loading it whole: each
    chunk is protein-sliced, imaged, and centered as it arrives, and only
    frames that pass ``start_ps``/``end_ps``/``frame_dt_ps`` filtering are
    ever kept in memory. ``fit`` superposes every frame against frame 0 of
    the raw stream as it's processed, matching the reference used when the
    whole trajectory was loaded at once -- since superposition is a
    per-frame operation, doing it per-chunk against that one persisted
    reference frame produces identical coordinates to fitting the full
    trajectory first and thinning afterward.
    """
    protein_idx = None
    reference_frame = None
    kept_frames: list = []
    last_kept_t: float | None = None

    for chunk in _iterload_rep_chunks(replicate, chunk_size=chunk_size):
        if protein_idx is None:
            protein_idx = chunk.topology.select("protein")
        if len(protein_idx) > 0:
            chunk = chunk.atom_slice(protein_idx)

        if apply_imaging:
            try:
                chunk.image_molecules(inplace=True)
            except Exception:
                # Best effort only; fallback keeps the chunk usable.
                pass

        if center:
            chunk.center_coordinates()

        if reference_frame is None:
            reference_frame = chunk.slice(0)

        if fit:
            chunk.superpose(reference_frame, frame=0, atom_indices=list(range(chunk.n_atoms)))

        for i in range(chunk.n_frames):
            t = float(chunk.time[i])
            if start_ps is not None and t < start_ps:
                continue
            if end_ps is not None and t > end_ps:
                continue
            if frame_dt_ps is not None and last_kept_t is not None and t < last_kept_t + float(frame_dt_ps) - 1e-9:
                continue
            kept_frames.append(chunk.slice(i))
            last_kept_t = t

    if not kept_frames:
        raise ValueError("No frames selected after start/end/frame_dt filtering.")

    traj = kept_frames[0] if len(kept_frames) == 1 else md.join(kept_frames)

    out = Path(out_pdb_path).expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    traj.save_pdb(str(out))
    return out


def prepare_set_from_dir(
    input_dir: str | Path,
    prepared_dir: str | Path,
    *,
    frame_dt_ps: float | None = 1000.0,
    start_ps: float | None = None,
    end_ps: float | None = None,
    apply_imaging: bool = True,
    center: bool = True,
    fit: bool = True,
) -> list[str]:
    """
    Prepare one set directory into normalized PDB replicates.
    Returns sorted list of prepared PDB paths.
    """
    reps = discover_replicate_dirs(input_dir)
    out_root = Path(prepared_dir).expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    out_paths: list[str] = []
    for rep_dir in reps:
        rep = detect_replicate_input(rep_dir)
        out_path = out_root / f"{rep_dir.name}_protein_centered_fit_nopbc.pdb"
        normalize_replicate_to_pdb(
            rep,
            out_path,
            frame_dt_ps=frame_dt_ps,
            start_ps=start_ps,
            end_ps=end_ps,
            apply_imaging=apply_imaging,
            center=center,
            fit=fit,
        )
        out_paths.append(str(out_path))
    return sorted(out_paths)
