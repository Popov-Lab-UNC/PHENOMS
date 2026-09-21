# Docker: kernel benchmark in container (demo / reproducibility)

This is **not** part of the published Python wheel; it lives in the GitHub repo so anyone can reproduce the containerized benchmark. Default helper settings run with a **4 CPU** cgroup cap and 4-thread library defaults.

**Platform note:** the Rust-vs-MDTraj comparison is thread-based and reproduces consistently across OSes. The MDAnalysis comparison uses `HydrogenBondAnalysis(..., backend="multiprocessing")`, whose per-call startup cost depends on the OS process-start method (cheap `fork` on Linux vs. expensive `spawn` on macOS/Windows) — running the same script outside this Linux container will very likely show a *larger* (not smaller) speedup over MDAnalysis than the numbers reported here, purely from that platform difference rather than any change in the underlying kernels.

## What it measures

Same script as a local run: **`scripts/benchmark_kernel.py`** — Rust `run_baker_hubbard`, one MDTraj `baker_hubbard` call per replicate, MDAnalysis `HydrogenBondAnalysis.run` on the same first *N* frames. After the run, open **`BENCHMARK_README.txt`** in the output folder for semantics of the three paths (they are not identical chemistry definitions).

## Requirements

- [Docker](https://docs.docker.com/get-docker/) (Desktop on macOS / Engine on Linux).
- Benchmark PDBs under:

  `phenoms/test_data/2_24_25_hbond_heatmaps/`

  with the three replicates expected by the script (see `scripts/benchmark_kernel.py`, `_default_replicates`).  
  If that directory is missing in your clone (large data often gitignored), add those files locally **before** `docker build` so they are copied into the build context.

## Recommended run (250 frames, default 4 CPU)

From the **repository root**:

```bash
export PHENOMS_BENCH_OUT_SUBDIR=docker_4cpu_250f
./scripts/run_benchmark_docker_1cpu.sh --sub-frames 250
```

By default the host directory mounted at `/out` is:

`./phenom_outputs/benchmarks/docker_4cpu_250f/` (under the repo root when you run the script from there)

You get: `kernel_summary.csv`, `kernel_per_replicate.csv`, `BENCHMARK_README.txt`, `benchmark_config.json`.

## Quick smoke (15 frames)

```bash
./scripts/run_benchmark_docker_1cpu.sh
```

## Manual commands

```bash
docker build -f docker/benchmark-one-cpu.Dockerfile -t phenoms-bench:1cpu .
mkdir -p ./phenom_outputs/benchmarks/docker_4cpu_250f
docker run --cpus=4 --rm \
  -e OMP_NUM_THREADS=4 -e OPENBLAS_NUM_THREADS=4 -e MKL_NUM_THREADS=4 \
  -e NUMEXPR_NUM_THREADS=4 -e VECLIB_MAXIMUM_THREADS=4 \
  -v "$(pwd)/phenom_outputs/benchmarks/docker_4cpu_250f:/out" \
  phenoms-bench:1cpu \
  python scripts/benchmark_kernel.py --sub-frames 250 --out-dir /out --rust-threads 4 --mda-workers 4
```

## PLCg2 truncation benchmark (preprocessed trajectory)

Get the data first (DOI: [10.5281/zenodo.22879722](https://doi.org/10.5281/zenodo.22879722), ~7.5GB — see `docs/data.rst`):

```bash
python scripts/fetch_zenodo_data.py --plcg2-only
```

Runs on PLCg2 AMBER inputs in `plcg2_data/` with preprocessing
(protein-only, imaged, centered, fitted), then benchmarks 10 windows by default:
`Full`, `Trunc1_1_1077`, `Trunc2_1_977`, ...

```bash
./scripts/run_benchmark_plcg2_docker.sh
```

Default run parameters:
- `--sub-frames 500`
- `--n-windows 10`
- `--truncation-step 100`
- `--rust-threads 4`
- `--mda-workers 4`

Preprocess on host first (recommended when Docker memory is tight), then benchmark
in container from the preprocessed PDB:

```bash
python scripts/preprocess_plcg2_for_benchmark.py --sub-frames 500
export PHENOMS_PLCG2_PREP_PDB="$(pwd)/phenom_outputs/benchmarks/plcg2_preprocessed/plcg2_preprocessed_500f.pdb"
./scripts/run_benchmark_plcg2_docker.sh
```

- **`--cpus=4`**: hard CPU cap for the container (change as needed).  
- **Image `ENV`**: `OMP_NUM_THREADS=4`, `OPENBLAS_NUM_THREADS=4`, etc.  
- **Script**: also applies the same thread defaults before importing NumPy / MDTraj / MDAnalysis.

To run a strict 1-CPU reference, set:

```bash
export PHENOMS_BENCH_CPUS=1
export PHENOMS_BENCH_OUT_SUBDIR=docker_1cpu_250f
./scripts/run_benchmark_docker_1cpu.sh --sub-frames 250 --rust-threads 1 --mda-workers 1
```

## Custom image tag

```bash
export PHENOMS_BENCH_IMAGE=my-phenoms-bench:latest
./scripts/run_benchmark_docker_1cpu.sh --sub-frames 250
```

## Local run (no Docker)

From repo root, with optional Rust extension (`maturin develop --release`) and `pip install phenoms[benchmark]` for MDAnalysis:

```bash
python scripts/benchmark_kernel.py --sub-frames 250
```

Default output: `$PHENOMS_OUTPUT_DIR/benchmarks/kernel` or `./phenom_outputs/benchmarks/kernel`.
