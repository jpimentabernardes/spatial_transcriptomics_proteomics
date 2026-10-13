#!/usr/bin/env python
"""03b_segment_segger.py -- run segger and convert its output, from Python / Jupyter.

Python version of 03b_segment_segger.sh: same three steps, same outputs.
  1. segger segment            transcript -> cell assignment (needs a CUDA GPU)
  2. segger export boundaries  cell polygons -> cell areas (optional)
  3. conversion to the common format read by R/helpers.R (03b_segger_to_common.py)

segger lives in its own pixi environment (Python 3.11 + CUDA + RAPIDS), so steps
1-2 are run as `pixi run -e cuda121 segger ...` from here, while this script /
notebook itself runs in exp1-imaging. For step 1 the notebook must run on a GPU
node (e.g. Jupyter started inside a GPU job). Step 2 needs the GPU too -- segger
loads its CUDA libraries (CuPy) at start-up, so on a node without GPU even
`segger --help` fails with "libcuda.so.1: cannot open shared object file". Step 3
runs anywhere: skip_segment=True + skip_export=True redoes only the conversion.

`pixi run --frozen` is used: the installed environment is taken as it is, never
re-solved -- compute nodes have no internet, and a re-solve there fails with
"Failed to fetch ... Network unreachable".

Large slides: cuSpatial can fail to build segger's tiling quadtree ("cuSpatial invalid
quadtree ... segger issue #40"). Step 1 therefore runs segger through py/segger_tiling.py,
which falls back to square tiles in that case (tiling="auto"; "square" to use them always,
tile_um to set their side length).

One-time install of segger (the login node is fine; it does not need the GPU):
  curl -fsSL https://pixi.sh/install.sh | sh  &&  source ~/.bashrc
  git clone https://github.com/dpeerlab/segger.git /work/.../segger
  cd /work/.../segger  &&  CONDA_OVERRIDE_CUDA=12.1 pixi install -e cuda121

Outputs
  <project>/data/segmentation/segger_raw/<slide_id>/segger_segmentation.parquet (+ export/)
  <project>/data/segmentation/segger/<slide_id>/   common format for 03c / load_segmentation()

Example
  python 03b_segment_segger.py --slide-id TMA_ORGAN --xenium-dir /path/to/output-XETG... \\
      --project-dir /work/Spatial_TMA --segger-repo /work/segger
"""

import argparse
import glob
import os
import shlex
import shutil
import subprocess
import sys
import time
from importlib import import_module

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slide-id", required=True)
    p.add_argument("--xenium-dir", required=True, help="Xenium output folder (output-XETG...) of this slide")
    p.add_argument("--project-dir", required=True, help="project root (EXP1_PROJECT_DIR)")
    p.add_argument("--segger-repo", default="~/segger", help="folder of the segger git clone (with pixi.toml)")
    p.add_argument("--pixi", default=None, help="path to the pixi program; default: on PATH or ~/.pixi/bin/pixi")
    p.add_argument("--pixi-env", default="cuda121", help="pixi environment of segger")
    p.add_argument("--segger-options", default="",
                   help="extra options for `segger segment` as one string, e.g. '--n-epochs 30' "
                        "(command line: --segger-options='--n-epochs 30'; list: pixi run -e cuda121 segger segment --help)")
    p.add_argument("--skip-segment", action="store_true",
                   help="reuse an existing segger_segmentation.parquet; only redo export + conversion")
    p.add_argument("--skip-export", action="store_true", help="no boundary export (cell areas stay NA)")
    p.add_argument("--cores-csv", default=None, help="default: <project>/tables/02_cores_all.csv")
    p.add_argument("--tiling", default="auto", choices=["auto", "adaptive", "square"],
                   help="how segger cuts the slide into tiles: auto = its quadtree, square tiles if cuSpatial "
                        "fails ('invalid quadtree', large slides); adaptive = quadtree only; square = square only")
    p.add_argument("--tile-um", type=float, default=None,
                   help="side length of square tiles in um (default: from the transcript density)")
    p.add_argument("--min-transcripts", type=int, default=1)
    p.add_argument("--min-qv", type=float, default=20.0)
    return p.parse_args()


def find_pixi(path=None):
    """pixi is often not on Jupyter's PATH (the installer only edits ~/.bashrc)."""
    for c in (path, shutil.which("pixi"), os.path.expanduser("~/.pixi/bin/pixi")):
        if c and os.path.isfile(os.path.expanduser(c)) and os.access(os.path.expanduser(c), os.X_OK):
            return os.path.expanduser(c)
    raise FileNotFoundError("pixi not found. Install it (curl -fsSL https://pixi.sh/install.sh | sh) "
                            "or give its full path as --pixi / pixi=...")


def cuda_env(repo, pixi_env):
    """Environment for segger with the CUDA runtime libraries of its own pixi environment on LD_LIBRARY_PATH.

    PyTorch's pip wheels bring libcudart.so.12, cuBLAS, ... as nvidia-* packages
    (site-packages/nvidia/<lib>/lib), but CuPy does not look there and fails with
    "libcudart.so.12: cannot open shared object file" unless they are on the path; Numba / cuDF
    also need the plain names (libcudart.so), linked in <repo>/.pixi/segger_lib_links/<env>."""
    env = dict(os.environ)
    ## Jupyter sets MPLBACKEND to its inline backend, which segger's own environment does not have
    ## ("'module://matplotlib_inline.backend_inline' is not a valid value for backend")
    env["MPLBACKEND"] = "Agg"
    libs = sorted(glob.glob(os.path.join(repo, ".pixi", "envs", pixi_env, "lib", "python3*",
                                         "site-packages", "nvidia", "*", "lib")))
    if libs:
        ## Numba (used by cuDF) opens the UNVERSIONED names (libcudart.so), the wheels only ship
        ## libcudart.so.12 -> "libcudart.so: cannot open shared object file". Link the plain names
        ## in a folder of our own, put first on the path.
        shim = os.path.join(repo, ".pixi", "segger_lib_links", pixi_env)
        os.makedirs(shim, exist_ok=True)
        n_links = 0
        for d in libs:
            for f in glob.glob(os.path.join(d, "lib*.so.*")):
                base = os.path.basename(f).split(".so")[0] + ".so"
                link = os.path.join(shim, base)
                if not os.path.exists(os.path.join(d, base)) and not os.path.lexists(link):
                    os.symlink(f, link)
                    n_links += 1
        env["LD_LIBRARY_PATH"] = ":".join([shim] + libs + [env.get("LD_LIBRARY_PATH", "")]).rstrip(":")
        print(f"CUDA libraries of the segger environment put on LD_LIBRARY_PATH ({len(libs)} folders"
              f"{', ' + str(n_links) + ' plain-name links added' if n_links else ''})")
    else:
        print("WARNING: no nvidia/*/lib folders in the segger environment -- is it installed (pixi install -e "
              f"{pixi_env})? If segger then fails on libcudart.so.12, load the cluster's CUDA 12 module.")
    return env


def segment_cmd(pixi, pixi_env, env, tiling="auto", tile_um=None):
    """`segger segment` started through py/segger_tiling.py (square-tile fallback for cuSpatial's
    'invalid quadtree' error on large slides). Sets its settings in env; returns the command start."""
    env["SEGGER_TILING"] = tiling
    env["SEGGER_TILE_UM"] = "" if tile_um is None else str(tile_um)
    return [pixi, "run", "--frozen", "-e", pixi_env, "python", os.path.join(HERE, "py", "segger_tiling.py"),
            "segment"]


def run(cmd, cwd, env=None, log=None, every_s=30, tail=25):
    """Run a command; returns its exit code.

    log=None: stream every output line (terminal). log=<file>: write ALL output to that file and
    show only one line every `every_s` seconds plus the last `tail` lines at the end -- segger
    prints thousands of progress lines, which make a notebook page so big that the browser
    reloads it. Follow the full output with: tail -f <log file>."""
    print("$", " ".join(cmd), flush=True)
    if log:
        print(f"  full output -> {log}   (follow it with: tail -f {log})", flush=True)
    t0 = time.time()
    last = 0.0
    recent = []
    proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                            bufsize=1)
    fh = open(log, "a") if log else None
    try:
        for line in proc.stdout:
            if fh is None:
                print(line, end="", flush=True)
                continue
            fh.write(line)
            fh.flush()
            recent = (recent + [line])[-tail:]
            if time.time() - last >= every_s and line.strip():
                last = time.time()
                print(f"  [{(last - t0) / 60:5.1f} min] {line.strip()[:200]}", flush=True)
        rc = proc.wait()
    finally:
        if fh:
            fh.close()
    if fh is not None:
        print("  --- last lines ---\n" + "".join(recent), end="", flush=True)
    print(f"[exit code {rc}, {(time.time() - t0) / 60:.1f} min]", flush=True)
    return rc


NO_GPU_HELP = ("segger needs a GPU, and this job has none it can use. On a GPU node a job only sees GPUs it "
               "asked for, so request one explicitly, e.g.\n"
               "  srun --pty --partition=gpu --gpus-per-node=1 --cpus-per-task=8 --mem=128G --time=12:00:00 /bin/bash\n"
               "(or #SBATCH --partition=gpu + #SBATCH --gpus-per-node=1), then check with: nvidia-smi")


def gpu_report():
    """Stop early, with the fix, when no GPU is visible (instead of segger's 'cudaErrorNoDevice')."""
    smi = shutil.which("nvidia-smi")
    if smi is None:
        raise RuntimeError("nvidia-smi not found -- this is not a GPU node.\n" + NO_GPU_HELP)
    out = subprocess.run([smi, "--query-gpu=name,memory.total,driver_version", "--format=csv,noheader"],
                         capture_output=True, text=True)
    gpus = [ln for ln in out.stdout.strip().splitlines() if ln.strip()]
    if out.returncode != 0 or not gpus:
        raise RuntimeError(f"no GPU visible ({(out.stdout + out.stderr).strip()}).\n" + NO_GPU_HELP)
    print("GPU:", "; ".join(gpus), "| CUDA_VISIBLE_DEVICES =", os.environ.get("CUDA_VISIBLE_DEVICES", "(not set)"))


def check_same_xenium(raw_out, slide_id, xenium_dir):
    """segger_raw/<slide_id>/ remembers which Xenium folder it was made from (xenium_dir.txt): a
    run of another folder under the same slide_id would overwrite that slide's segger output."""
    marker = os.path.join(raw_out, "xenium_dir.txt")
    new = os.path.normpath(os.path.abspath(xenium_dir))
    if os.path.isfile(marker):
        old = open(marker).read().strip()
        if old and os.path.normpath(old) != new:
            raise RuntimeError(
                f"{raw_out} holds segger output of slide_id={slide_id} made from\n    {old}\n"
                f"but this run uses\n    {new}\nCheck slide_id / xenium_dir (slides.csv). If the new folder "
                f"is right, delete {marker} (and the old output) first.")
    with open(marker, "w") as fh:
        fh.write(new + "\n")


def main():
    a = parse_args()
    repo = os.path.expanduser(a.segger_repo)
    raw_out = os.path.join(a.project_dir, "data", "segmentation", "segger_raw", a.slide_id)
    export_dir = os.path.join(raw_out, "export")
    common_out = os.path.join(a.project_dir, "data", "segmentation", "segger", a.slide_id)
    cores_csv = a.cores_csv or os.path.join(a.project_dir, "tables", "02_cores_all.csv")
    seg_parquet = os.path.join(raw_out, "segger_segmentation.parquet")
    for path, what in ((a.xenium_dir, "Xenium folder"), (os.path.join(repo, "pixi.toml"), "segger clone"),
                       (cores_csv, "core table from 02")):
        if not os.path.exists(path):
            raise FileNotFoundError(f"{what} not found: {path}")
    for d in (raw_out, export_dir, common_out):
        os.makedirs(d, exist_ok=True)
    pixi = find_pixi(a.pixi)
    env = cuda_env(repo, a.pixi_env)
    print(f"pixi: {pixi}\nsegger: {repo} (env {a.pixi_env})\noutput: {raw_out}")

    # ---- 1. segger segment (GPU)
    if a.skip_segment:
        print(f"skip_segment: reusing {seg_parquet}")
    else:
        check_same_xenium(raw_out, a.slide_id, a.xenium_dir)
        gpu_report()
        rc = run([*segment_cmd(pixi, a.pixi_env, env, a.tiling, a.tile_um), "-i", a.xenium_dir, "-o", raw_out,
                  *shlex.split(a.segger_options)], cwd=repo, env=env,
                 log=os.path.join(raw_out, "segger_segment.log"))
        if rc != 0:
            raise RuntimeError("segger segment failed -- see its output above")
    if not os.path.isfile(seg_parquet):
        raise FileNotFoundError(f"{seg_parquet} was not written. Files in {raw_out}: {sorted(os.listdir(raw_out))}")

    # ---- 2. segger export boundaries (cell polygons -> cell areas)
    if not a.skip_export:
        if a.skip_segment:
            gpu_report()                                # export loads segger's CUDA libraries too
        rc = run([pixi, "run", "--frozen", "-e", a.pixi_env, "segger", "export", "boundaries",
                  "-s", seg_parquet, "-o", export_dir], cwd=repo, env=env,
                 log=os.path.join(raw_out, "segger_export.log"))
        if rc != 0:
            print("boundary export failed -- continuing without cell areas")
    boundaries = os.path.join(export_dir, "cell_boundaries.parquet")
    if not os.path.isfile(boundaries):
        boundaries = None

    # ---- 3. common format (same code as 03b_segger_to_common.py)
    conv = import_module("03b_segger_to_common")
    argv = ["--slide-id", a.slide_id, "--xenium-dir", a.xenium_dir, "--segger-parquet", seg_parquet,
            "--cores-csv", cores_csv, "--out-dir", common_out,
            "--min-transcripts", str(a.min_transcripts), "--min-qv", str(a.min_qv)]
    if boundaries:
        argv += ["--boundaries", boundaries]
    conv.main(argv)
    print(f"Done: {common_out}\nNext: the other slide, then 03c_segmentation_benchmark.R")


if __name__ == "__main__":
    main()
