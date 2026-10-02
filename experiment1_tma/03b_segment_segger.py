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
    "libcudart.so.12: cannot open shared object file" unless they are on the path."""
    env = dict(os.environ)
    libs = sorted(glob.glob(os.path.join(repo, ".pixi", "envs", pixi_env, "lib", "python3*",
                                         "site-packages", "nvidia", "*", "lib")))
    if libs:
        env["LD_LIBRARY_PATH"] = ":".join(libs + [env.get("LD_LIBRARY_PATH", "")]).rstrip(":")
        print(f"CUDA libraries of the segger environment put on LD_LIBRARY_PATH ({len(libs)} folders)")
    else:
        print("WARNING: no nvidia/*/lib folders in the segger environment -- is it installed (pixi install -e "
              f"{pixi_env})? If segger then fails on libcudart.so.12, load the cluster's CUDA 12 module.")
    return env


def run(cmd, cwd, env=None):
    """Run a command and stream its output line by line (terminal and Jupyter). Returns the exit code."""
    print("$", " ".join(cmd), flush=True)
    t0 = time.time()
    proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                            bufsize=1)
    for line in proc.stdout:
        print(line, end="", flush=True)
    rc = proc.wait()
    print(f"[exit code {rc}, {(time.time() - t0) / 60:.1f} min]", flush=True)
    return rc


def gpu_report():
    smi = shutil.which("nvidia-smi")
    if smi is None:
        print("WARNING: nvidia-smi not found -- this does not look like a GPU node, and segger "
              "needs a CUDA GPU. Start Jupyter (or this script) inside a GPU job.")
        return False
    out = subprocess.run([smi, "--query-gpu=name,memory.total,driver_version", "--format=csv,noheader"],
                         capture_output=True, text=True)
    print("GPU:", out.stdout.strip() or out.stderr.strip())
    return out.returncode == 0


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
        gpu_report()
        rc = run([pixi, "run", "-e", a.pixi_env, "segger", "segment", "-i", a.xenium_dir, "-o", raw_out,
                  *shlex.split(a.segger_options)], cwd=repo, env=env)
        if rc != 0:
            raise RuntimeError("segger segment failed -- see its output above")
    if not os.path.isfile(seg_parquet):
        raise FileNotFoundError(f"{seg_parquet} was not written. Files in {raw_out}: {sorted(os.listdir(raw_out))}")

    # ---- 2. segger export boundaries (cell polygons -> cell areas)
    if not a.skip_export:
        if a.skip_segment:
            gpu_report()                                # export loads segger's CUDA libraries too
        rc = run([pixi, "run", "-e", a.pixi_env, "segger", "export", "boundaries",
                  "-s", seg_parquet, "-o", export_dir], cwd=repo, env=env)
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
