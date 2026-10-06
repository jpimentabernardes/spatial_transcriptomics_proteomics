#!/usr/bin/env python
"""Generate Jupyter notebooks (.ipynb) from the pipeline's Python scripts.

The .py scripts stay the source: they are what runs as batch jobs on the
cluster. The notebooks are the same code laid out for
interactive work:

  * the module docstring becomes the intro;
  * the command-line options become ONE "Parameters" cell (`--slide-id` ->
    `slide_id`), with the defaults and help texts of the script;
  * every helper function gets its own cell;
  * the body of main() is unrolled into top-level cells, split at the
    script's section comments, so every intermediate variable (images,
    transforms, tables) can be inspected after a cell has run.

Scripts already written as `# %%` cells (examples/) are converted cell by cell.

    python tools/make_notebooks.py            # (re)writes all notebooks

Re-run after changing a script. It OVERWRITES the notebooks, so edits made
in a notebook are lost -- copy a notebook under a new name to keep your own
version.
"""

import ast
import os
import re
import sys
import textwrap

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# script -> folder the notebook goes to (both relative to the repo root)
SCRIPTS = [
    ("experiment1_tma/03a_segment_cellpose.py", "experiment1_tma/notebooks"),
    ("experiment1_tma/03b_segment_segger.py", "experiment1_tma/notebooks"),
    ("experiment1_tma/03b_segger_to_common.py", "experiment1_tma/notebooks"),
    ("experiment1_tma/03d_segger_priors.py", "experiment1_tma/notebooks"),
    ("experiment1_tma/07_register_modalities.py", "experiment1_tma/notebooks"),
    ("tma_toolkit/scripts/make_core_metadata_template.py", "tma_toolkit/notebooks"),
    ("tma_toolkit/scripts/build_tma_object.py", "tma_toolkit/notebooks"),
    ("tma_toolkit/examples/core_selection_walkthrough.py", "tma_toolkit/notebooks"),
]

# packages each notebook needs: import name -> pip name (checked in the first cell, so a
# missing one is reported before anything runs, with the command to install it)
BASE = {"numpy": "numpy", "pandas": "pandas", "scipy": "scipy", "pyarrow": "pyarrow"}
IMG = {**BASE, "tifffile": "tifffile", "zarr": "zarr", "imagecodecs": "imagecodecs", "skimage": "scikit-image",
       "matplotlib": "matplotlib"}
PACKAGES = {
    "experiment1_tma/03a_segment_cellpose.py": {**IMG, "cellpose": "cellpose"},
    "experiment1_tma/03b_segment_segger.py": {**BASE},            # segger itself runs in its pixi env
    "experiment1_tma/03b_segger_to_common.py": {**BASE},          # + geopandas only with boundaries
    "experiment1_tma/03d_segger_priors.py": {**BASE, "skimage": "scikit-image"},
    "experiment1_tma/07_register_modalities.py": {**IMG, "SimpleITK": "SimpleITK"},
    "tma_toolkit/scripts/make_core_metadata_template.py": {**BASE},
    "tma_toolkit/scripts/build_tma_object.py": {**BASE, "anndata": "anndata", "sklearn": "scikit-learn",
                                                "h5py": "h5py", "matplotlib": "matplotlib"},
}

# conda environment each notebook needs: (name, environment file)
IMAGING = ("exp1-imaging", "experiment1_tma/env/py_imaging.yml")
TOOLKIT = ("exp1-imaging", "experiment1_tma/env/py_imaging.yml + pip install -e tma_toolkit")
ENVS = {
    "experiment1_tma/03a_segment_cellpose.py": IMAGING,
    "experiment1_tma/03b_segment_segger.py": IMAGING,
    "experiment1_tma/03b_segger_to_common.py": IMAGING,
    "experiment1_tma/03d_segger_priors.py": IMAGING,
    "experiment1_tma/07_register_modalities.py": IMAGING,
    "tma_toolkit/scripts/make_core_metadata_template.py": TOOLKIT,
    "tma_toolkit/scripts/build_tma_object.py": TOOLKIT,
}

BANNER = re.compile(r"^\s*#\s*([-=])\1{2,}\s*(.*?)\s*[-=]*\s*$")   # "# ----- title" or "# =====" lines


def notebook(cells):
    for i, c in enumerate(cells):
        c.id = f"cell-{i:03d}"            # stable ids: regenerating gives no diff unless code changed
    nb = new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nb.metadata["language_info"] = {"name": "python"}
    return nb


def strip_blank(lines):
    while lines and not lines[0].strip():
        lines = lines[1:]
    while lines and not lines[-1].strip():
        lines = lines[:-1]
    return lines


def split_heading(lines):
    """Leading section banner of a code chunk -> (heading text or None, remaining lines).

    Only banners (comment lines made of dashes / equals signs, optionally with a
    title on the same line, or a title framed by two dash lines) are taken out;
    ordinary comments stay with the code they describe."""
    lines = strip_blank(lines)
    i, titles, banner = 0, [], False
    while i < len(lines) and lines[i].lstrip().startswith("#"):
        m = BANNER.match(lines[i])
        if m:
            banner = True
            if m.group(2):
                titles.append(m.group(2))
        elif banner and len(titles) == 0:
            titles.append(lines[i].strip().lstrip("#").strip())
        else:
            break
        i += 1
    if not banner:
        return None, lines
    return (" ".join(titles) or None), strip_blank(lines[i:])


def segments(nodes, lines, start_line, dedent=0):
    """Source text of each statement, including the comments / blank lines above it."""
    out, prev_end = [], start_line
    for node in nodes:
        seg = lines[prev_end:node.end_lineno]
        if dedent:
            seg = [ln[dedent:] if ln[:dedent].strip() == "" else ln for ln in seg]
        out.append((node, seg))
        prev_end = node.end_lineno
    return out


def is_call_to(node, names):
    if not isinstance(node, ast.Call):
        return False
    f = node.func
    return (isinstance(f, ast.Name) and f.id in names) or (isinstance(f, ast.Attribute) and f.attr in names)


def collect_options(fn, text):
    """add_argument(...) calls -> list of (name, default_src, required, comment)."""
    opts = []
    for node in ast.walk(fn):
        if not (is_call_to(node, {"add_argument"}) and node.args):
            continue
        flag = node.args[0].value
        if not flag.startswith("--"):
            continue
        kw = {k.arg: k.value for k in node.keywords}
        name = kw["dest"].value if "dest" in kw else flag[2:].replace("-", "_")
        action = kw["action"].value if "action" in kw else None
        if "default" in kw:
            default = ast.get_source_segment(text, kw["default"])
        elif action == "store_true":
            default = "False"
        else:
            default = "None"
        required = "required" in kw and kw["required"].value is True
        notes = []
        if required:
            notes.append("REQUIRED")
        if "help" in kw:
            notes.append(kw["help"].value)
        if "choices" in kw:
            notes.append("one of " + ", ".join(repr(c.value) for c in kw["choices"].elts))
        if "nargs" in kw:
            notes.append("a list")
        opts.append((name, default, required, ". ".join(notes)))
    return opts


def parameters_cell(opts, var, script):
    width = max(len(f"    {n}={d},") for n, d, _, _ in opts) + 2
    body = []
    for n, d, _, note in opts:
        line = f"    {n}={d},"
        body.append(line + (" " * (width - len(line)) + "# " + note if note else ""))
    req = [n for n, _, r, _ in opts if r]
    src = ["from types import SimpleNamespace", "",
           f"# Same options as `python {script} --help` (--slide-id -> slide_id).", "",
           "# Options marked REQUIRED have no default: replace their None, e.g. slide_id=\"TMA_ORGAN\".",
           "# Paths are plain strings in quotes.", "",
           f"{var} = SimpleNamespace(", *body, ")"]
    if req:
        src += ["", f"missing = [k for k in {tuple(req)!r} if getattr({var}, k) in (None, '')]",
                "if missing:", "    raise ValueError(f'Fill in these options in the cell above (replace None with a value): {missing}')"]
    return new_code_cell("\n".join(src))


def defaults_cell(opts, var):
    """Separate cell: any option missing from the Parameters cell (e.g. one copied from an
    older version of the notebook) gets the script's default, with a message."""
    body = [f"    {n}={d}," for n, d, _, _ in opts]
    src = ["# Options not set in the Parameters cell (e.g. a cell copied from an older version",
           "# of this notebook) get the script's default. Nothing to edit here.",
           "_defaults = dict(", *body, ")",
           "for _k, _v in _defaults.items():",
           f"    if not hasattr({var}, _k):",
           f"        setattr({var}, _k, _v)",
           "        print(f'{_k} not in the Parameters cell -> default {_v!r}')"]
    req = [n for n, _, r, _ in opts if r]
    if req:
        # checked again here: catches a Parameters cell that was edited but not re-run
        src += [f"_missing = [k for k in {tuple(req)!r} if getattr({var}, k, None) in (None, '')]",
                "if _missing:",
                "    raise ValueError(f'Still empty: {_missing}. Fill them in the Parameters cell, run THAT cell '",
                "                     '(Shift+Enter), then this one.')",
                f"print('Options in use:', {{k: v for k, v in vars({var}).items() if v not in (None, False)}})"]
    return new_code_cell("\n".join(src))


def code_or_heading(cells, seg_lines, level="###"):
    heading, code = split_heading(seg_lines)
    if heading:
        cells.append(new_markdown_cell(f"{level} {heading[0].upper() + heading[1:]}"))
    if code:
        cells.append(new_code_cell("\n".join(code)))


def unroll_main(fn, lines, text):
    """Body of main() -> list of cells at top level. Returns (cells, args variable name)."""
    todo = list(fn.body)                       # look for `return` in main itself, not in nested defs
    while todo:
        node = todo.pop()
        if isinstance(node, ast.Return):
            sys.exit(f"main() has a return statement (line {node.lineno}); rewrite it as if/else first")
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda, ast.ClassDef)):
            todo.extend(ast.iter_child_nodes(node))
    var, kept = "args", []
    for node, seg in segments(fn.body, lines, fn.body[0].lineno - 1, dedent=4):
        if isinstance(node, ast.Assign) and is_call_to(node.value, {"parse_args"}):
            var = node.targets[0].id
            continue
        if isinstance(node, ast.Assign) and is_call_to(node.value, {"ArgumentParser"}):
            continue
        if isinstance(node, ast.Expr) and is_call_to(node.value, {"add_argument"}):
            continue
        kept.append(seg)

    has_banners = any(split_heading(seg)[0] is not None or
                      (strip_blank(seg) and BANNER.match(strip_blank(seg)[0])) for seg in kept)
    groups = []
    for seg in kept:
        if has_banners:
            new = split_heading(seg)[1] != strip_blank(seg) or not groups
        else:
            new = not groups or (seg and not seg[0].strip())      # blank line above = new group
        if new:
            groups.append(list(seg))
        else:
            groups[-1] += seg
    if not has_banners:   # merge very short groups into the next one
        merged = []
        for g in groups:
            if merged and len(strip_blank(merged[-1])) < 3:
                merged[-1] += g
            else:
                merged.append(g)
        groups = merged

    cells = []
    for g in groups:
        src = "\n".join(g)
        src = re.sub(r"\bsys\.exit\(", "raise RuntimeError(", src)       # stop with a readable error
        code_or_heading(cells, src.split("\n"))
    return cells, var



def script_dir_code(base, rel_dir, needs_py):
    """First cell: find the folder holding the script (and py/xenium_io.py), wherever the notebook was put."""
    req = [base] + (["py/xenium_io.py"] if needs_py else [])
    return (
        f"# Folder with {base}{' and py/xenium_io.py' if needs_py else ''}. Normally this notebook sits in\n"
        f"# notebooks/ next to the scripts ({rel_dir!r}); it is also found if the notebook is in the script\n"
        "# folder itself or deeper. If not found, set SCRIPT_DIR yourself, e.g.\n"
        "# SCRIPT_DIR = '/work_ikmb/<user>/Spatial_TMA/Scripts_tma'\n"
        "SCRIPT_DIR = None\n"
        f"_required = {req!r}\n"
        f"_cands = [os.path.abspath({rel_dir!r})] + [os.path.abspath(os.path.join(*['..'] * k) if k else '.') for k in range(4)]\n"
        "if SCRIPT_DIR is None:\n"
        "    SCRIPT_DIR = next((d for d in _cands if all(os.path.exists(os.path.join(d, f)) for f in _required)), None)\n"
        "if SCRIPT_DIR is None or not all(os.path.exists(os.path.join(SCRIPT_DIR, f)) for f in _required):\n"
        "    _near = next((d for d in _cands if os.path.exists(os.path.join(d, _required[0]))), None)\n"
        "    raise FileNotFoundError(\n"
        "        (f'Found {_required[0]} in {_near}, but not {_required[1:]} there -- copy the py/ folder of the '\n"
        "         'repository next to the scripts.' if _near else\n"
        "         f'Could not find {_required[0]} in {sorted(set(_cands))}. Put this notebook in the notebooks/ '\n"
        "         'folder next to the scripts, or set SCRIPT_DIR at the top of this cell.')\n"
        "        + f'\\n(this notebook runs in {os.getcwd()})')\n"
        "print('Scripts folder:', SCRIPT_DIR)")

def convert_script(src, nb_dir):
    path = os.path.join(REPO, src)
    text = open(path).read()
    lines = text.split("\n")
    tree = ast.parse(text)
    base = os.path.basename(src)
    rel_dir = os.path.relpath(os.path.dirname(path), os.path.join(REPO, nb_dir))

    doc = ast.get_docstring(tree, clean=True) or base
    first, _, rest = doc.partition("\n")
    title = first.replace(f"{base} --", "").replace(base, "").strip(" -") or base
    cells = [new_markdown_cell(f"# {base[:-3]}: {title}\n\n```text\n{rest.strip()}\n```")]
    cells.append(new_markdown_cell(
        f"**Notebook version of `{os.path.relpath(path, REPO)}`.** Same code, laid out for interactive work:\n\n"
        "1. Set the options in the **Parameters** cell (the script's command-line options: "
        "`--slide-id` becomes `slide_id`).\n"
        "2. Run the cells in order. After each cell you can inspect its variables, tables and images.\n\n"
        "Use the `.py` script for batch jobs on the cluster; both produce the same outputs. "
        "This notebook is generated by `tools/make_notebooks.py`, which overwrites it, so save your "
        "own edits under another name."))
    env = ENVS.get(src)
    check = ""
    if env:
        e = env[0]
        check = ("\n\n# Which Python runs this notebook? It must be the conda environment, not the system one.\n"
                 "import sys\nprint('Python:', sys.executable)\n"
                 "try:\n    import numpy  # noqa: F401\n"
                 "except ModuleNotFoundError:\n"
                 "    # find the environment's own python, so the fix below names the exact path\n"
                 f"    cands = [os.path.expanduser('~/.conda/envs/{e}/bin/python'),\n"
                 f"             os.path.join(sys.prefix, 'envs', '{e}', 'bin', 'python')]\n"
                 f"    cands += [os.path.join(d, '{e}', 'bin', 'python')\n"
                 "              for d in os.environ.get('CONDA_ENVS_PATH', '').split(os.pathsep) if d]\n"
                 "    found = [c for c in cands if os.path.exists(c)]\n"
                 f"    py = found[0] if found else '\"$CONDA_PREFIX/bin/python\"'\n"
                 "    raise ModuleNotFoundError(\n"
                 f"        'This kernel runs ' + sys.executable + ', not the {e} environment.\\n'\n"
                 f"        'If you already picked \"{e}\" under Kernel > Change kernel, that kernel was registered\\n'\n"
                 "        'with the wrong Python. Fix it once in a terminal (re-registers it with the right one):\\n'\n"
                 f"        '  conda activate {e}\\n'\n"
                 "        f'  {py} -m pip install ipykernel\\n'\n"
                 f"        f'  {{py}} -m ipykernel install --user --name {e} --display-name {e}\\n'\n"
                 f"        'then reload the Jupyter page, Kernel > Change kernel > {e}, run from the top.\\n'\n"
                 f"        '(environment file: {env[1]})')")
    pkgs = PACKAGES.get(src)
    if pkgs:
        check += ("\n\n# Packages this step needs (import name: pip name). Anything missing is installed INTO THIS\n"
                  "# kernel with %pip in a notebook cell -- not with pip in a terminal, which may be another Python.\n"
                  "import importlib.util\n"
                  f"needed = {pkgs!r}\n"
                  "missing = [pip for mod, pip in needed.items() if importlib.util.find_spec(mod) is None]\n"
                  "if missing:\n"
                  "    raise ModuleNotFoundError(\n"
                  "        f'Missing in this kernel: {\", \".join(missing)}\\n'\n"
                  "        f'Run this in a new cell:  %pip install {\" \".join(missing)}\\n'\n"
                  "        'then Kernel > Restart, and run the notebook again from the top.')\n"
                  "print('All packages found.')")
        if "cellpose" in pkgs:
            check += ("\n\n# cellpose runs on PyTorch. On a GPU node check that it sees the GPU (else set gpu=False,\n"
                      "# or install a CUDA build of torch first, see pytorch.org):\n"
                      "try:\n    import torch\n    print('torch', torch.__version__, '| GPU available:', torch.cuda.is_available())\n"
                      "except ImportError:\n    print('torch not found -- cellpose needs it:  %pip install torch')")
    cells.append(new_code_cell(
        "import os\n\n"
        + script_dir_code(base, rel_dir, "xenium_io" in text or "03b_segger_to_common" in text) + check))

    fns = {n.name: n for n in tree.body if isinstance(n, ast.FunctionDef)}
    main = fns.get("main")
    opt_fn = fns.get("parse_args") or main
    opts = collect_options(opt_fn, text) if opt_fn else []

    body_cells, block, helpers_started, params_done = [], [], False, False
    main_cells, var = unroll_main(main, lines, text) if main else ([], "args")

    def flush():
        nonlocal block, params_done
        if strip_blank(block):
            code_or_heading(body_cells, block)
            if not params_done and opts:
                body_cells.append(new_markdown_cell("## Parameters"))
                body_cells.append(parameters_cell(opts, var, os.path.relpath(path, REPO)))
                body_cells.append(defaults_cell(opts, var))
                params_done = True
        block = []

    for node, seg in segments(tree.body, lines, 0):
        if isinstance(node, ast.Expr) and isinstance(getattr(node, "value", None), ast.Constant) \
                and node is tree.body[0]:
            continue                                                     # module docstring
        if isinstance(node, ast.If) and "__name__" in ast.get_source_segment(text, node.test):
            continue
        if isinstance(node, ast.FunctionDef) and node.name in ("parse_args", "main"):
            continue
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            flush()
            if not helpers_started:
                body_cells.append(new_markdown_cell("## Helper functions"))
                helpers_started = True
            code_or_heading(body_cells, seg)
        else:
            block += seg
    flush()

    cells += body_cells
    if main_cells:
        cells.append(new_markdown_cell("## Run"))
        cells += main_cells
    for c in cells:
        if c.cell_type == "code":
            s = c.source.replace("os.path.dirname(os.path.abspath(__file__))", "SCRIPT_DIR")
            s = s.replace("os.path.abspath(__file__)", f"os.path.join(SCRIPT_DIR, {base!r})")
            # scripts force the file-only "Agg" backend; in Jupyter keep the inline one
            c.source = re.sub(r"\n[ \t]*matplotlib\.use\(\"Agg\"\)[^\n]*", "", s)
    return notebook(cells)


def convert_percent(src):
    """A script already split with `# %%` markers -> one notebook cell per block."""
    text = open(os.path.join(REPO, src)).read()
    cells = []
    for chunk in re.split(r"^# %%", text, flags=re.M)[1:]:
        header, _, body = chunk.partition("\n")
        header = header.strip()
        if header.startswith("[markdown]"):
            md = [re.sub(r"^# ?", "", ln) for ln in body.strip().split("\n")]
            cells.append(new_markdown_cell("\n".join(md)))
        else:
            code = strip_blank(body.split("\n"))
            if header:
                code = [f"# {header}"] + code
            cells.append(new_code_cell("\n".join(code)))
    return notebook(cells)


def main():
    for src, nb_dir in SCRIPTS:
        text = open(os.path.join(REPO, src)).read()
        nb = convert_percent(src) if re.search(r"^# %%", text, flags=re.M) else convert_script(src, nb_dir)
        if src in ENVS:
            # open with the conda kernel registered by `python -m ipykernel install --name <env>`;
            # "Python 3" would be whatever Python started Jupyter (often the base install, no numpy)
            name = ENVS[src][0]
            nb.metadata["kernelspec"] = {"name": name, "display_name": name, "language": "python"}
        nbformat.validate(nb)
        out = os.path.join(REPO, nb_dir, os.path.basename(src)[:-3] + ".ipynb")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        nbformat.write(nb, out)
        print(f"{src} -> {os.path.relpath(out, REPO)} ({len(nb.cells)} cells)")


if __name__ == "__main__":
    main()
