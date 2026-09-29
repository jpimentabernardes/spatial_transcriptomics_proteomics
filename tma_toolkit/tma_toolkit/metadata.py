"""Core-level metadata: one row per TMA core.

The core metadata sheet is the single source of truth for what each core IS
(patient, organ, tissue type, ...). Everything else -- selecting cores,
comparing groups, pseudobulk -- reads from it. Keep it in a spreadsheet you
edit, export to CSV, and let `validate_core_metadata` catch mistakes before
they silently split "Colon" and "colon " into two organs.

Columns
-------
Required    core_id, slide_id, core_row, core_col
Recommended patient_id, block_id, organ, site, tissue_type, diagnosis, include
Anything else you add (grade, stage, sex, age, treatment, MSI status, ...) is
kept and becomes selectable/comparable automatically.
"""

from __future__ import annotations

import string
import warnings

import numpy as np
import pandas as pd

REQUIRED = ["core_id", "slide_id", "core_row", "core_col"]
RECOMMENDED = ["patient_id", "block_id", "organ", "site", "tissue_type", "diagnosis",
               "is_control", "replicate", "include", "exclusion_reason", "notes"]
#: columns that describe the PATIENT, so they must agree across all cores of one patient
PATIENT_LEVEL_DEFAULT = ["sex", "age", "treatment"]


def row_labels(n: int) -> list[str]:
    """A..Z, AA, AB, ... (TMA row naming)."""
    out = []
    for i in range(n):
        s, k = "", i
        while True:
            s = string.ascii_uppercase[k % 26] + s
            k = k // 26 - 1
            if k < 0:
                break
        out.append(s)
    return out


def make_template(slide_id: str, n_rows: int, n_cols: int, tma_id: str | None = None,
                  extra_columns: list[str] | None = None, col_width: int = 2) -> pd.DataFrame:
    """Empty metadata sheet for an n_rows x n_cols TMA (e.g. 10 x 20 = 200 cores).

    core_id = <tma_id>_<row><col>, e.g. TMA2_A01 ... TMA2_J20.
    """
    tma_id = tma_id or slide_id
    rows = []
    for r in row_labels(n_rows):
        for c in range(1, n_cols + 1):
            rows.append({"core_id": f"{tma_id}_{r}{c:0{col_width}d}", "slide_id": slide_id,
                         "core_row": r, "core_col": c})
    df = pd.DataFrame(rows)
    for col in RECOMMENDED + PATIENT_LEVEL_DEFAULT + (extra_columns or []):
        df[col] = ""
    df["include"] = True
    df["is_control"] = False
    df["replicate"] = 1
    return df


def _clean_str(s: pd.Series) -> pd.Series:
    return s.map(lambda v: v.strip() if isinstance(v, str) else v).replace({"": np.nan})


def read_core_metadata(path: str, **validate_kw) -> pd.DataFrame:
    """Read CSV / TSV / Excel and validate. Returns the cleaned table indexed by core_id."""
    if path.endswith((".xlsx", ".xls")):
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path, sep=None, engine="python")
    return validate_core_metadata(df, **validate_kw)


def validate_core_metadata(df: pd.DataFrame, allowed_values: dict | None = None,
                           patient_level: list[str] | None = None, fix_case: bool = True,
                           strict: bool = False) -> pd.DataFrame:
    """Checks the sheet and returns a cleaned copy indexed by core_id.

    - required columns present, core_id unique, one core per (slide, row, col)
    - whitespace stripped; empty strings -> NaN
    - categorical values that differ only by case/spaces ("Colon", "colon ") are
      unified (fix_case) and reported
    - optional controlled vocabulary: allowed_values={"organ": [...], "tissue_type": [...]}
    - patient-level columns (sex, age, ...) must agree across a patient's cores
    Problems are collected and raised together (strict) or warned.
    """
    df = df.copy()
    problems = []
    missing = [c for c in REQUIRED if c not in df]
    if missing:
        raise ValueError(f"core metadata lacks required columns {missing}; has {list(df.columns)}")

    for c in df.columns:
        if df[c].dtype == object:
            df[c] = _clean_str(df[c])
    df["core_id"] = df["core_id"].astype(str)
    df["slide_id"] = df["slide_id"].astype(str)
    df["core_row"] = df["core_row"].astype(str)
    df["core_col"] = df["core_col"].map(lambda v: str(int(v)) if isinstance(v, (int, float, np.integer)) and not pd.isna(v) else str(v))

    dup = df["core_id"][df["core_id"].duplicated()].unique()
    if len(dup):
        problems.append(f"duplicated core_id: {list(dup)[:10]}")
    pos_dup = df[df.duplicated(["slide_id", "core_row", "core_col"], keep=False)]
    if len(pos_dup):
        problems.append(f"several cores share one grid position: {pos_dup.core_id.tolist()[:10]}")

    for c in df.columns:
        if c in ("core_id", "notes", "exclusion_reason") or df[c].dtype != object:
            continue
        vals = df[c].dropna().astype(str)
        norm = vals.str.lower().str.replace(r"[\s_\-]+", "_", regex=True)
        groups = vals.groupby(norm).unique()
        variants = {k: list(v) for k, v in groups.items() if len(v) > 1}
        if variants:
            msg = f"column '{c}': spelling variants {variants}"
            if fix_case:
                canon = {x: v[0] for v in variants.values() for x in v}
                df[c] = df[c].map(lambda x: canon.get(x, x))
                msg += " -> unified to the first spelling"
            problems.append(msg)

    for c, allowed in (allowed_values or {}).items():
        if c in df:
            bad = sorted(set(df[c].dropna()) - set(allowed))
            if bad:
                problems.append(f"column '{c}': values not in allowed list: {bad}")

    if "include" in df:
        df["include"] = df["include"].map(
            lambda v: True if pd.isna(v) else str(v).strip().lower() in ("true", "1", "yes", "y"))
    else:
        df["include"] = True
    if "is_control" in df:
        df["is_control"] = df["is_control"].map(lambda v: str(v).strip().lower() in ("true", "1", "yes", "y"))

    if "patient_id" in df:
        for c in (patient_level if patient_level is not None else PATIENT_LEVEL_DEFAULT):
            if c in df:
                n = df.dropna(subset=["patient_id", c]).groupby("patient_id")[c].nunique()
                bad = n[n > 1].index.tolist()
                if bad:
                    problems.append(f"patient-level column '{c}' differs between cores of patients {bad[:10]}")

    for c in df.columns:                              # numbers typed as text -> numbers
        if df[c].dtype == object and c not in ("core_id", "slide_id", "core_row", "core_col", "patient_id"):
            conv = pd.to_numeric(df[c], errors="coerce")
            if conv.notna().sum() == df[c].notna().sum() and df[c].notna().any():
                df[c] = conv

    if problems:
        text = "Core metadata issues:\n  - " + "\n  - ".join(problems)
        if strict:
            raise ValueError(text)
        warnings.warn(text, stacklevel=2)
    return df.set_index("core_id", drop=False).rename_axis(None)


def attach_core_metadata(adata, core_meta: pd.DataFrame, core_key: str = "core_id",
                         columns: list[str] | None = None) -> None:
    """Copy core-level columns onto every cell (adata.obs), in place. Categorical
    dtype for text columns keeps the h5ad small with 200 cores x millions of cells."""
    cols = columns or [c for c in core_meta.columns if c not in ("core_id", "slide_id", "core_row", "core_col")]
    cm = core_meta.set_index("core_id", drop=False) if "core_id" in core_meta else core_meta
    ids = adata.obs[core_key].astype(str)
    for c in cols:
        v = cm[c].reindex(ids.to_numpy()).to_numpy()
        adata.obs[c] = pd.Categorical(v) if cm[c].dtype == object else v
