#!/usr/bin/env python
"""edited_csv_to_geojson.py -- a core grid CSV back into QuPath, to check / correct it again.

Turns <slide>_edited.csv (from export_core_edits.groovy) or <slide>.csv (from 02b) into
a GeoJSON of named circles that QuPath imports (File > Import objects from file, or drag
the file onto the image) on the Xenium morphology image (morphology_focus/ch0000_dapi.ome.tif).

Classes: missing (missing = TRUE), both (has tissue), and check: cores that are not where
their row / column puts them (more than --max-dev x core spacing from a grid fitted to
all cores of their sector -- e.g. two names swapped, a circle dragged onto a neighbour),
plus cores that overlap another. Fix those, set their class back to both / missing, and
export again with export_core_edits.groovy (a core left as 'check' counts as tissue).

Core names like A_1-a, C_3-f (sector _ number - letter) are checked; other names are
converted without the check.

Example
  python qupath/edited_csv_to_geojson.py config/core_grid/TMA_UC_A_edited.csv
  -> config/core_grid/TMA_UC_A_edited_qupath.geojson
"""

import argparse
import json
import os

import numpy as np
import pandas as pd

COLOURS = {"both": [242, 194, 0], "missing": [154, 154, 154], "check": [230, 30, 200]}


def flag_cores(g, max_dev):
    """Cores far from the position their sector / number / letter predicts, and overlapping cores."""
    m = g["core_name"].astype(str).str.extract(r"^([A-Za-z]+)_(\d+)-([A-Za-z])$")
    why = pd.Series("", index=g.index)
    if m[0].isna().all():
        return why
    g = g.assign(sec=m[0], num=pd.to_numeric(m[1]), let=m[2].str.lower().map(lambda c: ord(c) - 96 if isinstance(c, str) else np.nan))
    for sec, s in g.dropna(subset=["sec"]).groupby("sec"):
        if len(s) < 6:
            continue
        X = np.c_[s["num"], s["let"], np.ones(len(s))]
        Y = s[["x_um", "y_um"]].to_numpy(float)
        keep = np.ones(len(s), bool)
        for _ in range(6):                      # robust affine: name -> position
            B, *_ = np.linalg.lstsq(X[keep], Y[keep], rcond=None)
            r = np.linalg.norm(X @ B - Y, axis=1)
            keep = r < max(3 * np.median(r[keep]), 150)
        pitch = min(np.linalg.norm(B[0]), np.linalg.norm(B[1]))
        pred = X @ B
        for k, idx in enumerate(s.index):
            if r[k] > max_dev * pitch:
                near = s.index[int(np.argmin(np.linalg.norm(pred - Y[k], axis=1)))]
                why[idx] = f"{r[k]:.0f} um from its row/column slot (nearest slot: {g.at[near, 'core_name']})"
    P = g[["x_um", "y_um"]].to_numpy(float)
    d = g["diameter_um"].to_numpy(float) if "diameter_um" in g else np.full(len(g), np.nan)
    lim = 0.6 * np.nanmedian(d) if np.isfinite(d).any() else 0
    for i in range(len(g)):
        for j in range(i + 1, len(g)):
            dist = float(np.hypot(*(P[i] - P[j])))
            if dist < lim:
                for a, b in ((i, j), (j, i)):
                    idx = g.index[a]
                    why[idx] = (why[idx] + "; " if why[idx] else "") + f"overlaps {g['core_name'].iloc[b]} ({dist:.0f} um)"
    return why


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", nargs="+", help="<slide>_edited.csv or <slide>.csv (columns core_name, x_um, y_um, diameter_um, missing)")
    p.add_argument("--pixel-size", type=float, default=0.2125, help="um per pixel of the image opened in QuPath "
                   "(Xenium morphology image: 0.2125)")
    p.add_argument("--max-dev", type=float, default=0.45, help="flag cores further than this x core spacing from their slot")
    p.add_argument("--out-dir", default=None, help="default: next to the CSV")
    a = p.parse_args(argv)
    th = np.linspace(0, 2 * np.pi, 48, endpoint=False)
    for f in a.csv:
        g = pd.read_csv(f)
        need = {"core_name", "x_um", "y_um"}
        if not need <= set(g.columns):
            raise SystemExit(f"{f}: needs columns {sorted(need)}; found {list(g.columns)}")
        if "diameter_um" not in g:
            g["diameter_um"] = 1200.0
        miss = g["missing"].astype(str).str.lower().isin(["true", "t", "1", "yes"]) if "missing" in g else False
        g["cls"] = np.where(miss, "missing", "both")
        why = flag_cores(g, a.max_dev)
        g.loc[why != "", "cls"] = "check"
        feats = []
        for r, w in zip(g.itertuples(index=False), why):
            rad = r.diameter_um / 2 / a.pixel_size
            cx, cy = r.x_um / a.pixel_size, r.y_um / a.pixel_size
            ring = [[float(cx + rad * np.cos(t)), float(cy + rad * np.sin(t))] for t in th]
            ring.append(ring[0])
            props = {"objectType": "annotation", "name": str(r.core_name),
                     "classification": {"name": r.cls, "color": COLOURS[r.cls]}, "isLocked": False}
            if w:
                props["description"] = w
            feats.append({"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [ring]},
                          "properties": props})
        out_dir = a.out_dir or os.path.dirname(os.path.abspath(f))
        out = os.path.join(out_dir, os.path.splitext(os.path.basename(f))[0] + "_qupath.geojson")
        with open(out, "w") as fh:
            json.dump({"type": "FeatureCollection", "features": feats}, fh)
        print(f"{os.path.basename(f)}: {len(g)} cores ({int((g.cls == 'missing').sum())} missing, "
              f"{int((g.cls == 'check').sum())} to check) -> {out}")
        for name, w in zip(g["core_name"], why):
            if w:
                print(f"   check {name}: {w}")


if __name__ == "__main__":
    main()
