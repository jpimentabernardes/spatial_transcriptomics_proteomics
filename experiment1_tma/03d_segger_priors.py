#!/usr/bin/env python
"""03d_segger_priors.py -- segger started from two segmentations: 10x Xenium and cellpose.

segger does not draw new cells. It takes an existing segmentation (the
"prior"), learns from the transcripts inside its NUCLEI what each cell looks
like, and then re-assigns every transcript to one of the prior's cells (or to
none). Its output cell id IS the prior's cell id. So running it twice, once
from the 10x cells and once from the cellpose cells, on the same cores, answers:

  * how much does segger change each prior? (transcripts kept in their cell,
    moved to another cell, removed from cells, added from outside cells)
  * which prior is worse? -- the one segger has to correct more, and whose
    contamination (MECR, 03e) drops more after segger
  * do both runs end up at the same answer? (agreement of the two segger
    results vs. agreement of the two priors)
  * for every transcript: its source cell in each prior and its segger cell

Steps (all per slide, on the cores cellpose segmented):
  1. build   two segger input folders in Xenium format, over the SAME transcripts:
               input_xenium    10x cells / nuclei, cut to the cores
               input_cellpose  cellpose cells / nuclei (03a masks) as polygons,
                               transcripts re-labelled with the cellpose cells
  2. segment `segger segment` on each (GPU) [+ `segger export boundaries`]
  3. convert to the common format -> data/segmentation/segger_xenium/<slide>/ and
             segger_cellpose/<slide>/ (03c benchmarks them next to xenium/cellpose)
  4. attribute  one row per transcript: source cell in both priors + segger cell
             of both runs -> transcript_attribution.parquet, plus the summary
             tables 03e_segger_prior_comparison.R plots

Needs 03a run with save_masks (per-core .npz). Re-run 03a with the current
version so the .npz also hold the cellpose NUCLEI; with older .npz the nuclei
are approximated by shrinking each cell by expand_um (a warning says so).

Outputs
  <project>/data/segmentation/segger_priors/<slide>/
      input_xenium/, input_cellpose/       segger inputs (Xenium format)
      raw_xenium/, raw_cellpose/           segger output (segger_segmentation.parquet, export/)
      cores_used.csv, transcript_attribution.parquet, cell_fate_<prior>.csv.gz
  <project>/data/segmentation/segger_xenium/<slide>/, segger_cellpose/<slide>/   common format
  <project>/tables/03d_transcript_fate_<slide>.csv, 03d_gene_fate_<slide>.csv,
                   03d_partition_agreement_<slide>.csv

Example (on a GPU node, see 03b_segment_segger.py for the srun line)
  python 03d_segger_priors.py --slide-id TMA_ORGAN --xenium-dir /path/to/output-XETG... \\
      --project-dir /work/Spatial_TMA/Trial --segger-repo /work/segger
"""

import argparse
import glob
import json
import os
import shlex
import shutil
import sys
import time
from importlib import import_module

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "py"))
from xenium_io import _as_str  # noqa: E402


def as_text(col):
    """Column as strings (cell ids are binary in old XOA, integers in XOA 1.x)."""
    return pc.cast(_as_str(col), pa.string())

PRIORS = ("xenium", "cellpose")
FATES = ["same cell", "moved to another cell", "removed from cells", "added to a cell", "unassigned in both"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slide-id", required=True)
    p.add_argument("--xenium-dir", required=True, help="Xenium output folder (output-XETG...) of this slide")
    p.add_argument("--project-dir", required=True, help="project root (EXP1_PROJECT_DIR)")
    p.add_argument("--segger-repo", default="~/segger", help="folder of the segger git clone (with pixi.toml)")
    p.add_argument("--pixi", default=None, help="path to the pixi program; default: on PATH or ~/.pixi/bin/pixi")
    p.add_argument("--pixi-env", default="cuda121", help="pixi environment of segger")
    p.add_argument("--segger-options", default="",
                   help="extra options for `segger segment` as one string, e.g. '--n-epochs 30'")
    p.add_argument("--priors", nargs="+", default=["xenium", "cellpose"], choices=["xenium", "cellpose"],
                   help="which priors to run (default both)")
    p.add_argument("--cellpose-dir", default=None,
                   help="03a output of this slide; default <project>/data/segmentation/cellpose/<slide>")
    p.add_argument("--cores", nargs="*", default=None,
                   help="only these core_ids (default: every core 03a segmented and saved masks for)")
    p.add_argument("--pad-um", type=float, default=20.0,
                   help="transcripts this far around each core box are included, so cells at the box edge are whole")
    p.add_argument("--build-only", action="store_true", help="only step 1 (no GPU needed)")
    p.add_argument("--skip-build", action="store_true", help="reuse the input folders of an earlier run")
    p.add_argument("--skip-segment", action="store_true", help="reuse segger output; redo convert + attribute")
    p.add_argument("--skip-export", action="store_true", help="no boundary export (cell areas stay NA)")
    p.add_argument("--min-qv", type=float, default=20.0)
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def box_of(x, y, cores, pad=0.0):
    """core_id per point (first matching box, grown by pad), None outside all boxes."""
    out = np.full(len(x), None, dtype=object)
    for c in cores.itertuples(index=False):
        m = ((x >= c.xmin - pad) & (x < c.xmax + pad) & (y >= c.ymin - pad) & (y < c.ymax + pad)
             & (out == None))  # noqa: E711
        out[m] = c.core_id
    return out


def xenium_version(xenium_dir):
    with open(os.path.join(xenium_dir, "experiment.xenium")) as fh:
        v = json.load(fh).get("analysis_sw_version", "xenium-2.0.0")
    try:
        return int(v.split("-")[-1].split(".")[0])
    except ValueError:
        return 2


def read_xenium_cells(xenium_dir):
    for f in ("cells.parquet", "cells.csv.gz"):
        path = os.path.join(xenium_dir, f)
        if os.path.exists(path):
            d = pd.read_parquet(path) if f.endswith(".parquet") else pd.read_csv(path)
            d["cell_id"] = d["cell_id"].astype(str)
            return d[["cell_id", "x_centroid", "y_centroid"]]
    raise FileNotFoundError(f"no cells.parquet / cells.csv.gz in {xenium_dir}")


def read_transcripts_subset(xenium_dir, cores, pad):
    """All columns of transcripts.parquet for the transcripts in the (padded) core boxes, file order kept."""
    dataset = ds.dataset(os.path.join(xenium_dir, "transcripts.parquet"))
    xmin, xmax = cores.xmin.min() - pad, cores.xmax.max() + pad
    ymin, ymax = cores.ymin.min() - pad, cores.ymax.max() + pad
    flt = ((ds.field("x_location") >= xmin) & (ds.field("x_location") < xmax)
           & (ds.field("y_location") >= ymin) & (ds.field("y_location") < ymax))
    parts = []
    for batch in dataset.to_batches(filter=flt, batch_size=5_000_000):
        x = batch.column("x_location").to_numpy(zero_copy_only=False)
        y = batch.column("y_location").to_numpy(zero_copy_only=False)
        keep = box_of(x, y, cores, pad) != None  # noqa: E711
        if keep.any():
            parts.append(pa.Table.from_batches([batch]).filter(pa.array(keep)))
    if not parts:
        raise RuntimeError("no transcripts inside the selected cores")
    tab = pa.concat_tables(parts)
    print(f"{tab.num_rows:,} transcripts in {len(cores)} cores (+{pad:g} um)")
    return tab


def read_boundaries(xenium_dir, fname, keep_ids):
    t = pq.read_table(os.path.join(xenium_dir, fname))
    ids = as_text(t.column("cell_id"))
    t = t.set_column(t.schema.get_field_index("cell_id"), "cell_id", ids)
    return t.filter(pc.is_in(ids, value_set=pa.array(sorted(keep_ids))))


def mask_polygons(labels, keep, offset_rc, px):
    """Outline (closed ring, microns) of each label in `keep` -> dict label -> (x, y)."""
    from skimage.measure import approximate_polygon, find_contours, regionprops
    r0, c0 = offset_rc
    out = {}
    for rp in regionprops(labels):
        if rp.label not in keep:
            continue
        sr, sc = rp.slice
        sub = np.pad(labels[sr, sc] == rp.label, 1).astype(np.float32)
        cs = find_contours(sub, 0.5)
        if not cs:
            continue
        ring = approximate_polygon(max(cs, key=len), tolerance=0.5)
        if len(ring) < 4:
            continue
        if not np.array_equal(ring[0], ring[-1]):
            ring = np.vstack([ring, ring[:1]])
        rows = r0 + sr.start + ring[:, 0] - 1
        cols = c0 + sc.start + ring[:, 1] - 1
        out[rp.label] = ((cols + 0.5) * px, (rows + 0.5) * px)
    return out


def approx_nuclei(masks, expand_px):
    """Nuclei of nuclei_expand cells (old .npz without nuclei): expand_labels grew each nucleus only into
    background, so remove the cell pixels within expand_px of background; touching cells keep their shared side."""
    from scipy import ndimage
    dist_bg = ndimage.distance_transform_edt(masks != 0)
    nuc = masks.copy()
    nuc[dist_bg <= expand_px] = 0
    return nuc


def polygons_table(polys_by_id):
    ids, xs, ys, lab = [], [], [], []
    for i, (cid, (x, y)) in enumerate(polys_by_id.items()):
        ids.append(np.repeat(cid, len(x)))
        xs.append(x)
        ys.append(y)
        lab.append(np.full(len(x), i + 1, dtype=np.int64))
    return pa.table({"cell_id": pa.array(np.concatenate(ids).astype(str)),
                     "vertex_x": pa.array(np.concatenate(xs).astype(np.float32)),
                     "vertex_y": pa.array(np.concatenate(ys).astype(np.float32)),
                     "label_id": pa.array(np.concatenate(lab))})


def write_input(folder, xenium_dir, tx, cell_ids, nuclear, cell_bd, nuc_bd, null_id):
    """A Xenium-format folder segger reads: transcripts with this prior's cell_id / overlaps_nucleus."""
    os.makedirs(folder, exist_ok=True)
    cid = np.where(pd.isna(cell_ids), null_id, cell_ids).astype(str)
    for name, col in (("cell_id", pa.array(cid)),
                      ("overlaps_nucleus", pa.array(np.asarray(nuclear, dtype=np.uint8)))):
        i = tx.schema.get_field_index(name)
        tx = tx.set_column(i, name, col) if i >= 0 else tx.append_column(name, col)
    pq.write_table(tx, os.path.join(folder, "transcripts.parquet"))
    pq.write_table(cell_bd, os.path.join(folder, "cell_boundaries.parquet"))
    pq.write_table(nuc_bd, os.path.join(folder, "nucleus_boundaries.parquet"))
    shutil.copy(os.path.join(xenium_dir, "experiment.xenium"), os.path.join(folder, "experiment.xenium"))
    n_cells = len(set(cell_bd.column("cell_id").to_pylist()))
    print(f"  {folder}: {n_cells:,} cells, {np.mean(~pd.isna(cell_ids)):.1%} of transcripts in cells, "
          f"{np.mean(nuclear):.1%} nuclear")


# ---------------------------------------------------------------------------
# step 1: inputs
# ---------------------------------------------------------------------------
def select_cores(a, cp_dir, cores_csv):
    cores = pd.read_csv(cores_csv)
    cores = cores[(cores["slide_id"] == a.slide_id) & cores["core_id"].notna()].copy()
    cp_cells = pd.read_csv(os.path.join(cp_dir, "cells.csv.gz"))
    have = set(cp_cells["core_id"].dropna())
    masks = {os.path.basename(f)[:-4] for f in glob.glob(os.path.join(cp_dir, "masks", "*.npz"))}
    if not masks:
        raise FileNotFoundError(f"no cellpose masks in {os.path.join(cp_dir, 'masks')} -- run 03a with "
                                "save_masks=True (--save-masks) for this slide")
    no_mask = sorted(have - masks)
    if no_mask:
        print(f"WARNING: {len(no_mask)} cellpose cores without saved masks are left out: {no_mask[:10]}")
    use = have & masks
    if a.cores:
        missing = sorted(set(a.cores) - use)
        if missing:
            print(f"WARNING: no cellpose cells/masks for {missing} -- left out")
        use &= set(a.cores)
    cores = cores[cores["core_id"].isin(use)].reset_index(drop=True)
    if cores.empty:
        raise RuntimeError("no core with cellpose cells + masks among the selected cores")
    print(f"{len(cores)} cores: {', '.join(cores['core_id'])}")
    return cores, cp_cells


def build_inputs(a, work, cp_dir, cores, cp_cells):
    t0 = time.time()
    null_id = "UNASSIGNED" if xenium_version(a.xenium_dir) >= 2 else "-1"
    tx = read_transcripts_subset(a.xenium_dir, cores, a.pad_um)
    x = tx.column("x_location").to_numpy(zero_copy_only=False)
    y = tx.column("y_location").to_numpy(zero_copy_only=False)
    core_id = box_of(x, y, cores)                       # unpadded box: where the transcript is counted
    core_pad = box_of(x, y, cores, a.pad_um)            # padded box: which core's masks cover it

    # --- 10x prior: cells whose centroid lies in a core box (same rule as 03a)
    xc = read_xenium_cells(a.xenium_dir)
    xc = xc[box_of(xc["x_centroid"].to_numpy(), xc["y_centroid"].to_numpy(), cores) != None]  # noqa: E711
    xen_keep = set(xc["cell_id"])
    xid = as_text(tx.column("cell_id")).to_numpy(zero_copy_only=False).astype(object)
    xid[~pd.Series(xid).isin(xen_keep).to_numpy()] = None
    if "overlaps_nucleus" in tx.schema.names:
        xnuc = (tx.column("overlaps_nucleus").to_numpy(zero_copy_only=False).astype(int) == 1) & pd.notna(xid)
    else:
        print("WARNING: transcripts.parquet has no overlaps_nucleus; 10x nuclear transcripts unknown")
        xnuc = np.zeros(len(xid), dtype=bool)

    # --- cellpose prior: 03a cells, masks -> transcript labels + polygons
    params = {}
    if os.path.exists(os.path.join(cp_dir, "run_params.json")):
        params = json.load(open(os.path.join(cp_dir, "run_params.json")))
    cid = np.full(len(x), None, dtype=object)
    cnuc = np.zeros(len(x), dtype=bool)
    cell_polys, nuc_polys = {}, {}
    approximated = False
    for c in cores.itertuples(index=False):
        z = np.load(os.path.join(cp_dir, "masks", f"{c.core_id}.npz"))
        masks, (r0, c0), px = z["masks"], z["offset_rc"], float(z["pixel_size"])
        if "nuclei" in z.files:
            nuclei = z["nuclei"]
        elif params.get("mode", "nuclei_expand") == "nuclei_expand":
            nuclei = approx_nuclei(masks, params.get("expand_um", 5.0) / px)
            approximated = True
        else:
            raise RuntimeError("cellpose masks from cyto mode have no nuclei; segger needs nuclei. Re-run 03a "
                               "(current version saves them) or use mode nuclei_expand")
        kept = cp_cells.loc[cp_cells["core_id"] == c.core_id, "cell_id"]
        keep_lab = {int(s.rsplit("-", 1)[1]) for s in kept}
        idx = np.flatnonzero(core_pad == c.core_id)
        cc = np.floor(x[idx] / px).astype(np.int64) - c0
        rr = np.floor(y[idx] / px).astype(np.int64) - r0
        ok = (cc >= 0) & (cc < masks.shape[1]) & (rr >= 0) & (rr < masks.shape[0])
        idx, cc, rr = idx[ok], cc[ok], rr[ok]
        lab = masks[rr, cc]
        hit = np.isin(lab, list(keep_lab))
        cid[idx[hit]] = [f"{c.core_id}-{v}" for v in lab[hit]]
        cnuc[idx[hit]] = nuclei[rr[hit], cc[hit]] == lab[hit]
        for lab_, xy in mask_polygons(masks, keep_lab, (r0, c0), px).items():
            cell_polys[f"{c.core_id}-{lab_}"] = xy
        for lab_, xy in mask_polygons(nuclei, keep_lab, (r0, c0), px).items():
            nuc_polys[f"{c.core_id}-{lab_}"] = xy
        print(f"  {c.core_id}: {len(keep_lab):,} cellpose cells, {hit.sum():,} transcripts in them", flush=True)
    if approximated:
        print("WARNING: these .npz files have no nuclei -> cellpose nuclei approximated by shrinking each cell by "
              f"{params.get('expand_um', 5.0)} um. Re-run 03a (current version) for the real nuclei.")
    nuc_polys = {k: v for k, v in nuc_polys.items() if k in cell_polys}

    # --- write both inputs over the SAME transcripts (same rows, same order)
    write_input(os.path.join(work, "input_xenium"), a.xenium_dir, tx, xid, xnuc,
                read_boundaries(a.xenium_dir, "cell_boundaries.parquet", xen_keep),
                read_boundaries(a.xenium_dir, "nucleus_boundaries.parquet", xen_keep), null_id)
    write_input(os.path.join(work, "input_cellpose"), a.xenium_dir, tx, cid, cnuc,
                polygons_table(cell_polys), polygons_table(nuc_polys), null_id)

    # --- key table: row_index = row in the input transcripts.parquet = segger's row_index
    key = pd.DataFrame({
        "row_index": np.arange(len(x), dtype=np.int64),
        "x": x, "y": y,
        "feature_name": as_text(tx.column("feature_name")).to_numpy(zero_copy_only=False),
        "qv": tx.column("qv").to_numpy(zero_copy_only=False),
        "core_id": core_id,
        "xenium_cell": xid, "xenium_nuclear": xnuc,
        "cellpose_cell": cid, "cellpose_nuclear": cnuc,
    })
    if "transcript_id" in tx.schema.names:
        key.insert(1, "transcript_id", tx.column("transcript_id").to_numpy(zero_copy_only=False))
    key.to_parquet(os.path.join(work, "transcripts_key.parquet"), index=False)
    print(f"inputs built in {(time.time() - t0) / 60:.1f} min")


# ---------------------------------------------------------------------------
# step 4: attribution and summaries
# ---------------------------------------------------------------------------
def fate(src, seg):
    s, g = pd.notna(src), pd.notna(seg)
    out = np.full(len(src), "unassigned in both", dtype=object)
    out[s & g & (src == seg)] = "same cell"
    out[s & g & (src != seg)] = "moved to another cell"
    out[s & ~g] = "removed from cells"
    out[~s & g] = "added to a cell"
    return out


def agreement(a, b):
    """Over transcripts in a cell in both: share that sit with the bulk of their cell in the other partition."""
    m = pd.notna(a) & pd.notna(b)
    n = int(m.sum())
    if n == 0:
        return n, np.nan, np.nan
    d = pd.DataFrame({"a": a[m], "b": b[m]})
    ab = d.groupby(["a", "b"], sort=False).size()
    a_in_b = ab.groupby(level=0).max().sum() / n      # each a-cell: its largest piece of one b-cell
    b_in_a = ab.groupby(level=1).max().sum() / n
    return n, a_in_b, b_in_a


def attribute(a, work, tables_dir, priors):
    key = pd.read_parquet(os.path.join(work, "transcripts_key.parquet"))
    for p in priors:
        seg = pd.read_parquet(os.path.join(work, f"raw_{p}", "segger_segmentation.parquet"),
                              columns=["row_index", "segger_cell_id", "segger_similarity", "filtered"])
        seg["segger_cell_id"] = seg["segger_cell_id"].astype("string")
        seg.loc[~seg["filtered"].fillna(False).astype(bool), "segger_cell_id"] = pd.NA
        seg = seg.rename(columns={"segger_cell_id": f"segger_{p}_cell",
                                  "segger_similarity": f"segger_{p}_similarity"}).drop(columns="filtered")
        key = key.merge(seg, on="row_index", how="left", indicator=f"in_segger_{p}")
        key[f"in_segger_{p}"] = key[f"in_segger_{p}"] == "both"
        key[f"segger_{p}_cell"] = key[f"segger_{p}_cell"].astype(object).where(key[f"segger_{p}_cell"].notna(), None)
    key.to_parquet(os.path.join(work, "transcript_attribution.parquet"), index=False)
    print(f"transcript attribution -> {os.path.join(work, 'transcript_attribution.parquet')}")

    # transcripts segger scored in every run, inside the core boxes (not the padding)
    use = key["core_id"].notna()
    for p in priors:
        use &= key[f"in_segger_{p}"]
    k = key[use].copy()
    print(f"{len(k):,} transcripts compared (in the core boxes, scored by segger in every run)")

    fates, genes = [], []
    for p in priors:
        k[f"fate_{p}"] = fate(k[f"{p}_cell"].to_numpy(object), k[f"segger_{p}_cell"].to_numpy(object))
        f = k.groupby(["core_id", f"fate_{p}"]).size().rename("n").reset_index().rename(columns={f"fate_{p}": "fate"})
        f.insert(0, "prior", p)
        fates.append(f)
        g = k.groupby(["feature_name", f"fate_{p}"]).size().rename("n").reset_index() \
             .rename(columns={f"fate_{p}": "fate"})
        g.insert(0, "prior", p)
        genes.append(g)

        # per source cell: what happened to its transcripts, and what it gained
        src, sgc = f"{p}_cell", f"segger_{p}_cell"
        own = k[k[src].notna()]
        cell = own.groupby(src).agg(core_id=("core_id", "first"), n_tx_prior=("row_index", "size"),
                                    n_nuclear_prior=(f"{p}_nuclear", "sum"))
        fc = own.groupby([src, f"fate_{p}"]).size().unstack(fill_value=0)
        for col, name in (("same cell", "n_kept"), ("moved to another cell", "n_moved_out"),
                          ("removed from cells", "n_removed")):
            cell[name] = fc[col] if col in fc.columns else 0
        got = k[k[sgc].notna()]
        cell["n_tx_segger"] = got.groupby(sgc).size().reindex(cell.index).fillna(0).astype(int)
        moved_in = got[got[src].notna() & (got[src] != got[sgc])].groupby(sgc).size()
        added = got[got[src].isna()].groupby(sgc).size()
        cell["n_gained_from_cells"] = moved_in.reindex(cell.index).fillna(0).astype(int)
        cell["n_gained_unassigned"] = added.reindex(cell.index).fillna(0).astype(int)
        cell.index.name = "cell_id"
        cell.reset_index().assign(slide_id=a.slide_id, prior=p).to_csv(
            os.path.join(work, f"cell_fate_{p}.csv.gz"), index=False, compression="gzip")

    fates = pd.concat(fates)
    fates.insert(0, "slide_id", a.slide_id)
    fates.to_csv(os.path.join(tables_dir, f"03d_transcript_fate_{a.slide_id}.csv"), index=False)
    genes = pd.concat(genes)
    genes.insert(0, "slide_id", a.slide_id)
    genes.to_csv(os.path.join(tables_dir, f"03d_gene_fate_{a.slide_id}.csv"), index=False)

    # agreement of partitions, per core and over all compared transcripts
    parts = {"xenium": "xenium_cell", "cellpose": "cellpose_cell"}
    parts.update({f"segger_{p}": f"segger_{p}_cell" for p in priors})
    names = list(parts)
    rows = []
    for core, d in list(k.groupby("core_id")) + [("all", k)]:
        for i, p1 in enumerate(names):
            for p2 in names[i + 1:]:
                n, ab, ba = agreement(d[parts[p1]].to_numpy(object), d[parts[p2]].to_numpy(object))
                rows.append({"slide_id": a.slide_id, "core_id": core, "a": p1, "b": p2, "n_tx_both": n,
                             "a_in_b": ab, "b_in_a": ba, "agreement": np.nanmean([ab, ba]) if n else np.nan})
    agr = pd.DataFrame(rows)
    agr.to_csv(os.path.join(tables_dir, f"03d_partition_agreement_{a.slide_id}.csv"), index=False)

    tot = fates.groupby(["prior", "fate"])["n"].sum().unstack(fill_value=0)
    print("\nTranscript fate per prior (% of compared transcripts):")
    print((100 * tot.div(tot.sum(axis=1), axis=0)).round(1).reindex(columns=[f for f in FATES if f in tot.columns])
          .to_string())
    print("\nAgreement between segmentations (all cores; 1 = same cells):")
    print(agr[agr["core_id"] == "all"][["a", "b", "n_tx_both", "agreement"]].round(3).to_string(index=False))
    print(f"\nTables -> {tables_dir}/03d_*_{a.slide_id}.csv ; next: 03c (benchmark) and 03e (comparison figures)")


# ---------------------------------------------------------------------------
def main(argv=None):
    a = parse_args(argv)
    seg = import_module("03b_segment_segger")         # pixi / GPU / run helpers
    repo = os.path.expanduser(a.segger_repo)
    seg_dir = os.path.join(a.project_dir, "data", "segmentation")
    work = os.path.join(seg_dir, "segger_priors", a.slide_id)
    tables_dir = os.path.join(a.project_dir, "tables")
    cp_dir = a.cellpose_dir or os.path.join(seg_dir, "cellpose", a.slide_id)
    cores_csv = os.path.join(a.project_dir, "tables", "02_cores_all.csv")
    for path, what in ((a.xenium_dir, "Xenium folder"), (cores_csv, "core table from 02"),
                       (os.path.join(cp_dir, "cells.csv.gz"), "cellpose output from 03a")):
        if not os.path.exists(path):
            raise FileNotFoundError(f"{what} not found: {path}")
    os.makedirs(work, exist_ok=True)
    os.makedirs(tables_dir, exist_ok=True)

    # ---- 1. inputs
    cores_used = os.path.join(work, "cores_used.csv")
    if a.skip_build or a.skip_segment:
        if not os.path.exists(cores_used):
            raise FileNotFoundError(f"{cores_used} missing -- run once without skip_build")
        print(f"reusing the inputs in {work}")
    else:
        cores, cp_cells = select_cores(a, cp_dir, cores_csv)
        cores.to_csv(cores_used, index=False)
        build_inputs(a, work, cp_dir, cores, cp_cells)
    if a.build_only:
        print("build_only: done. Next, on a GPU node: the same call with skip_build=True")

    # ---- 2. segger on each prior (GPU)
    if not (a.skip_segment or a.build_only):
        seg.gpu_report()
        pixi = seg.find_pixi(a.pixi)
        env = seg.cuda_env(repo, a.pixi_env)
        for p in a.priors:
            raw = os.path.join(work, f"raw_{p}")
            os.makedirs(raw, exist_ok=True)
            print(f"\n===== segger from the {p} segmentation =====")
            rc = seg.run([pixi, "run", "-e", a.pixi_env, "segger", "segment", "-i", os.path.join(work, f"input_{p}"),
                          "-o", raw, *shlex.split(a.segger_options)], cwd=repo, env=env)
            if rc != 0:
                raise RuntimeError(f"segger segment failed for the {p} prior -- see its output above")
            if not a.skip_export:
                os.makedirs(os.path.join(raw, "export"), exist_ok=True)
                rc = seg.run([pixi, "run", "-e", a.pixi_env, "segger", "export", "boundaries",
                              "-s", os.path.join(raw, "segger_segmentation.parquet"),
                              "-o", os.path.join(raw, "export")], cwd=repo, env=env)
                if rc != 0:
                    print("boundary export failed -- continuing without cell areas")

    # ---- 3. common format for 03c (segger_xenium, segger_cellpose)
    conv = import_module("03b_segger_to_common")
    for p in ([] if a.build_only else a.priors):
        raw = os.path.join(work, f"raw_{p}")
        argv_c = ["--slide-id", a.slide_id, "--xenium-dir", a.xenium_dir,
                  "--segger-parquet", os.path.join(raw, "segger_segmentation.parquet"),
                  "--cores-csv", cores_used, "--out-dir", os.path.join(seg_dir, f"segger_{p}", a.slide_id),
                  "--min-qv", str(a.min_qv), "--method-name", f"segger_{p}"]
        bd = os.path.join(raw, "export", "cell_boundaries.parquet")
        if os.path.isfile(bd):
            argv_c += ["--boundaries", bd]
        conv.main(argv_c)

    # ---- 4. transcript attribution + summaries
    if not a.build_only:
        attribute(a, work, tables_dir, a.priors)


if __name__ == "__main__":
    main()
