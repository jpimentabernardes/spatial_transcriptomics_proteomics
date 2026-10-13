"""Run the segger command line with fallbacks for cuSpatial's "invalid quadtree" error.

On large slides (tens of millions of transcripts) cuSpatial can fail to build the quadtrees segger
uses, at every size segger tries:
    RuntimeError: cuSpatial invalid quadtree after sizes [...] (see segger issue #40)
segger builds them in two places, and this file patches both (segger itself is not changed):
  1. cutting the slide into tiles (QuadTreeTiling). segger also has a plain square-grid tiling,
     but its command line does not expose it: the quadtree is tried first and, when it fails,
     square tiles are used. Empty squares are dropped and, unless given, the side length is
     chosen so that the densest squares hold about as many nodes as a quadtree tile would
     (--max-nodes-per-tile, default 50k).
  2. finding which tile / region each transcript lies in (point-in-polygon join). When the
     quadtree over all points fails, the points are cut into vertical stripes (4, 16, 64, 256)
     and each stripe gets its own, smaller quadtree. Same result, just in parts.

Started by 03b/03d inside segger's pixi environment, with the same arguments as `segger`:
    pixi run --frozen -e cuda121 python <this file> segment -i <xenium_dir> -o <out>
Settings (environment variables):
    SEGGER_TILING   auto (quadtree, square tiles if it fails; default) | square (square tiles
                    only) | adaptive (plain segger, no fallbacks at all)
    SEGGER_TILE_UM  side length of the square tiles in um (default: from the transcript density)
"""
import math
import os
import sys

## this folder (py/) holds the pipeline's own modules; keep it off the import path of segger
_here = os.path.dirname(os.path.abspath(__file__))
sys.path = [p for p in sys.path if os.path.abspath(p or ".") != _here]

import geopandas as gpd
import numpy as np
from shapely import box

import segger.data.data_module as data_module
import segger.geometry.query as query
from segger.data.tiling import QuadTreeTiling, Tiling

MODE = os.environ.get("SEGGER_TILING", "auto").strip().lower()
TILE_UM = os.environ.get("SEGGER_TILE_UM", "").strip()


def say(msg):
    print(f"[segger_tiling] {msg}", file=sys.stderr, flush=True)


class GridTiling(Tiling):
    """Square tiles over the extent of the positions; tiles without any position are dropped."""

    def __init__(self, xy, side_length):
        x0, y0 = xy.min(0)
        x1, y1 = xy.max(0)
        nx = max(1, math.ceil((x1 - x0) / side_length))
        ny = max(1, math.ceil((y1 - y0) / side_length))
        ex = np.linspace(x0, max(x1, x0 + 1e-3), nx + 1)
        ey = np.linspace(y0, max(y1, y0 + 1e-3), ny + 1)
        counts, _, _ = np.histogram2d(xy[:, 0], xy[:, 1], bins=[ex, ey])
        ii, jj = np.nonzero(counts)
        self._tiles = gpd.GeoSeries([box(ex[i], ey[j], ex[i + 1], ey[j + 1]) for i, j in zip(ii, jj)])
        self.side_length = float(ex[1] - ex[0])
        self.counts = counts[ii, jj]
        super().__init__()

    @property
    def tiles(self):
        return self._tiles


def side_from_density(xy, max_tile_size, bin_um=50.0):
    """Side length (um) at which the densest squares (99.9th percentile) hold ~max_tile_size nodes."""
    span = xy.max(0) - xy.min(0)
    nb = np.maximum(1, np.ceil(span / bin_um)).astype(int)
    h, _, _ = np.histogram2d(xy[:, 0], xy[:, 1], bins=nb)
    area = max(span[0] / nb[0], 1e-6) * max(span[1] / nb[1], 1e-6)
    dense = np.percentile(h[h > 0], 99.9) / area           # nodes per um^2
    return float(np.clip(math.sqrt(max_tile_size / dense), 100.0, 1000.0))


def square_tiling(positions, max_tile_size):
    xy = positions.detach().cpu().numpy().astype(np.float64)
    side = float(TILE_UM) if TILE_UM else side_from_density(xy, max_tile_size)
    t = GridTiling(xy, side)
    say(f"square tiling: {len(t.tiles)} tiles of {t.side_length:.0f} um, nodes per tile: "
        f"median {np.median(t.counts):.0f}, max {t.counts.max():.0f} (quadtree target {max_tile_size})")
    return t


def tiling_with_fallback(positions, max_tile_size):
    if MODE == "square":
        return square_tiling(positions, max_tile_size)
    try:
        return QuadTreeTiling(positions=positions, max_tile_size=max_tile_size)
    except RuntimeError as e:
        if MODE == "adaptive" or "quadtree" not in str(e):
            raise
        say(f"quadtree tiling failed ({e}) -> using square tiles instead")
        return square_tiling(positions, max_tile_size)


_contains = query._points_in_polygons_contains


def contains_in_stripes(points, polygons, max_size=None, batches=None):
    """segger's point-in-polygon join; if cuSpatial cannot build the quadtree over all points, the
    points are split into vertical stripes, each joined with its own quadtree."""
    try:
        return _contains(points, polygons, max_size=max_size, batches=batches)
    except RuntimeError as e:
        if "quadtree" not in str(e):
            raise
        err = e
    import cudf
    import cupy as cp
    import cuspatial
    xy = cp.asarray(points.points.xy).reshape(-1, 2)
    order = cp.argsort(xy[:, 0])
    for k in (4, 16, 64, 256):
        try:
            parts = []
            for idx in cp.array_split(order, k):
                if len(idx) == 0:
                    continue
                sub = cuspatial.GeoSeries.from_points_xy(cp.ascontiguousarray(xy[idx]).ravel())
                r = _contains(sub, polygons, max_size=max_size, batches=batches).reset_index(drop=True)
                r["index_query"] = idx[r["index_query"].values]      # stripe position -> original point
                parts.append(r)
            say(f"point-in-polygon join on {len(xy)} points: quadtree failed, done in {k} stripes")
            return cudf.concat(parts, ignore_index=True)
        except RuntimeError as e:
            if "quadtree" not in str(e):
                raise
            err = e
    raise err


if MODE not in ("auto", "adaptive", "square"):
    raise SystemExit(f"SEGGER_TILING must be auto, adaptive or square, not {MODE!r}")
if MODE != "adaptive":
    ## both are looked up in their own module at run time, so replacing the module attribute is enough
    data_module.QuadTreeTiling = tiling_with_fallback          # ISTDataModule: tiling the graph
    query._points_in_polygons_contains = contains_in_stripes   # points_in_polygons + 'intersects' path
say(f"tiling mode: {MODE}" + (f", square tiles of {TILE_UM} um" if TILE_UM else ""))

if __name__ == "__main__":
    from segger.cli.main import app
    sys.argv[0] = "segger"
    app()
