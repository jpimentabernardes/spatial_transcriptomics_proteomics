"""Run the segger command line with a fallback for cuSpatial's "invalid quadtree" error.

On large slides (tens of millions of transcripts) cuSpatial can fail to build the quadtree that
segger uses to cut the slide into tiles:
    RuntimeError: cuSpatial invalid quadtree after sizes [...] (see segger issue #40)
segger also has a plain square-grid tiling, but its command line does not expose it. This file
replaces segger's QuadTreeTiling by a function that tries the quadtree first and, when it fails,
uses square tiles. Empty squares are dropped and, unless given, the side length is chosen so that
the densest squares hold about as many nodes as a quadtree tile would (--max-nodes-per-tile,
default 50k).

Started by 03b/03d inside segger's pixi environment, with the same arguments as `segger`:
    pixi run --frozen -e cuda121 python <this file> segment -i <xenium_dir> -o <out>
Settings (environment variables):
    SEGGER_TILING   auto (quadtree, square tiles if it fails; default) | adaptive (quadtree only,
                    = plain segger) | square (square tiles only)
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


if MODE not in ("auto", "adaptive", "square"):
    raise SystemExit(f"SEGGER_TILING must be auto, adaptive or square, not {MODE!r}")
## ISTDataModule looks QuadTreeTiling up in its own module when it tiles the graph
data_module.QuadTreeTiling = tiling_with_fallback
say(f"tiling mode: {MODE}" + (f", square tiles of {TILE_UM} um" if TILE_UM else ""))

if __name__ == "__main__":
    from segger.cli.main import app
    sys.argv[0] = "segger"
    app()
