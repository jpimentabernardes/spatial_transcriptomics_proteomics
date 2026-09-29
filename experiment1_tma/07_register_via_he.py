#!/usr/bin/env python
"""07_register_via_he.py -- align CellScape to Xenium THROUGH their H&E stains.

Both slides get an H&E after their run: the Xenium slide (post-Xenium H&E) and
the CellScape slide (CellScape is non-destructive). That gives two easy
same-section registrations and one consecutive-section registration between
images that look alike:

    CellScape cells --B--> H&E(CellScape) --D--> H&E(Xenium) --A--> Xenium um
      (CellScape DAPI)        same section    consecutive     same section   (Xenium DAPI)

  A  H&E(Xenium) -> Xenium DAPI       same section: haematoxylin ~ DAPI, near-exact
  B  H&E(CellScape) -> CellScape DAPI same section: near-exact
  D  H&E(CellScape) -> H&E(Xenium)    CONSECUTIVE sections, same stain on both sides.
                                      This is the hard step, and H&E vs H&E is the best
                                      possible contrast for it (much better than DAPI vs DAPI
                                      from two instruments, or RNA vs protein labels).

Each registration = global affine, then a rigid correction per TMA core (cores
shift/rotate independently between sections). A core correction is kept only if
it improves the image match (normalised cross-correlation) within the shift /
rotation limits; otherwise the core keeps the global transform.

Global initialisation:
  A, B  8 orientations (0/90/180/270 deg, with / without mirror) are tried at low
        resolution -- scanners often store H&E rotated or flipped.
  D     from matched TMA core centres (core ids from 02 and 06, mapped through A and B);
        falls back to the orientation search if < 3 cores match.
  Any of them can be forced with a landmark CSV (moving_x, moving_y, fixed_x, fixed_y; um).

All work is in microns: every image needs its pixel size (read from OME
metadata when present, else pass --*-px).

OUTPUTS (data/registration/)
  cellscape_aligned_<slide>.csv.gz  CellScape cells in Xenium um (x_aligned, y_aligned;
                                    x_global, y_global = globals only) -- same format as
                                    07_register_modalities.py, so 07b / 08 work unchanged
  transforms_via_he_<slide>.json    all transforms (A, B, D; global + per core)
  transforms_he_<slide>.json        H&E(Xenium) px -> Xenium um, for 09_he_regions.R
  registration_qc_via_he_<slide>.csv  per registration x core: NCC before/after, shift, accepted
  figures/07_via_he_<slide>_{A,B,D}.png  overlays (green = fixed, magenta = moving)

Example
  python 07_register_via_he.py --project-dir /work/Spatial_TMA --slide-id TMA_ORGAN \\
      --xenium-dir /work/.../output-XETG... \\
      --he-xenium /work/he/TMA_ORGAN_xenium_HE.ome.tif \\
      --cellscape-image /work/cellscape/TMA_ORGAN.ome.tif --cellscape-dapi DAPI \\
      --he-cellscape /work/he/TMA_ORGAN_cellscape_HE.ome.tif
"""

import argparse
import json
import os
import sys
from importlib import import_module

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "py"))
sys.path.insert(0, HERE)
from xenium_io import PyramidImage, find_morphology_channels, pixel_size_um  # noqa: E402

_reg = import_module("07_register_modalities")        # shared: apply_tf, robust_fit, rot_deg, lineage_agreement
apply_tf, robust_fit, rot_deg, lineage_agreement = _reg.apply_tf, _reg.robust_fit, _reg.rot_deg, _reg.lineage_agreement
icp_label_aware = _reg.icp_label_aware


# ---------------------------------------------------------------------------
# images as "nuclear signal" in microns
# ---------------------------------------------------------------------------
class NuclearView:
    """One image (DAPI or H&E) read as a nuclear-signal map at any resolution in um.
    H&E -> haematoxylin optical density (nuclei bright, like DAPI)."""

    def __init__(self, path, px0=None, kind="dapi", channel=0, name=""):
        self.path, self.kind, self.name = path, kind, name
        probe = PyramidImage(path, level=0, channel=channel)
        self.px0 = px0 or probe.pixel_size_um()
        if not self.px0:
            raise ValueError(f"{name}: no pixel size in OME metadata of {path}; pass it explicitly")
        ch = None if kind == "he" else channel
        self.levels = [PyramidImage(path, level=l, channel=ch) for l in range(probe.n_levels)]
        self.extent_um = (probe.width * self.px0, probe.height * self.px0)

    def read_um(self, x0, y0, x1, y1, res):
        """Nuclear signal on a grid of `res` um covering the box; returns (array, origin_um).
        Downsampling is area-averaging (a stride would alias single nuclei into noise)."""
        from skimage.transform import rescale
        best = self.levels[0]
        for lv in self.levels:                      # coarsest level still finer than res
            if self.px0 * lv.downsample <= res * 1.001:
                best = lv
        px = self.px0 * best.downsample
        step = max(1, int((res / px) // 4))          # only for huge single-level scans: keep 4x oversampling
        px_eff = px * step
        c0, r0 = int(np.floor(max(0, x0) / px)), int(np.floor(max(0, y0) / px))
        c1, r1 = int(np.ceil(x1 / px)), int(np.ceil(y1 / px))
        a = best.read(r0, r1, c0, c1, step=step).astype(np.float32)
        if a.size == 0 or min(a.shape[:2]) < 2:
            return np.zeros((1, 1), np.float32), (x0, y0)
        a = self._nuclear(a)
        f = px_eff / res
        if abs(f - 1) > 0.01:
            a = rescale(a, f, order=1, anti_aliasing=f < 1, preserve_range=True)
        return self._normalise(a).astype(np.float32), (c0 * px, r0 * px)

    def _nuclear(self, a):
        if self.kind == "he":
            if a.ndim == 3 and a.shape[-1] >= 3:
                from skimage.color import rgb2hed
                rgb = np.clip(a[..., :3] / (255.0 if a.max() > 1.5 else 1.0), 1 / 255, 1)
                h = rgb2hed(rgb)[..., 0]
                glass = (-np.log(rgb)).sum(-1) < 0.15            # unstained background
                bg = np.median(h[glass]) if glass.sum() > 100 else np.percentile(h, 5)
                return np.clip(h - bg, 0, None)
            return np.clip(a, 0, None)
        a = a if a.ndim == 2 else a[..., 0]
        return np.log1p(np.clip(a, 0, None))

    @staticmethod
    def _normalise(a):
        hi = np.percentile(a, 99.5) if a.size > 10 else 1
        return np.clip(a / max(hi, 1e-6), 0, 1)


def warp(moving: NuclearView, M, shape, origin, res):
    """Moving image resampled onto a fixed-frame grid (origin, res, shape) via M (moving um -> fixed um)."""
    from scipy.ndimage import map_coordinates
    h, w = shape
    yy, xx = np.mgrid[0:h, 0:w]
    pts = np.c_[origin[0] + (xx.ravel() + 0.5) * res, origin[1] + (yy.ravel() + 0.5) * res]
    mp = apply_tf(np.linalg.inv(M), pts)
    pad = 5 * res
    arr, mo = moving.read_um(mp[:, 0].min() - pad, mp[:, 1].min() - pad, mp[:, 0].max() + pad, mp[:, 1].max() + pad, res)
    out = map_coordinates(arr, [(mp[:, 1] - mo[1]) / res - 0.5, (mp[:, 0] - mo[0]) / res - 0.5],
                          order=1, mode="constant", cval=0)
    return out.reshape(h, w).astype(np.float32)


def ncc(a, b):
    m = (a > 0.05) | (b > 0.05)
    if m.sum() < 50:
        return np.nan
    a, b = a[m] - a[m].mean(), b[m] - b[m].mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum() + 1e-12))


# ---------------------------------------------------------------------------
# SimpleITK registration in physical um
# ---------------------------------------------------------------------------
def _sitk_image(arr, origin, res):
    import SimpleITK as sitk
    im = sitk.GetImageFromArray(arr.astype(np.float32))
    im.SetSpacing((res, res))
    im.SetOrigin((origin[0] + res / 2, origin[1] + res / 2))
    return im


def _tf_to_matrix(tf):
    import SimpleITK as sitk
    if tf.GetName() == "CompositeTransform":
        ct = sitk.CompositeTransform(tf)
        M = np.eye(3)
        for i in range(ct.GetNumberOfTransforms()):
            M = M @ _tf_to_matrix(ct.GetNthTransform(i))
        return M
    kinds = {"AffineTransform": sitk.AffineTransform, "Euler2DTransform": sitk.Euler2DTransform,
             "Similarity2DTransform": sitk.Similarity2DTransform}
    t = kinds[tf.GetName()](tf)
    A = np.array(t.GetMatrix()).reshape(2, 2)
    c, tr = np.array(t.GetCenter()), np.array(t.GetTranslation())
    M = np.eye(3)
    M[:2, :2] = A
    M[:2, 2] = c + tr - A @ c
    return M


def register_pair(F, fo, Mv, mo, res, init_m2f=np.eye(3), kind="affine", iters=200, shrink=(4, 2, 1)):
    """Returns M (moving um -> fixed um). SimpleITK optimises fixed -> moving; we invert."""
    import SimpleITK as sitk
    f, m = _sitk_image(F, fo, res), _sitk_image(Mv, mo, res)
    T0 = np.linalg.inv(init_m2f)
    if kind == "affine":
        tx = sitk.AffineTransform(2)
        tx.SetMatrix(T0[:2, :2].ravel().tolist())
        tx.SetTranslation(T0[:2, 2].tolist())
    else:                                            # rigid, around the image centre
        tx = sitk.Euler2DTransform()
        c = np.array(fo) + np.array(F.shape[::-1]) * res / 2
        tx.SetCenter(c.tolist())
        tx.SetAngle(float(np.arctan2(T0[1, 0], T0[0, 0])))
        tx.SetTranslation((T0[:2, :2] @ c + T0[:2, 2] - c).tolist())
    reg = sitk.ImageRegistrationMethod()
    reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins=32)
    reg.SetMetricSamplingStrategy(reg.RANDOM)
    reg.SetMetricSamplingPercentage(0.3, seed=7)
    reg.SetInterpolator(sitk.sitkLinear)
    reg.SetOptimizerAsRegularStepGradientDescent(learningRate=1.0, minStep=1e-4, numberOfIterations=iters,
                                                 relaxationFactor=0.6)
    reg.SetOptimizerScalesFromPhysicalShift()
    small = min(F.shape + Mv.shape)
    sh = [s for s in shrink if small / s >= 24] or [1]
    reg.SetShrinkFactorsPerLevel(sh)
    reg.SetSmoothingSigmasPerLevel([max(0, s - 1) for s in sh])
    reg.SetInitialTransform(tx, inPlace=False)
    try:
        out = reg.Execute(f, m)
    except RuntimeError as e:
        print(f"    registration failed: {e}")
        return init_m2f
    return np.linalg.inv(_tf_to_matrix(out))


def orientation_candidates(F, fo, Mv, mo, res):
    """Rigid initial guesses: 4 rotations x mirror, aligning signal-weighted centroids."""
    def centroid(a, o):
        w = np.clip(a - np.median(a), 0, None)
        yy, xx = np.mgrid[0:a.shape[0], 0:a.shape[1]]
        s = w.sum() + 1e-9
        return np.array([o[0] + ((xx * w).sum() / s + 0.5) * res, o[1] + ((yy * w).sum() / s + 0.5) * res])
    cf, cm = centroid(F, fo), centroid(Mv, mo)
    T = lambda v: np.array([[1, 0, v[0]], [0, 1, v[1]], [0, 0, 1.0]])  # noqa: E731
    out = []
    for mirror in (False, True):
        for k in range(4):
            t = np.radians(90 * k)
            R = np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1.0]])
            Mi = np.diag([-1.0 if mirror else 1.0, 1, 1])
            out.append((f"rot{90 * k}{'+mirror' if mirror else ''}", T(cf) @ R @ Mi @ T(-cm)))
    return out


def register_global(fixed: NuclearView, moving: NuclearView, res, init=None, label=""):
    F, fo = fixed.read_um(0, 0, *fixed.extent_um, res)
    Mv, mo = moving.read_um(0, 0, *moving.extent_um, res)
    cands = [("init", init)] if init is not None else orientation_candidates(F, fo, Mv, mo, res)
    best = None
    for name, c in cands:
        M = register_pair(F, fo, Mv, mo, res, c, kind="affine", iters=150)
        score = ncc(F, warp(moving, M, F.shape, fo, res))
        if best is None or (not np.isnan(score) and score > best[2]):
            best = (name, M, score)
    name, M, score = best
    # second pass at double resolution from the winner
    F2, fo2 = fixed.read_um(0, 0, *fixed.extent_um, res / 2)
    Mv2, mo2 = moving.read_um(0, 0, *moving.extent_um, res / 2)
    M2 = register_pair(F2, fo2, Mv2, mo2, res / 2, M, kind="affine", iters=150)
    s2 = ncc(F2, warp(moving, M2, F2.shape, fo2, res / 2))
    s1 = ncc(F2, warp(moving, M, F2.shape, fo2, res / 2))
    if not np.isnan(s2) and s2 >= s1:
        M, score = M2, s2
    det = np.linalg.det(M[:2, :2])
    print(f"  {label}: start={name}, NCC={score:.3f}, rotation={rot_deg(M):.1f} deg, "
          f"scale={np.sqrt(abs(det)):.3f}{' MIRRORED' if det < 0 else ''}")
    return M, score


def refine_cores(fixed: NuclearView, moving: NuclearView, G, boxes: dict, res, max_shift, max_rot, label="",
                 smooth_um: float = 0.0):
    """Rigid per-core correction C (fixed um -> fixed um) on top of G; kept only if NCC improves.

    smooth_um > 0 blurs both images first. Use it for CONSECUTIVE sections: their nuclei are
    different cells, so only tissue architecture (glands, crypts, vessels, lymphoid follicles)
    can be matched -- at nucleus scale the two images do not correspond."""
    from scipy.ndimage import gaussian_filter
    sm = (lambda a: gaussian_filter(a, smooth_um / res)) if smooth_um > 0 else (lambda a: a)  # noqa: E731
    corr, rows = {}, []
    for cid, (x0, y0, x1, y1) in boxes.items():
        F, fo = fixed.read_um(x0, y0, x1, y1, res)
        if F.size < 400 or F.max() == 0:
            continue
        F = sm(F)
        W = sm(warp(moving, G, F.shape, fo, res))
        n0 = ncc(F, W)
        C = register_pair(F, fo, W, fo, res, np.eye(3), kind="rigid", iters=200, shrink=(2, 1))
        centre = np.array([[(x0 + x1) / 2, (y0 + y1) / 2]])
        shift = float(np.linalg.norm(apply_tf(C, centre) - centre))
        rot = rot_deg(C)
        n1 = ncc(F, sm(warp(moving, C @ G, F.shape, fo, res)))
        ok = (shift <= max_shift and abs(rot) <= max_rot and not np.isnan(n1) and n1 > 0.2
              and (np.isnan(n0) or n1 > n0 + 0.005))
        if ok:
            corr[cid] = C
        rows.append(dict(registration=label, core_id=cid, ncc_global=n0, ncc_refined=n1,
                         shift_um=shift, rotation_deg=rot, accepted=ok))
    acc = sum(r["accepted"] for r in rows)
    print(f"  {label}: {acc}/{len(rows)} cores refined")
    return corr, rows


def overlay_png(path, fixed, moving, G, corr, boxes, res, title, n=12):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ids = list(boxes)[:n]
    if not ids:
        return
    nc = min(6, len(ids))
    nr = int(np.ceil(len(ids) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(2.6 * nc, 2.6 * nr + 0.4), squeeze=False)
    for ax, cid in zip(axes.ravel(), ids):
        F, fo = fixed.read_um(*boxes[cid], res)
        W = warp(moving, corr.get(cid, np.eye(3)) @ G, F.shape, fo, res)
        ax.imshow(np.dstack([W, F, W]))          # green = fixed, magenta = moving, white = overlap
        ax.set_title(f"{cid}{' (refined)' if cid in corr else ''}", fontsize=7)
        ax.axis("off")
    for ax in axes.ravel()[len(ids):]:
        ax.axis("off")
    fig.suptitle(title + "  --  green: fixed, magenta: moving, white: overlap", fontsize=9)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
def boxes_from_cores(df, factor=1.25):
    out = {}
    for r in df.itertuples(index=False):
        h = r.radius * factor
        out[str(r.core_id)] = (r.x_center - h, r.y_center - h, r.x_center + h, r.y_center + h)
    return out


def landmark_tf(path):
    lm = pd.read_csv(path)
    M, res, keep, model = robust_fit(lm[["moving_x", "moving_y"]].to_numpy(), lm[["fixed_x", "fixed_y"]].to_numpy())
    print(f"  landmarks {os.path.basename(path)}: {model}, median residual {np.median(res[keep]):.1f} um")
    return M


def full(G, corr, cid):
    return corr[cid] @ G if cid in corr else G


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-dir", required=True)
    p.add_argument("--slide-id", required=True, help="Xenium slide_id (slides.csv)")
    p.add_argument("--xenium-dir", required=True)
    p.add_argument("--he-xenium", required=True, help="H&E of the Xenium slide (OME-TIFF, pyramidal preferred)")
    p.add_argument("--he-xenium-px", type=float, help="um/px of that H&E (default: OME metadata)")
    p.add_argument("--cellscape-image", required=True, help="CellScape OME-TIFF with a DAPI channel")
    p.add_argument("--cellscape-dapi", default="DAPI", help="channel name or index")
    p.add_argument("--cellscape-px", type=float, help="um/px of the CellScape image (default: OME metadata)")
    p.add_argument("--he-cellscape", required=True, help="H&E of the CellScape slide")
    p.add_argument("--he-cellscape-px", type=float)
    p.add_argument("--landmarks-a", help="CSV to force A: H&E(Xenium) um -> Xenium um")
    p.add_argument("--landmarks-b", help="CSV to force B: H&E(CellScape) um -> CellScape um")
    p.add_argument("--landmarks-d", help="CSV to force D: H&E(CellScape) um -> H&E(Xenium) um")
    p.add_argument("--overview-res", type=float, default=20.0, help="um/px for global registration")
    p.add_argument("--core-res", type=float, default=2.0, help="um/px, per-core refinement of A and B (same section)")
    p.add_argument("--core-res-consecutive", type=float, default=6.0, help="um/px, per-core refinement of D")
    p.add_argument("--smooth-consecutive-um", type=float, default=12.0,
                   help="blur for D: match tissue architecture, not individual nuclei")
    p.add_argument("--no-core-refine", action="store_true")
    p.add_argument("--max-shift-same", type=float, default=40.0, help="per-core limit, same-section (A, B)")
    p.add_argument("--max-shift-consecutive", type=float, default=150.0, help="per-core limit, H&E->H&E (D)")
    p.add_argument("--max-rot-deg", type=float, default=15.0)
    p.add_argument("--polish-with-labels", action="store_true",
                   help="after the H&E chain, a small label-aware rigid ICP per core (as in 07_register_modalities "
                        "--refine points); kept only where it raises cross-modal lineage agreement")
    p.add_argument("--polish-max-shift", type=float, default=15.0)
    return p.parse_args()


def main():
    a = parse_args()
    reg_dir = os.path.join(a.project_dir, "data", "registration")
    tab_dir = os.path.join(a.project_dir, "tables")
    fig_dir = os.path.join(a.project_dir, "figures")
    for d in (reg_dir, fig_dir):
        os.makedirs(d, exist_ok=True)

    # ---- images
    xd = NuclearView(find_morphology_channels(a.xenium_dir)["dapi"], pixel_size_um(a.xenium_dir), "dapi", 0, "Xenium DAPI")
    hx = NuclearView(a.he_xenium, a.he_xenium_px, "he", name="H&E (Xenium slide)")
    probe = PyramidImage(a.cellscape_image)
    names = probe.channel_names()
    ch = names.index(a.cellscape_dapi) if a.cellscape_dapi in names else int(a.cellscape_dapi) if a.cellscape_dapi.isdigit() else 0
    cd = NuclearView(a.cellscape_image, a.cellscape_px, "dapi", ch, "CellScape DAPI")
    hc = NuclearView(a.he_cellscape, a.he_cellscape_px, "he", name="H&E (CellScape slide)")
    for v in (xd, hx, cd, hc):
        print(f"{v.name}: {v.px0:.4f} um/px, {v.extent_um[0]:.0f} x {v.extent_um[1]:.0f} um")

    # ---- cores (Xenium frame from 02, CellScape frame from 06)
    xc = pd.read_csv(os.path.join(tab_dir, "02_cores_all.csv"))
    xc = xc[(xc.slide_id == a.slide_id) & xc.core_id.notna()]
    cc = pd.read_csv(os.path.join(tab_dir, "06_cs_cores_all.csv"))
    cc = cc[(cc.slide_id == a.slide_id) & cc.core_id.notna()]
    x_boxes, c_boxes = boxes_from_cores(xc), boxes_from_cores(cc)
    qc_rows = []
    refine = not a.no_core_refine

    # ---- A: H&E(Xenium) -> Xenium DAPI (same section)
    print("A  H&E(Xenium) -> Xenium")
    GA = landmark_tf(a.landmarks_a) if a.landmarks_a else register_global(xd, hx, a.overview_res, label="A")[0]
    CA, rows = refine_cores(xd, hx, GA, x_boxes, a.core_res, a.max_shift_same, a.max_rot_deg, "A") if refine else ({}, [])
    qc_rows += rows

    # ---- B: H&E(CellScape) -> CellScape DAPI (same section)
    print("B  H&E(CellScape) -> CellScape")
    GB = landmark_tf(a.landmarks_b) if a.landmarks_b else register_global(cd, hc, a.overview_res, label="B")[0]
    CB, rows = refine_cores(cd, hc, GB, c_boxes, a.core_res, a.max_shift_same, a.max_rot_deg, "B") if refine else ({}, [])
    qc_rows += rows

    # ---- D: H&E(CellScape) -> H&E(Xenium) (consecutive sections)
    print("D  H&E(CellScape) -> H&E(Xenium)")
    shared = sorted(set(xc.core_id) & set(cc.core_id))
    hx_centres = {cid: apply_tf(np.linalg.inv(full(GA, CA, cid)), xc.set_index("core_id").loc[[cid], ["x_center", "y_center"]].to_numpy())[0] for cid in shared}
    hc_centres = {cid: apply_tf(np.linalg.inv(full(GB, CB, cid)), cc.set_index("core_id").loc[[cid], ["x_center", "y_center"]].to_numpy())[0] for cid in shared}
    if a.landmarks_d:
        GD = landmark_tf(a.landmarks_d)
    else:
        init = None
        if len(shared) >= 3:
            init, res_, keep, model = robust_fit(np.array([hc_centres[c] for c in shared]),
                                                 np.array([hx_centres[c] for c in shared]))
            print(f"  D init from {len(shared)} matched cores ({model}), median residual {np.median(res_[keep]):.1f} um")
        GD = register_global(hx, hc, a.overview_res, init=init, label="D")[0]
    radius = xc.set_index("core_id").radius
    d_boxes = {cid: (c[0] - 1.25 * radius[cid], c[1] - 1.25 * radius[cid], c[0] + 1.25 * radius[cid], c[1] + 1.25 * radius[cid])
               for cid, c in hx_centres.items()}
    CD, rows = refine_cores(hx, hc, GD, d_boxes, a.core_res_consecutive, a.max_shift_consecutive, a.max_rot_deg, "D",
                            smooth_um=a.smooth_consecutive_um) if refine else ({}, [])
    qc_rows += rows

    # ---- compose: CellScape um -> H&E(C) um -> H&E(X) um -> Xenium um, per core
    mov = pd.read_csv(os.path.join(reg_dir, f"cellscape_cells_{a.slide_id}.csv.gz"))
    xy = mov[["x_um", "y_um"]].to_numpy(float)
    glob = GA @ GD @ np.linalg.inv(GB)
    g = apply_tf(glob, xy)
    mov["x_global"], mov["y_global"] = g[:, 0], g[:, 1]
    out = g.copy()
    per_core = {}
    for cid in mov.core_id.dropna().unique():
        T = full(GA, CA, cid) @ full(GD, CD, cid) @ np.linalg.inv(full(GB, CB, cid))
        m = (mov.core_id == cid).to_numpy()
        out[m] = apply_tf(T, xy[m])
        per_core[str(cid)] = T.tolist()
    mov["x_aligned"], mov["y_aligned"] = out[:, 0], out[:, 1]
    mov["refine_method"] = "via_he"

    fix = pd.read_csv(os.path.join(reg_dir, f"xenium_cells_{a.slide_id}.csv.gz"))
    if a.polish_with_labels:
        n_pol = 0
        for cid in sorted(set(mov.core_id.dropna()) & set(fix.core_id.dropna())):
            mm, ff = (mov.core_id == cid).to_numpy(), (fix.core_id == cid).to_numpy()
            mxy = mov.loc[mm, ["x_aligned", "y_aligned"]].to_numpy()
            fxy = fix.loc[ff, ["x_um", "y_um"]].to_numpy()
            ml, fl = mov.loc[mm, "lineage"].astype(str).to_numpy(), fix.loc[ff, "lineage"].astype(str).to_numpy()
            P = icp_label_aware(mxy, ml, fxy, fl, max_dist=2 * a.polish_max_shift)
            c = fxy.mean(0, keepdims=True)
            if np.linalg.norm(apply_tf(P, c) - c) > a.polish_max_shift:
                continue
            new = apply_tf(P, mxy)
            if lineage_agreement(new, ml, fxy, fl)[1] > lineage_agreement(mxy, ml, fxy, fl)[1]:
                mov.loc[mm, ["x_aligned", "y_aligned"]] = new
                mov.loc[mm, "refine_method"] = "via_he+labels"
                per_core[str(cid)] = (P @ np.asarray(per_core[str(cid)])).tolist()
                n_pol += 1
        print(f"label polish accepted in {n_pol} cores")
    # ---- QC: cross-modal lineage agreement per core (independent of the images)
    for cid in sorted(set(mov.core_id.dropna()) & set(fix.core_id.dropna())):
        mm, ff = (mov.core_id == cid).to_numpy(), (fix.core_id == cid).to_numpy()
        f_xy, f_lab = fix.loc[ff, ["x_um", "y_um"]].to_numpy(), fix.loc[ff, "lineage"].astype(str).to_numpy()
        m_lab = mov.loc[mm, "lineage"].astype(str).to_numpy()
        _, lift_g, nn_g = lineage_agreement(mov.loc[mm, ["x_global", "y_global"]].to_numpy(), m_lab, f_xy, f_lab)
        _, lift_r, nn_r = lineage_agreement(mov.loc[mm, ["x_aligned", "y_aligned"]].to_numpy(), m_lab, f_xy, f_lab)
        qc_rows.append(dict(registration="final", core_id=cid, lift_global=lift_g, lift_refined=lift_r,
                            median_nn_global=nn_g, median_nn_refined=nn_r))
    mov.to_csv(os.path.join(reg_dir, f"cellscape_aligned_{a.slide_id}.csv.gz"), index=False)
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(os.path.join(reg_dir, f"registration_qc_via_he_{a.slide_id}.csv"), index=False)
    fin = qc[qc.registration == "final"]
    if len(fin):
        print(f"final: median lineage-agreement lift {fin.lift_global.median():.2f} (globals only) -> "
              f"{fin.lift_refined.median():.2f} (per-core); median NN distance {fin.median_nn_refined.median():.1f} um")

    # ---- transforms (json) + transforms_he for 09
    tj = lambda d: {k: np.asarray(v).tolist() for k, v in d.items()}  # noqa: E731
    with open(os.path.join(reg_dir, f"transforms_via_he_{a.slide_id}.json"), "w") as fh:
        json.dump({"units": "um", "chain": "cellscape -> inv(B) -> D -> A -> xenium",
                   "A_he_xenium_to_xenium": {"global": GA.tolist(), "cores": tj(CA)},
                   "B_he_cellscape_to_cellscape": {"global": GB.tolist(), "cores": tj(CB)},
                   "D_he_cellscape_to_he_xenium": {"global": GD.tolist(), "cores": tj(CD)},
                   "cellscape_to_xenium": {"global": glob.tolist(), "cores": per_core}}, fh, indent=1)
    S = np.diag([hx.px0, hx.px0, 1.0])
    with open(os.path.join(reg_dir, f"transforms_he_{a.slide_id}.json"), "w") as fh:
        json.dump({"moving": "he_px_level0", "fixed": "xenium_um", "global_model": "via_he_A",
                   "global": (GA @ S).tolist(), "cores": tj(CA)}, fh, indent=1)

    for lab, fixed, moving, G, C, boxes in (("A", xd, hx, GA, CA, x_boxes), ("B", cd, hc, GB, CB, c_boxes),
                                            ("D", hx, hc, GD, CD, d_boxes)):
        overlay_png(os.path.join(fig_dir, f"07_via_he_{a.slide_id}_{lab}.png"), fixed, moving, G, C, boxes,
                    a.core_res * 2, f"{a.slide_id} {lab}: {moving.name} -> {fixed.name}")
    print("Done. Next: 07b_nonrigid_refine.py (optional) or 08_integrate_rna_protein.R")


if __name__ == "__main__":
    main()
