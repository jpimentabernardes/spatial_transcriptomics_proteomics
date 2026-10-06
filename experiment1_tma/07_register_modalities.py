#!/usr/bin/env python
"""07_register_modalities.py -- bring CellScape (consecutive section) and H&E
into the Xenium coordinate frame (microns, Xenium = fixed reference).

Used now for --mode he (H&E of the Xenium slide -> Xenium, for 09). The CellScape
steps (06, 07 via H&E, 07b, 08) are removed for now; --mode cellscape needs the 06 export.

Strategy for TMAs (coarse -> fine, every step checked):

 1. GLOBAL affine from matched TMA cores. Both modalities were dearrayed with
    the same TMA map (02 / 06), so every core id present in both gives one
    landmark pair (core centres) for free. Robust least squares with outlier
    rejection. If a section is mirrored the affine has a negative
    determinant -- reported, and handled.
    Fallback: a landmark CSV (moving_x, moving_y, fixed_x, fixed_y).
 2. PER-CORE refinement. Cores in serial sections shift/rotate a little
    independently (they are separate plugs of tissue), so one global affine
    is not enough at cell scale. Options (--refine):
      points  label-aware rigid ICP on cell centroids: each CellScape cell is
              pulled towards the nearest Xenium cell OF THE SAME LINEAGE.
              Needs no images. Default.
      image   intensity registration (SimpleITK, mutual information) of the
              CellScape DAPI (or H&E haematoxylin) to the Xenium DAPI crop.
              Most precise when both images are available.
      none    global only.
    A refinement is accepted only if it improves cross-modal lineage
    agreement and stays within --max-shift-um / --max-rot-deg. Otherwise the
    core keeps the global transform (and is listed in the QC table).

Serial sections are ~5 um apart: even a perfect registration aligns TISSUE
STRUCTURE, not the same cells. 08 treats the match accordingly
(neighbourhood-level integration, distance cut-offs).

Modes
-----
  --mode cellscape  moving = CellScape cells (06 export). Output:
                    data/registration/cellscape_aligned_<slide>.csv.gz,
                    transforms_cellscape_<slide>.json, registration_qc_cellscape_<slide>.csv
  --mode he         moving = H&E image (pixels). If H&E was done on the Xenium
                    slide itself after the run (standard post-Xenium H&E) this
                    is SAME-section: the global affine from landmarks or from a
                    Xenium Explorer image-alignment matrix is usually enough.
                    Output: transforms_he_<slide>.json (H&E px -> Xenium um),
                    used by 09_he_regions.R.

Examples
--------
  python 07_register_modalities.py --project-dir /work/Spatial_TMA --slide-id TMA_ORGAN --mode cellscape --refine points
  python 07_register_modalities.py --project-dir /work/Spatial_TMA --slide-id TMA_ORGAN --mode cellscape --refine image \
      --xenium-dir /work/.../output-XETG... --moving-image /work/.../image.ome.tif --moving-pixel-size-um 0.325 --moving-channel DAPI
  python 07_register_modalities.py --project-dir /work/Spatial_TMA --slide-id TMA_ORGAN --mode he \
      --xenium-dir /work/.../output-XETG... --moving-image he.ome.tif --xe-alignment-csv he_imagealignment.csv --refine image
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "py"))
from xenium_io import PyramidImage, find_morphology_channels, pixel_size_um  # noqa: E402


# ---------------------------------------------------------------------------
# transforms (3x3 homogeneous, applied to column vectors [x, y, 1])
# ---------------------------------------------------------------------------
def apply_tf(M, xy):
    xy = np.asarray(xy, dtype=float)
    h = np.c_[xy, np.ones(len(xy))] @ M.T
    return h[:, :2]


def fit_affine(src, dst):
    A = np.c_[src, np.ones(len(src))]
    P, *_ = np.linalg.lstsq(A, dst, rcond=None)       # 3x2
    M = np.eye(3)
    M[:2, :] = P.T
    return M


def fit_similarity(src, dst, allow_reflection=True):
    """Umeyama: rotation + uniform scale + translation (+ reflection if it fits better)."""
    mu_s, mu_d = src.mean(0), dst.mean(0)
    s, d = src - mu_s, dst - mu_d
    U, S, Vt = np.linalg.svd(d.T @ s / len(src))
    D = np.eye(2)
    if np.linalg.det(U @ Vt) < 0 and not allow_reflection:
        D[1, 1] = -1
    R = U @ D @ Vt
    scale = np.trace(np.diag(S) @ D) / s.var(0).sum()
    M = np.eye(3)
    M[:2, :2] = scale * R
    M[:2, 2] = mu_d - scale * R @ mu_s
    return M


def fit_rigid(src, dst):
    mu_s, mu_d = src.mean(0), dst.mean(0)
    U, _, Vt = np.linalg.svd((dst - mu_d).T @ (src - mu_s))
    D = np.diag([1, np.sign(np.linalg.det(U @ Vt))])
    R = U @ D @ Vt
    M = np.eye(3)
    M[:2, :2] = R
    M[:2, 2] = mu_d - R @ mu_s
    return M


def robust_fit(src, dst, model="auto", n_iter=5):
    """Least squares with iterative rejection of pairs > 3 MAD (min 50 um)."""
    src, dst = np.asarray(src, float), np.asarray(dst, float)
    if len(src) < 3:
        raise ValueError(f"need >= 3 landmark pairs, got {len(src)}")
    if model == "auto":
        model = "affine" if len(src) >= 6 else "similarity"
    fit = fit_affine if model == "affine" else fit_similarity
    keep = np.ones(len(src), bool)
    for _ in range(n_iter):
        M = fit(src[keep], dst[keep])
        res = np.linalg.norm(apply_tf(M, src) - dst, axis=1)
        thr = max(50.0, np.median(res[keep]) + 3 * 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep]))))
        new_keep = res <= thr
        if new_keep.sum() < 3 or (new_keep == keep).all():
            break
        keep = new_keep
    M = fit(src[keep], dst[keep])
    res = np.linalg.norm(apply_tf(M, src) - dst, axis=1)
    return M, res, keep, model


def rot_deg(M):
    return float(np.degrees(np.arctan2(M[1, 0], M[0, 0])))


# ---------------------------------------------------------------------------
# QC metric: cross-modal lineage agreement
# ---------------------------------------------------------------------------
def lineage_agreement(mov_xy, mov_lab, fix_xy, fix_lab, radius=20.0):
    """Fraction of moving cells whose nearest fixed cell (within radius) has the
    same lineage, divided by what random pairing would give. 1 = no better
    than chance; well-aligned tissue with consistent annotation >> 1."""
    if len(mov_xy) == 0 or len(fix_xy) == 0:
        return np.nan, np.nan, np.nan
    d, idx = cKDTree(fix_xy).query(mov_xy, k=1, distance_upper_bound=radius)
    ok = np.isfinite(d)
    if ok.sum() < 20:
        return np.nan, np.nan, float(np.median(d[ok])) if ok.any() else np.nan
    agree = np.mean(mov_lab[ok] == fix_lab[idx[ok]])
    labs = np.union1d(mov_lab, fix_lab)
    pm = np.array([np.mean(mov_lab[ok] == l) for l in labs])
    pf = np.array([np.mean(fix_lab == l) for l in labs])
    expected = float((pm * pf).sum())
    return float(agree), float(agree / expected) if expected > 0 else np.nan, float(np.median(d[ok]))


# ---------------------------------------------------------------------------
# per-core refinement: label-aware rigid ICP
# ---------------------------------------------------------------------------
def icp_label_aware(mov_xy, mov_lab, fix_xy, fix_lab, max_dist=40.0, n_iter=40, trim=0.7):
    trees = {l: (cKDTree(fix_xy[fix_lab == l]), fix_xy[fix_lab == l])
             for l in np.unique(fix_lab) if (fix_lab == l).sum() >= 5}
    T = np.eye(3)
    cur = mov_xy.copy()
    for it in range(n_iter):
        src, dst, dist = [], [], []
        for l, (tree, pts) in trees.items():
            m = mov_lab == l
            if m.sum() < 5:
                continue
            d, i = tree.query(cur[m], k=1, distance_upper_bound=max_dist)
            ok = np.isfinite(d)
            src.append(cur[m][ok]); dst.append(pts[i[ok]]); dist.append(d[ok])
        if not src:
            break
        src, dst, dist = np.vstack(src), np.vstack(dst), np.concatenate(dist)
        if len(src) < 30:
            break
        keep = dist <= np.quantile(dist, trim)          # trimmed ICP: serial sections have no 1:1 partners
        step = fit_rigid(src[keep], dst[keep])
        cur = apply_tf(step, cur)
        T = step @ T
        if np.hypot(*step[:2, 2]) < 0.05 and abs(rot_deg(step)) < 0.01:
            break
    return T


# ---------------------------------------------------------------------------
# per-core refinement: image-based (SimpleITK, mutual information)
# ---------------------------------------------------------------------------
def choose_level(path, base_px, target_um, channel=0):
    best, best_img = None, None
    lvl = 0
    while True:
        img = PyramidImage(path, level=lvl, channel=channel)
        px = base_px * img.downsample
        if best is None or abs(px - target_um) < abs(best - target_um):
            best, best_img = px, img
        if lvl >= img.n_levels - 1:
            break
        lvl += 1
    return best_img, best


def moving_intensity(img, y0, y1, x0, x1, he=False):
    a = img.read(y0, y1, x0, x1).astype(np.float32)
    if he:
        from skimage.color import rgb2hed
        if a.ndim == 3 and a.shape[-1] >= 3:
            a = rgb2hed(np.clip(a[..., :3] / 255.0, 0, 1))[..., 0]   # haematoxylin ~ nuclei ~ DAPI
    return a


def refine_core_image(box, G, fixed_img, fixed_px, mov_img, mov_unit_to_px0, he=False, max_shift_um=100.0):
    """Returns 3x3 correction C (Xenium um -> Xenium um) so that C @ G maps
    moving -> fixed, or None if registration failed."""
    import SimpleITK as sitk
    from scipy.ndimage import map_coordinates

    xmin, xmax, ymin, ymax = box
    fy0, fy1 = int(ymin / fixed_px), int(np.ceil(ymax / fixed_px))
    fx0, fx1 = int(xmin / fixed_px), int(np.ceil(xmax / fixed_px))
    fixed = fixed_img.read(fy0, fy1, fx0, fx1).astype(np.float32)
    if fixed.size == 0 or fixed.max() == 0:
        return None
    # physical coords of fixed pixel centres (Xenium um)
    yy, xx = np.mgrid[fy0:fy0 + fixed.shape[0], fx0:fx0 + fixed.shape[1]]
    pts = np.c_[(xx.ravel() + 0.5) * fixed_px, (yy.ravel() + 0.5) * fixed_px]
    mov_units = apply_tf(np.linalg.inv(G), pts)
    mov_pix = mov_units * mov_unit_to_px0 / mov_img.downsample   # moving pyramid-level pixels
    my0, mx0 = np.floor(mov_pix[:, 1].min()) - 2, np.floor(mov_pix[:, 0].min()) - 2
    my1, mx1 = np.ceil(mov_pix[:, 1].max()) + 2, np.ceil(mov_pix[:, 0].max()) + 2
    mov = moving_intensity(mov_img, my0, my1, mx0, mx1, he=he)
    if mov.size == 0:
        return None
    warped = map_coordinates(mov, [mov_pix[:, 1] - max(0, my0), mov_pix[:, 0] - max(0, mx0)],
                             order=1, mode="constant", cval=0).reshape(fixed.shape).astype(np.float32)
    if warped.max() == 0:
        return None

    f = sitk.GetImageFromArray(fixed)
    m = sitk.GetImageFromArray(warped)
    for im in (f, m):
        im.SetSpacing((fixed_px, fixed_px))
        im.SetOrigin(((fx0 + 0.5) * fixed_px, (fy0 + 0.5) * fixed_px))
    init = sitk.CenteredTransformInitializer(f, m, sitk.Euler2DTransform(),
                                             sitk.CenteredTransformInitializerFilter.GEOMETRY)
    reg = sitk.ImageRegistrationMethod()
    reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins=32)
    reg.SetMetricSamplingStrategy(reg.RANDOM)
    reg.SetMetricSamplingPercentage(0.2, seed=42)
    reg.SetInterpolator(sitk.sitkLinear)
    reg.SetOptimizerAsRegularStepGradientDescent(learningRate=2.0, minStep=1e-3, numberOfIterations=300)
    reg.SetOptimizerScalesFromPhysicalShift()
    reg.SetShrinkFactorsPerLevel([4, 2, 1])
    reg.SetSmoothingSigmasPerLevel([2, 1, 0])
    reg.SetInitialTransform(init, inPlace=False)
    try:
        tf = reg.Execute(f, m)
    except RuntimeError:
        return None
    e = sitk.Euler2DTransform(sitk.CompositeTransform(tf).GetNthTransform(0)) \
        if tf.GetName() == "CompositeTransform" else sitk.Euler2DTransform(tf)
    R = np.array(e.GetMatrix()).reshape(2, 2)
    c = np.array(e.GetCenter())
    t = np.array(e.GetTranslation())
    MT = np.eye(3)
    MT[:2, :2] = R
    MT[:2, 2] = c + t - R @ c          # T(x) = R(x - c) + c + t : fixed -> warped-moving
    C = np.linalg.inv(MT)              # warped-moving (Xenium frame) -> fixed
    centre = np.array([[(xmin + xmax) / 2, (ymin + ymax) / 2]])
    if np.linalg.norm(apply_tf(C, centre) - centre) > max_shift_um:
        return None
    return C


# ---------------------------------------------------------------------------
# QC plots
# ---------------------------------------------------------------------------
def plot_core_grid(path, cores, fix, mov, xcol, ycol, title, max_cores=30):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cores = list(cores)[:max_cores]
    n = len(cores)
    if n == 0:
        return
    nc = min(6, n)
    nr = int(np.ceil(n / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(3.2 * nc, 3.2 * nr), squeeze=False)
    labs = sorted(set(mov["lineage"].dropna()) | set(fix["lineage"].dropna()))
    cmap = plt.get_cmap("tab10")
    colors = {l: cmap(i % 10) for i, l in enumerate(labs)}
    for ax, cid in zip(axes.ravel(), cores):
        f = fix[fix.core_id == cid]
        m = mov[mov.core_id == cid]
        size = max(0.3, min(4.0, 3000 / max(len(f), 1)))     # readable for small and huge cores
        ax.scatter(f.x_um, f.y_um, s=size * 1.6, c="#b0b0b0", linewidths=0)
        ax.scatter(m[xcol], m[ycol], s=size, c=[colors.get(l, (0, 0, 0, 1)) for l in m.lineage], linewidths=0)
        ax.set_title(cid, fontsize=7)
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    handles = [plt.Line2D([], [], marker="o", ls="", color=colors[l], label=l) for l in labs]
    fig.legend(handles=handles, loc="lower center", ncol=min(6, len(labs)), fontsize=7)
    fig.suptitle(title, fontsize=9)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-dir", required=True)
    p.add_argument("--slide-id", required=True, help="Xenium slide_id (slides.csv)")
    p.add_argument("--mode", choices=["cellscape", "he"], default="cellscape")
    p.add_argument("--refine", choices=["none", "points", "image"], default="points")
    p.add_argument("--global-model", choices=["auto", "similarity", "affine"], default="auto")
    p.add_argument("--landmarks", help="CSV moving_x,moving_y,fixed_x,fixed_y (moving: CellScape um or H&E px)")
    p.add_argument("--xe-alignment-csv", help="H&E mode: Xenium Explorer image-alignment 3x3 matrix")
    p.add_argument("--xenium-dir", help="needed for --refine image and --mode he")
    p.add_argument("--moving-image", help="CellScape OME-TIFF (DAPI) or H&E image")
    p.add_argument("--moving-channel", default="0", help="channel name or index of DAPI in the moving image")
    p.add_argument("--moving-pixel-size-um", type=float, help="CellScape / H&E level-0 pixel size")
    p.add_argument("--work-res-um", type=float, default=1.0, help="resolution for image refinement")
    p.add_argument("--max-shift-um", type=float, default=100.0)
    p.add_argument("--max-rot-deg", type=float, default=10.0)
    p.add_argument("--agree-radius-um", type=float, default=20.0)
    return p.parse_args()


def load_cells(path):
    df = pd.read_csv(path)
    df["lineage"] = df["lineage"].astype(str)
    return df


def main():
    a = parse_args()
    reg_dir = os.path.join(a.project_dir, "data", "registration")
    tab_dir = os.path.join(a.project_dir, "tables")
    fig_dir = os.path.join(a.project_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    xcores = pd.read_csv(os.path.join(tab_dir, "02_cores_all.csv"))
    xcores = xcores[(xcores.slide_id == a.slide_id) & xcores.core_id.notna()].set_index("core_id")
    fix = load_cells(os.path.join(reg_dir, f"xenium_cells_{a.slide_id}.csv.gz"))

    fixed_img = fixed_px = None
    if a.refine == "image" or a.mode == "he":
        if not a.xenium_dir:
            sys.exit("--xenium-dir is required for image refinement / H&E mode")
        base_px = pixel_size_um(a.xenium_dir)
        fixed_img, fixed_px = choose_level(find_morphology_channels(a.xenium_dir)["dapi"], base_px, a.work_res_um)
        print(f"Xenium DAPI at {fixed_px:.3f} um/px for refinement")

    mov_img = None
    if a.refine == "image":
        if not (a.moving_image and a.moving_pixel_size_um):
            sys.exit("--moving-image and --moving-pixel-size-um are required for --refine image")
        probe = PyramidImage(a.moving_image)
        names = probe.channel_names()
        if a.moving_channel in names:
            ch = names.index(a.moving_channel)
        elif a.moving_channel.isdigit():
            ch = int(a.moving_channel)
        else:
            sys.exit(f"channel '{a.moving_channel}' not in {names}")
        mov_img, mov_px = choose_level(a.moving_image, a.moving_pixel_size_um, a.work_res_um, channel=ch)
        print(f"moving image at {mov_px:.3f} um/px (channel {ch})")

    # ===================================================================
    if a.mode == "cellscape":
        mov = load_cells(os.path.join(reg_dir, f"cellscape_cells_{a.slide_id}.csv.gz"))
        ccores = pd.read_csv(os.path.join(tab_dir, "06_cs_cores_all.csv"))
        ccores = ccores[(ccores.slide_id == a.slide_id) & ccores.core_id.notna()].set_index("core_id")

        # ---- 1. global transform: CellScape um -> Xenium um
        if a.landmarks:
            lm = pd.read_csv(a.landmarks)
            src, dst = lm[["moving_x", "moving_y"]].to_numpy(), lm[["fixed_x", "fixed_y"]].to_numpy()
            names = [f"landmark{i}" for i in range(len(lm))]
        else:
            shared = sorted(set(xcores.index) & set(ccores.index))
            src = ccores.loc[shared, ["x_center", "y_center"]].to_numpy()
            dst = xcores.loc[shared, ["x_center", "y_center"]].to_numpy()
            names = shared
            print(f"{len(shared)} cores matched between Xenium and CellScape")
        G, res, keep, model = robust_fit(src, dst, a.global_model)
        det = np.linalg.det(G[:2, :2])
        print(f"global {model}: median residual {np.median(res[keep]):.1f} um, "
              f"{(~keep).sum()} outlier(s) {[n for n, k in zip(names, keep) if not k]}, "
              f"rotation {rot_deg(G):.1f} deg, scale {np.sqrt(abs(det)):.3f}"
              + ("  ** MIRRORED section (det<0) **" if det < 0 else ""))
        if np.median(res[keep]) > 200:
            print("WARNING: large residuals -- core labels probably disagree between 02 and 06 (check both dearray plots).")

        xy_g = apply_tf(G, mov[["x_um", "y_um"]].to_numpy())
        mov["x_global"], mov["y_global"] = xy_g[:, 0], xy_g[:, 1]
        mov["x_aligned"], mov["y_aligned"] = mov["x_global"], mov["y_global"]

        # ---- 2. per-core refinement
        qc_rows, core_tf = [], {}
        for cid in sorted(set(mov.core_id.dropna()) & set(fix.core_id.dropna())):
            mm = (mov.core_id == cid).to_numpy()
            ff = (fix.core_id == cid).to_numpy()
            mxy = mov.loc[mm, ["x_global", "y_global"]].to_numpy()
            fxy = fix.loc[ff, ["x_um", "y_um"]].to_numpy()
            mlab, flab = mov.loc[mm, "lineage"].to_numpy(), fix.loc[ff, "lineage"].to_numpy()
            agree_g, lift_g, nn_g = lineage_agreement(mxy, mlab, fxy, flab, a.agree_radius_um)
            C = None
            if a.refine == "points":
                C = icp_label_aware(mxy, mlab, fxy, flab)
            elif a.refine == "image" and cid in xcores.index:
                r = xcores.loc[cid]
                C = refine_core_image((r.xmin, r.xmax, r.ymin, r.ymax), G, fixed_img, fixed_px, mov_img,
                                      1.0 / a.moving_pixel_size_um, he=False, max_shift_um=a.max_shift_um)
            accepted, lift_r, nn_r, shift, rot = False, np.nan, np.nan, 0.0, 0.0
            if C is not None:
                centre = fxy.mean(0)
                shift = float(np.linalg.norm(apply_tf(C, centre[None])[0] - centre))
                rot = rot_deg(C)
                rxy = apply_tf(C, mxy)
                _, lift_r, nn_r = lineage_agreement(rxy, mlab, fxy, flab, a.agree_radius_um)
                accepted = (shift <= a.max_shift_um and abs(rot) <= a.max_rot_deg
                            and (np.isnan(lift_g) or (not np.isnan(lift_r) and lift_r >= lift_g)))
                if accepted:
                    mov.loc[mm, ["x_aligned", "y_aligned"]] = rxy
                    core_tf[cid] = C.tolist()
            qc_rows.append(dict(core_id=cid, n_cellscape=int(mm.sum()), n_xenium=int(ff.sum()),
                                agreement_global=agree_g, lift_global=lift_g, median_nn_global=nn_g,
                                lift_refined=lift_r, median_nn_refined=nn_r, shift_um=shift,
                                rotation_deg=rot, refinement_accepted=accepted))

        qc = pd.DataFrame(qc_rows)
        qc.to_csv(os.path.join(reg_dir, f"registration_qc_cellscape_{a.slide_id}.csv"), index=False)
        print(qc[["core_id", "lift_global", "lift_refined", "shift_um", "refinement_accepted"]].to_string(index=False))
        print("lift = cross-modal lineage agreement / chance. ~1 = not aligned (or annotations disagree).")

        mov["refine_method"] = np.where(mov.core_id.isin(core_tf.keys()), a.refine, "global")
        mov.to_csv(os.path.join(reg_dir, f"cellscape_aligned_{a.slide_id}.csv.gz"), index=False)
        with open(os.path.join(reg_dir, f"transforms_cellscape_{a.slide_id}.json"), "w") as fh:
            json.dump({"moving": "cellscape_um", "fixed": "xenium_um", "global_model": model,
                       "global": G.tolist(), "cores": core_tf,
                       "global_residuals_um": dict(zip(map(str, names), np.round(res, 1).tolist()))}, fh, indent=1)
        cores_plot = qc.sort_values("lift_global").core_id
        plot_core_grid(os.path.join(fig_dir, f"07_registration_cellscape_{a.slide_id}_global.png"),
                       cores_plot, fix, mov, "x_global", "y_global",
                       f"{a.slide_id}: global transform only (grey = Xenium, colour = CellScape lineage)")
        plot_core_grid(os.path.join(fig_dir, f"07_registration_cellscape_{a.slide_id}_refined.png"),
                       cores_plot, fix, mov, "x_aligned", "y_aligned",
                       f"{a.slide_id}: after per-core refinement ({a.refine})")

    # ===================================================================
    else:  # H&E
        if not a.moving_image:
            sys.exit("--moving-image (the H&E) is required in he mode")
        base_px = pixel_size_um(a.xenium_dir)
        if a.landmarks:
            lm = pd.read_csv(a.landmarks)
            G, res, keep, model = robust_fit(lm[["moving_x", "moving_y"]].to_numpy(),
                                             lm[["fixed_x", "fixed_y"]].to_numpy(), a.global_model)
            print(f"H&E global {model} from landmarks: median residual {np.median(res[keep]):.1f} um")
        elif a.xe_alignment_csv:
            M = np.loadtxt(a.xe_alignment_csv, delimiter=",")
            S = np.diag([base_px, base_px, 1.0])
            # Direction of the Xenium Explorer matrix: pick the candidate that maps the
            # H&E image into the Xenium tissue extent.
            he = PyramidImage(a.moving_image)
            corners = np.array([[0, 0], [he.width, 0], [0, he.height], [he.width, he.height],
                                [he.width / 2, he.height / 2]], float)
            lo, hi = fix[["x_um", "y_um"]].min().to_numpy(), fix[["x_um", "y_um"]].max().to_numpy()
            span = hi - lo

            def score(Gc):
                c = apply_tf(Gc, corners)[-1]
                return np.abs((c - (lo + hi) / 2) / span).sum()
            cands = {"M": S @ M, "inv(M)": S @ np.linalg.inv(M)}
            best = min(cands, key=lambda k: score(cands[k]))
            G, model = cands[best], "xenium_explorer"
            print(f"Xenium Explorer matrix interpreted as {best} (H&E px -> Xenium px), scaled to um")
        else:
            sys.exit("H&E mode needs --landmarks or --xe-alignment-csv for the global transform")

        core_tf, qc_rows = {}, []
        if a.refine == "image":
            he_img, he_px = choose_level(a.moving_image, a.moving_pixel_size_um or 1.0, a.work_res_um)
            for cid, r in xcores.iterrows():
                C = refine_core_image((r.xmin, r.xmax, r.ymin, r.ymax), G, fixed_img, fixed_px, he_img,
                                      1.0, he=True, max_shift_um=a.max_shift_um)
                ok = C is not None and abs(rot_deg(C)) <= a.max_rot_deg
                if ok:
                    core_tf[cid] = C.tolist()
                qc_rows.append(dict(core_id=cid, refinement_accepted=ok,
                                    rotation_deg=rot_deg(C) if C is not None else np.nan))
            pd.DataFrame(qc_rows).to_csv(os.path.join(reg_dir, f"registration_qc_he_{a.slide_id}.csv"), index=False)
        with open(os.path.join(reg_dir, f"transforms_he_{a.slide_id}.json"), "w") as fh:
            json.dump({"moving": "he_px_level0", "fixed": "xenium_um", "global_model": model,
                       "global": G.tolist(), "cores": core_tf}, fh, indent=1)
        print(f"Wrote transforms_he_{a.slide_id}.json ({len(core_tf)} per-core refinements)")


if __name__ == "__main__":
    main()
