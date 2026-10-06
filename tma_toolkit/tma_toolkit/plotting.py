"""Plots for TMA cores (matplotlib).

Colour rules: categories get 8 fixed hues in a fixed order; a 9th+ category is
folded into "Other" (grey) instead of generating more hues. Magnitudes use one
blue ramp (light -> dark), signed values a blue <-> grey <-> red diverging ramp.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
OTHER = "#b4b2a9"
SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
DIVERGING = ["#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f4a3a3", "#e34948", "#a52828"]


def _cmap(colors, name):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name, colors)


def category_colors(values, max_categories: int = 8, order=None) -> dict:
    """Map categories -> colours. The most frequent 7 (or all, if <= 8) keep a hue;
    the rest become 'Other'. Pass `order` to keep colours stable across plots."""
    s = pd.Series(values).astype(str)
    cats = list(order) if order is not None else s.value_counts().index.tolist()
    if len(cats) > max_categories:
        keep = cats[: max_categories - 1]
        cmap = {c: CATEGORICAL[i] for i, c in enumerate(keep)}
        cmap["Other"] = OTHER
        return cmap
    return {c: CATEGORICAL[i] for i, c in enumerate(cats)}


def _fold(values, cmap):
    s = pd.Series(values).astype(str)
    return s.where(s.isin(cmap.keys()), "Other").to_numpy()


def plot_cores(adata, cores: list[str], color: str, core_key="core_id", ncols: int = 5,
               panel_size: float = 3.0, point_size: float | None = None, order=None, title=None):
    """Small-multiple gallery: one panel per core, cells coloured by an obs column."""
    import matplotlib.pyplot as plt
    cores = list(cores)
    n = len(cores)
    nr = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nr, min(ncols, n), figsize=(panel_size * min(ncols, n), panel_size * nr + 0.6),
                             squeeze=False)
    obs = adata.obs
    sel = obs[core_key].isin(cores).to_numpy()
    vals = obs.loc[sel, color]
    numeric = pd.api.types.is_numeric_dtype(vals) and not pd.api.types.is_bool_dtype(vals)
    if not numeric:
        cmap = category_colors(vals, order=order)
        vmin = vmax = None
    else:
        vmin, vmax = np.nanquantile(vals, [0.01, 0.99])
    for ax, cid in zip(axes.ravel(), cores):
        m = (obs[core_key] == cid).to_numpy()
        xy = adata.obsm["spatial"][m]
        s = point_size or max(0.5, min(6.0, 4000 / max(m.sum(), 1)))
        if numeric:
            sc = ax.scatter(xy[:, 0], xy[:, 1], c=obs.loc[m, color], s=s, cmap=_cmap(SEQ_BLUE, "seq"),
                            vmin=vmin, vmax=vmax, linewidths=0)
        else:
            lab = _fold(obs.loc[m, color], cmap)
            ax.scatter(xy[:, 0], xy[:, 1], c=[cmap[v] for v in lab], s=s, linewidths=0)
        ax.set_title(cid, fontsize=8)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.axis("off")
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    if numeric:
        fig.colorbar(sc, ax=axes.ravel().tolist(), shrink=0.6, label=color)
    else:
        handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=k) for k, c in cmap.items()]
        fig.legend(handles=handles, loc="lower center", ncol=min(8, len(handles)), fontsize=8, frameon=False)
    fig.suptitle(title or color, fontsize=10)
    return fig


def plot_composition(comp: pd.DataFrame, cores: pd.DataFrame, group_col: str | None = None, title=None):
    """Stacked bars, one per core, grouped by a metadata column. At most 8 colours:
    beyond that, the least abundant cell types are summed into 'Other'."""
    import matplotlib.pyplot as plt
    comp = comp.copy()
    order = comp.mean().sort_values(ascending=False).index.tolist()
    if len(order) > 8:
        keep = order[:7]
        comp["Other"] = comp[order[7:]].sum(1)
        comp = comp[keep + ["Other"]]
        order = keep + ["Other"]
    colors = {c: (OTHER if c == "Other" else CATEGORICAL[i]) for i, c in enumerate(order)}
    if group_col:
        g = cores[group_col].reindex(comp.index).astype(str)
        comp = comp.assign(_g=g.to_numpy()).sort_values("_g")
        groups = comp.pop("_g")
    fig, ax = plt.subplots(figsize=(max(6, 0.12 * len(comp) + 2), 4))
    bottom = np.zeros(len(comp))
    x = np.arange(len(comp))
    for c in order:
        ax.bar(x, comp[c], bottom=bottom, color=colors[c], width=0.9, edgecolor="white", linewidth=0.3, label=c)
        bottom += comp[c].to_numpy()
    if group_col:
        starts = np.flatnonzero(np.r_[True, groups.to_numpy()[1:] != groups.to_numpy()[:-1]])
        for s, e in zip(starts, np.r_[starts[1:], len(groups)]):
            ax.text((s + e - 1) / 2, 1.02, groups.iloc[s], ha="center", va="bottom", fontsize=8)
            if s > 0:
                ax.axvline(s - 0.5, color="#555", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(comp.index, rotation=90, fontsize=5 if len(comp) > 60 else 7)
    ax.set_ylabel("fraction of cells")
    ax.set_ylim(0, 1.08)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
    ax.set_title(title or "Composition per core", fontsize=10)
    fig.tight_layout()
    return fig


def plot_similarity(corr: pd.DataFrame, cores: pd.DataFrame, annotate=("organ", "patient_id"), title=None):
    """Clustered core x core correlation heatmap with metadata strips."""
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform
    d = 1 - corr.fillna(0).to_numpy()
    np.fill_diagonal(d, 0)
    order = leaves_list(linkage(squareform(d, checks=False), "average")) if len(corr) > 2 else np.arange(len(corr))
    c = corr.iloc[order, order]
    ann = [a for a in annotate if a in cores and cores[a].reindex(c.index).nunique() <= 8]
    skipped = [a for a in annotate if a in cores and a not in ann]
    if skipped:
        print(f"plot_similarity: no colour strip for {skipped} (> 8 values); use similarity_by_relation() instead")
    fig = plt.figure(figsize=(9, 8))
    h = 0.03
    ax = fig.add_axes([0.12, 0.1, 0.7, 0.7])
    im = ax.imshow(c.to_numpy(), cmap=_cmap(SEQ_BLUE, "seq"), aspect="auto")
    ax.set_xticks([]); ax.set_yticks([])
    for i, a in enumerate(ann):
        vals = cores[a].reindex(c.index).astype(str)
        cmap = category_colors(vals)
        rgb = [[cmap.get(v, OTHER) for v in _fold(vals, cmap)]]
        from matplotlib.colors import to_rgb
        axs = fig.add_axes([0.12, 0.81 + i * (h + 0.005), 0.7, h])
        axs.imshow([[to_rgb(x) for x in rgb[0]]], aspect="auto")
        axs.set_xticks([]); axs.set_yticks([])
        axs.set_ylabel(a, rotation=0, ha="right", va="center", fontsize=8)
        if len(cmap) <= 8:
            axs.legend(handles=[plt.Line2D([], [], marker="s", ls="", color=v, label=k) for k, v in cmap.items()],
                       bbox_to_anchor=(1.01, 0.5), loc="center left", fontsize=6, frameon=False, ncol=4)
    fig.colorbar(im, cax=fig.add_axes([0.84, 0.1, 0.02, 0.3]), label="correlation")
    fig.suptitle(title or "Core similarity (pseudobulk)", fontsize=10)
    return fig


def plot_group_values(values: pd.DataFrame, cores: pd.DataFrame, group_col: str, columns: list[str],
                      ncols: int = 4, ylabel="value", title=None):
    """One panel per column (cell type / gene): box + points per group, one point per core."""
    import matplotlib.pyplot as plt
    d = values.join(cores[[group_col]], how="inner")
    groups = sorted(d[group_col].dropna().astype(str).unique())
    cmap = category_colors(groups, order=groups)
    n = len(columns)
    nr = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nr, min(ncols, n), figsize=(3.2 * min(ncols, n), 2.8 * nr), squeeze=False)
    rng = np.random.default_rng(0)
    for ax, col in zip(axes.ravel(), columns):
        data = [d.loc[d[group_col].astype(str) == g, col].dropna().to_numpy() for g in groups]
        ax.boxplot(data, widths=0.5, showfliers=False, medianprops={"color": "#333"})
        for i, (g, v) in enumerate(zip(groups, data)):
            c = cmap.get(g, OTHER)
            ax.scatter(np.full(len(v), i + 1) + rng.uniform(-0.15, 0.15, len(v)), v, s=10, color=c, linewidths=0)
        ax.set_xticks(range(1, len(groups) + 1))
        ax.set_xticklabels(groups, rotation=45, ha="right", fontsize=7)
        ax.set_title(col, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    axes[0, 0].set_ylabel(ylabel)
    fig.suptitle(title or f"by {group_col} (one point per core)", fontsize=10)
    fig.tight_layout()
    return fig


def plot_dearray(adata, cores: pd.DataFrame, slide_id: str, max_cells: int = 300_000):
    """Cells of one slide + each detected core labelled 'core_id (row col, #detected)'.
    THE plot to compare against the physical TMA map."""
    import matplotlib.pyplot as plt
    m = (adata.obs["slide_id"] == slide_id).to_numpy()
    xy = adata.obsm["spatial"][m]
    assigned = adata.obs.loc[m, "core_id"].notna().to_numpy()
    if len(xy) > max_cells:
        idx = np.random.default_rng(0).choice(len(xy), max_cells, replace=False)
        xy, assigned = xy[idx], assigned[idx]
    w, h = np.ptp(xy[:, 0]), np.ptp(xy[:, 1])
    fig, ax = plt.subplots(figsize=(min(24, 4 + w / 800), min(18, 3 + h / 800)))
    s = float(np.clip(2e5 / max(len(xy), 1), 0.2, 3.0))     # visible for 10k and 3M cells alike
    ax.scatter(xy[~assigned, 0], xy[~assigned, 1], s=s, c=CATEGORICAL[7], linewidths=0, label="not in a core")
    ax.scatter(xy[assigned, 0], xy[assigned, 1], s=s, c=OTHER, linewidths=0, label="in a core")
    for r in cores.itertuples(index=False):
        lab = f"{r.core_id if isinstance(r.core_id, str) else '??'}\n{r.core_row}{r.core_col} #{r.detected_core}"
        ax.text(r.x_center, r.y_center, lab, ha="center", va="center", fontsize=6, fontweight="bold",
                color="#000" if isinstance(r.core_id, str) else CATEGORICAL[7])
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_title(f"{slide_id}: core assignment -- compare every label with the TMA map "
                 f"(rotation corrected: {cores['grid_rotation_deg'].iloc[0]:.1f} deg)", fontsize=9)
    ax.legend(loc="upper right", markerscale=20, fontsize=7, frameon=False)
    return fig
