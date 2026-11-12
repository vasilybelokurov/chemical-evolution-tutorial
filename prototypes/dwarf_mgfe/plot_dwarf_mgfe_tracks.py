"""All five dwarf [Mg/Fe]-[Fe/H] median tracks on one panel; Fornax extended to low [Fe/H].

LMC, SMC, Sgr, GSE: running median in 0.1-dex [Fe/H] bins (step 0.05), bins with >= 20 stars.
Fornax: same binned median where bins have >= 20 stars; below [Fe/H] = -1.1 (sparse tail) a
  median in sliding windows of K = 10 stars ordered by [Fe/H] (step 2); x = window median [Fe/H]; drawn solid, joined to the core.
Bands: bootstrap 16-84% of the median (300 resamples). Same quality cuts as plot_dwarf_mgfe.py.

Usage (from this folder): ~/Work/venvs/gce/bin/python plot_dwarf_mgfe_tracks.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_dwarf_mgfe import COLORS, FIGS, GALAXIES, GRID, INK, INK2, SURF, load, running_median

K, STEP, TAIL_EDGE = 10, 2, -1.1


def equal_count_median(x, y, k=K, step=STEP, nboot=300, seed=1):
    """Median [Mg/Fe] in sliding windows of k stars ordered by x. Returns xc, med, lo, hi, width."""
    rng = np.random.default_rng(seed)
    o = np.argsort(x)
    x, y = x[o], y[o]
    rows = []
    for i in range(0, len(x) - k + 1, step):
        xs, ys = x[i:i + k], y[i:i + k]
        boot = np.median(rng.choice(ys, (nboot, k)), axis=1)
        rows.append((np.median(xs), np.median(ys), *np.percentile(boot, [16, 84]), xs[-1] - xs[0]))
    return np.array(rows).T


def main():
    fig, ax = plt.subplots(figsize=(8.5, 5.6), facecolor=SURF)
    ax.set_facecolor(SURF)
    for (key, name), col in zip(GALAXIES, COLORS):
        x, y = load(key)
        if key == "fornax":
            # core: same binned median as the other galaxies
            c, m, lo, hi, _ = running_median(x, y)
            ax.fill_between(c, lo, hi, color=col, alpha=0.22, lw=0)
            ax.plot(c, m, color=col, lw=2.2, label=f"{name} (N = {len(x)})")
            # metal-poor tail: equal-count windows of the stars below the core
            tail = x < TAIL_EDGE
            tc, tm, tlo, thi, tw = equal_count_median(x[tail], y[tail])
            tc, tm, tlo, thi = (np.r_[a, b] for a, b in ((tc, c[0]), (tm, m[0]), (tlo, lo[0]), (thi, hi[0])))
            ax.fill_between(tc, tlo, thi, color=col, alpha=0.22, lw=0)
            ax.plot(tc, tm, color=col, lw=2.2)
            print(f"Fornax core [Fe/H] {c.min():+.2f}..{c.max():+.2f}; tail: {tail.sum()} stars, "
                  f"{len(tw)} windows of {K}, centres {tc.min():+.2f}..{tc[-2]:+.2f}, widths {tw.min():.2f}-{tw.max():.2f} dex")
        else:
            c, m, lo, hi, _ = running_median(x, y)
            ax.fill_between(c, lo, hi, color=col, alpha=0.22, lw=0)
            ax.plot(c, m, color=col, lw=2.2, label=f"{name} (N = {len(x)})")
    ax.axhline(0, color=INK2, lw=0.8, zorder=0)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=10)
    ax.set_xlim(-2.3, -0.1)
    ax.set_ylim(-0.3, 0.4)
    ax.set_xlabel("[Fe/H]", color=INK, fontsize=12)
    ax.set_ylabel("median [Mg/Fe]", color=INK, fontsize=12)
    ax.legend(frameon=False, fontsize=10, labelcolor=INK2, loc="upper right")
    ax.set_title("[Mg/Fe] tracks of five dwarf galaxies (APOGEE DR17, H21 selection)",
                 color=INK, fontsize=12, loc="left")
    ax.text(0.02, 0.03, "bands: bootstrap 16–84% of the median\n"
            f"Fornax at [Fe/H] < {TAIL_EDGE}: sliding windows of {K} stars (sparse)",
            transform=ax.transAxes, fontsize=8.5, color=INK2, va="bottom")
    fig.tight_layout()
    out = FIGS / "dwarf_mgfe_tracks.png"
    fig.savefig(out, dpi=200, facecolor=SURF)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
