"""[Mg/Fe] vs [Fe/H] for the five Hasselquist+21 dwarfs (APOGEE DR17), from data_build/ (cerium_bump extracts).

Panels 1-5: one galaxy each (its stars coloured, the other four in light grey) with its running median.
Panel 6: all five running medians with bootstrap 16-84% bands of the median.
Quality: MG_FE_FLAG == 0 and FE_H_FLAG == 0, finite values (samples already carry the H21 selection).
Running median: 0.1-dex [Fe/H] bins, step 0.05 dex, shown where a bin has >= 20 stars.

Usage (from this folder): ~/Work/venvs/gce/bin/python plot_dwarf_mgfe.py [--out dwarf_mgfe_feh.png]
"""
import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

# H21 Table-2 members x APOGEE DR17 (copied from cerium_bump/results/samples; see data/MANIFEST.json).
SAMPLES = Path(os.environ.get("CE_SAMPLES", Path(__file__).resolve().parents[2] / "data_build"))
FIGS = Path(__file__).with_name("figures")
GALAXIES = [("lmc", "LMC"), ("smc", "SMC"), ("sgr", "Sgr"), ("gse", "GSE"), ("fornax", "Fornax")]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]  # validated categorical slots 1-5
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
XLIM, YLIM = (-2.3, 0.0), (-0.45, 0.55)


def load(key):
    t = Table.read(SAMPLES / f"{key}_h21.ecsv")
    ok = (t["mgfe_flag"] == 0) & (t["feh_flag"] == 0) & np.isfinite(t["mgfe"]) & np.isfinite(t["feh"])
    return np.asarray(t["feh"][ok], float), np.asarray(t["mgfe"][ok], float)


def running_median(x, y, width=0.1, step=0.05, nmin=20, nboot=300, seed=1):
    """Median of y in sliding x bins; bootstrap 16/84% of the median. Returns centres, med, lo, hi, n."""
    rng = np.random.default_rng(seed)
    centres = np.arange(x.min(), x.max() + step, step)
    rows = []
    for c in centres:
        s = y[np.abs(x - c) <= width / 2]
        if len(s) < nmin:
            continue
        boot = np.median(rng.choice(s, (nboot, len(s))), axis=1)
        rows.append((c, np.median(s), *np.percentile(boot, [16, 84]), len(s)))
    return np.array(rows).T if rows else np.empty((5, 0))


def style(ax):
    ax.set_facecolor(SURF)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(FIGS / "dwarf_mgfe_feh.png"))
    out = Path(ap.parse_args().out)

    data = {k: load(k) for k, _ in GALAXIES}
    meds = {k: running_median(*data[k]) for k in data}

    fig, axes = plt.subplots(2, 3, figsize=(13, 7.6), sharex=True, sharey=True, facecolor=SURF)
    axes = axes.ravel()
    for i, ((key, name), col) in enumerate(zip(GALAXIES, COLORS)):
        ax = axes[i]
        style(ax)
        for other, _ in GALAXIES:
            if other != key:
                ax.scatter(*data[other], s=3, c="#c9c8c3", lw=0, rasterized=True)
        x, y = data[key]
        ax.scatter(x, y, s=5 if len(x) > 1000 else 9, c=col, alpha=0.55, lw=0, rasterized=True)
        c, m, lo, hi, _ = meds[key]
        ax.plot(c, m, color=SURF, lw=4.5)  # surface ring under the median line
        ax.plot(c, m, color=INK, lw=2)
        ax.text(0.97, 0.95, f"{name}\nN = {len(x)}", transform=ax.transAxes, ha="right", va="top",
                fontsize=11, color=INK, fontweight="bold" if False else "normal")
    ax = axes[5]
    style(ax)
    for (key, name), col in zip(GALAXIES, COLORS):
        c, m, lo, hi, _ = meds[key]
        ax.fill_between(c, lo, hi, color=col, alpha=0.25, lw=0)
        ax.plot(c, m, color=col, lw=2, label=name)
    ax.legend(frameon=False, fontsize=9, loc="upper right", labelcolor=INK2)
    ax.text(0.03, 0.04, "running medians (0.1 dex bins, N ≥ 20)\nbands: bootstrap 16–84% of median",
            transform=ax.transAxes, fontsize=8, color=INK2, va="bottom")
    for ax in axes[3:]:
        ax.set_xlabel("[Fe/H]", color=INK, fontsize=11)
    for ax in axes[::3]:
        ax.set_ylabel("[Mg/Fe]", color=INK, fontsize=11)
    fig.suptitle("[Mg/Fe] vs [Fe/H] in five dwarf galaxies — APOGEE DR17, Hasselquist+21 selection",
                 color=INK, fontsize=12, x=0.07, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(out, dpi=200, facecolor=SURF)
    print(f"wrote {out}")
    for key, name in GALAXIES:
        c, m, *_ , n = meds[key]
        print(f"{name:7s} N={len(data[key][0]):5d}  [Fe/H] {c.min():+.2f}..{c.max():+.2f}  "
              f"median [Mg/Fe] {m.max():+.3f} -> {m.min():+.3f}")


if __name__ == "__main__":
    main()
