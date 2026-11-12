"""Tracks of H21's Table-4 models vs our grid best fits (VICE and flexCE), on the same APOGEE stars.

All tracks are taken from the full-grid runs (no new model runs):
  red      H21's published best-fit parameters, run in VICE
  black    our VICE best fit (lowest chi2 in the grid)
  purple   our flexCE best fit
Rows: [Si/Fe] vs [Fe/H] with stars; SFR(t) normalised to the H21 model's peak; [Fe/H](t).

Usage (repo root, GCE venv): python prototypes/h21_vice/plot_h21_vs_ours.py
Reads grid_runs/h21_grid_{vice,flexce}.npz and h21_grid_scores.csv; writes figures/h21_vs_ours.png.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.table import Table

import h21_vice as H

HERE = Path(__file__).parent
RUNS = HERE / "grid_runs"
KEYS = ("sfe0", "eta", "tau_b", "sigma_b", "f_b")
COLORS = {"LMC": "#2a78d6", "SMC": "#eb6834", "Sgr": "#1baf7a", "GSE": "#eda100", "Fornax": "#e87ba4"}
INK, GRID, SURF, RED, PURPLE = "#0b0b0b", "#e4e3df", "#fcfcfb", "#d62828", "#7b2cbf"

z = {e: np.load(RUNS / f"h21_grid_{e}.npz") for e in ("vice", "flexce")}
df = pd.read_csv(RUNS / "h21_grid_scores.csv")


def index_of(g, p):
    """Row in the grid file for galaxy g and parameters p (None = unused)."""
    zz = z["vice"]
    m = zz["galaxy"] == g
    for k in KEYS:
        m &= np.isnan(zz[k]) if p[k] is None else np.isclose(zz[k], p[k])
    (i,) = np.flatnonzero(m)
    return i


def label(row):
    return ", ".join(f"{k} {row[k]:g}" for k in KEYS if np.isfinite(row[k]))


fig, ax = plt.subplots(3, 5, figsize=(17, 10), facecolor=SURF)
for j, g in enumerate(H.H21_BEST):
    t = Table.read(H.ROOT / f"data_build/{g.lower()}_h21.ecsv")
    ok = (t["feh"] > -9000) & (t["sife"] > -9000)
    i_h21 = index_of(g, H.H21_BEST[g])
    s = df[df.galaxy == g]
    n = s.n.iloc[0]
    v = s[s.engine == "vice"].reset_index(drop=True)          # same row order as the grid file block
    off = np.flatnonzero(z["vice"]["galaxy"] == g)[0]
    curves = [("H21 Table 4", "vice", i_h21, RED, "-", v.chi2[i_h21 - off] / n)]
    for e, c, ls in (("vice", INK, "-"), ("flexce", PURPLE, "--")):
        se = s[s.engine == e].reset_index(drop=True)
        b = int(se.chi2.idxmin())
        curves.append((f"our {'VICE' if e == 'vice' else 'flexCE'} best", e, off + b, c, ls, se.chi2[b] / n))
    print(f"{g:7s} H21: {label(v.iloc[i_h21 - off])}  |  " +
          "  |  ".join(f"{name}: {label(s[s.engine == e].reset_index(drop=True).iloc[i - off])}"
                       for name, e, i, *_ in curves[1:]))

    a = ax[0, j]
    a.scatter(t["feh"][ok], t["sife"][ok], s=4 if len(t) > 500 else 10, color=COLORS[g], alpha=0.3, lw=0)
    peak = z["vice"]["sfr"][i_h21].max()
    for name, e, i, c, ls, chi in curves:
        zz = z[e]
        ax[0, j].plot(zz["feh"][i], zz["sife"][i], color=c, ls=ls, lw=1.6, label=f"{name} (χ²/N {chi:.2f})")
        ax[1, j].plot(zz["time"], zz["sfr"][i] / peak, color=c, ls=ls, lw=1.6)
        ax[2, j].plot(zz["time"], zz["feh"][i], color=c, ls=ls, lw=1.6)
    a.set(xlim=(-2.5, 0.1), ylim=(-0.5, 0.6), xlabel="[Fe/H]")
    a.set_title(g, loc="left", fontsize=11, color=INK)
    a.legend(frameon=False, fontsize=7.5, loc="lower left")
    ax[1, j].set(xlim=(0, 12), ylim=(0, None), xlabel="time [Gyr]")
    ax[2, j].set(xlim=(0, 12), ylim=(-2.5, 0.1), xlabel="time [Gyr]")
ax[0, 0].set_ylabel("[Si/Fe]"); ax[1, 0].set_ylabel("SFR / peak SFR of H21 model"); ax[2, 0].set_ylabel("[Fe/H]")
for a in ax.flat:
    a.set_facecolor(SURF); a.grid(color=GRID, lw=0.6)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
fig.suptitle("H21's best-fit parameters (red, run in VICE) vs our grid best fits: VICE (black), flexCE (purple dashed)",
             x=0.01, ha="left", color=INK)
fig.tight_layout()
out = HERE / "figures/h21_vs_ours.png"
fig.savefig(out, dpi=140, facecolor=SURF)
print("wrote", out)
