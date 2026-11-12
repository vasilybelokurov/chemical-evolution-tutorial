"""VICE vs flexCE for the five H21 best fits (H21 Table 4), on the same APOGEE stars.

flexCE tracks come from chemical_evolution's gce_tools.flexce_h21 (NumPy-1 venv), saved to an .npz with keys
<galaxy>_<time|sfr|mgas|feh|sife|mgfe>. VICE tracks are run here (dt = 0.03 Gyr, flexCE uses 30 Myr).

Usage (repo root, GCE venv): python prototypes/h21_vice/compare_vice_flexce.py FLEXCE.npz
"""
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

import h21_vice as H

warnings.filterwarnings("ignore")
fx = np.load(sys.argv[1])
COLORS = {"LMC": "#2a78d6", "SMC": "#eb6834", "Sgr": "#1baf7a", "GSE": "#eda100", "Fornax": "#e87ba4"}
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

fig, ax = plt.subplots(3, 5, figsize=(16, 9.5), facecolor=SURF)
print(f"{'galaxy':7s} {'chi2/N VICE':>11s} {'chi2/N flexCE':>13s}  max|d[Fe/H]|(t>1)  max|d[Si/Fe]| at same [Fe/H]  SFR peak t: VICE flexCE")
for j, (g, p) in enumerate(H.H21_BEST.items()):
    v = H.run(p, dt=0.03)
    f = {k: fx[f"{g}_{k}"] for k in ("time", "sfr", "mgas", "feh", "sife", "mgfe")}
    t = Table.read(H.ROOT / f"data_build/{g.lower()}_h21.ecsv")
    pts, err = np.c_[t["feh"], t["sife"]].astype(float), np.c_[t["feh_err"], t["sife_err"]].astype(float)
    cv, n = H.chi2_track(pts, err, np.c_[v["feh"], v["sife"]])
    cf, _ = H.chi2_track(pts, err, np.c_[f["feh"], f["sife"]])
    late = f["time"] > 1
    dfeh = np.max(np.abs(np.interp(f["time"][late], v["time"], v["feh"]) - f["feh"][late]))
    # [Si/Fe] difference at equal [Fe/H], over the [Fe/H] range both tracks cover, where [Fe/H] rises monotonically
    lo, hi = max(v["feh"][5], f["feh"][5]), min(v["feh"][-1], f["feh"][-1])
    grid = np.linspace(lo, hi, 50)
    dsife = np.max(np.abs(np.interp(grid, v["feh"], v["sife"]) - np.interp(grid, f["feh"], f["sife"])))
    print(f"{g:7s} {cv / n:11.2f} {cf / n:13.2f}  {dfeh:17.3f}  {dsife:28.3f}  "
          f"{v['time'][np.argmax(v['sfr'])]:5.2f} {f['time'][np.argmax(f['sfr'])]:5.2f}")
    ok = (t["feh"] > -9000) & (t["sife"] > -9000)
    a = ax[0, j]
    a.scatter(t["feh"][ok], t["sife"][ok], s=4 if len(t) > 500 else 10, color=COLORS[g], alpha=0.35, lw=0)
    a.plot(v["feh"], v["sife"], color=INK, lw=1.6, label=f"VICE (χ²/N {cv / n:.1f})")
    a.plot(f["feh"], f["sife"], color="#7b2cbf", lw=1.6, ls="--", label=f"flexCE (χ²/N {cf / n:.1f})")
    a.set(xlim=(-2.5, 0.1), ylim=(-0.5, 0.6), xlabel="[Fe/H]")
    a.set_title(g, loc="left", fontsize=11, color=INK); a.legend(frameon=False, fontsize=8, loc="lower left")
    a = ax[1, j]
    a.plot(v["time"], v["sfr"] / v["sfr"].max(), color=INK, lw=1.6)
    a.plot(f["time"], f["sfr"] / f["sfr"].max(), color="#7b2cbf", lw=1.6, ls="--")
    a.set(xlim=(0, 12), ylim=(0, 1.05), xlabel="time [Gyr]")
    a = ax[2, j]
    a.plot(v["time"], v["feh"], color=INK, lw=1.6)
    a.plot(f["time"], f["feh"], color="#7b2cbf", lw=1.6, ls="--")
    a.set(xlim=(0, 12), ylim=(-2.5, 0.1), xlabel="time [Gyr]")
ax[0, 0].set_ylabel("[Si/Fe]"); ax[1, 0].set_ylabel("SFR / peak SFR"); ax[2, 0].set_ylabel("[Fe/H]")
for a in ax.flat:
    a.set_facecolor(SURF); a.grid(color=GRID, lw=0.6)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
fig.suptitle("H21 best fits (Table 4): VICE 1.3.1 (black) vs flexCE 2450aa29 (purple dashed), same yields and stars",
             x=0.01, ha="left", color=INK)
fig.tight_layout()
out = Path(__file__).with_name("figures") / "h21_vice_vs_flexce.png"
fig.savefig(out, dpi=150, facecolor=SURF); print("wrote", out)
