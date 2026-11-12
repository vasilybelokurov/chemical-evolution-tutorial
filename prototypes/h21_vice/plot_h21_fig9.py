"""H21 Figure 9 layout, with the tutorial's VICE runs of the five H21 best fits (Table 4).

Row 1: [Si/Fe] vs [Fe/H], APOGEE stars (H21 Table 2 members) with the model track (black).
Row 2: star-formation rate normalised by its peak, vs time since the model start (dotted: time of peak).
Row 3: [Fe/H] vs time. Red crosses: values H21 quote in §5.3 for their flexCE fits.

Usage (repo root, GCE venv): python prototypes/h21_vice/plot_h21_fig9.py
"""
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

import h21_vice as H

warnings.filterwarnings("ignore")
COLORS = {"LMC": "#2a78d6", "SMC": "#eb6834", "Sgr": "#1baf7a", "GSE": "#eda100", "Fornax": "#e87ba4"}
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
# H21 §5.3, flexCE: (time [Gyr], [Fe/H]) pairs they quote
H21_QUOTES = {"LMC": [(5, -1.2), (10, -0.7), (9.5, -0.8), (12, -0.3)], "SMC": [(5, -1.5)],
              "GSE": [(5, -0.5), (4, -0.7)]}

fig, ax = plt.subplots(3, 5, figsize=(16, 9.5), facecolor=SURF)
for j, (g, p) in enumerate(H.H21_BEST.items()):
    tr = H.run(p, dt=0.03)
    t = Table.read(H.ROOT / f"data_build/{g.lower()}_h21.ecsv")
    ok = (t["feh"] > -9000) & (t["sife"] > -9000)
    chi, n = H.chi2_track(np.c_[t["feh"], t["sife"]].astype(float), np.c_[t["feh_err"], t["sife_err"]].astype(float),
                          np.c_[tr["feh"], tr["sife"]])
    a = ax[0, j]
    a.scatter(t["feh"][ok], t["sife"][ok], s=4 if len(t) > 500 else 10, color=COLORS[g], alpha=0.4, lw=0)
    a.plot(tr["feh"], tr["sife"], color=INK, lw=1.6)
    a.set(xlim=(-2.5, 0.1), ylim=(-0.5, 0.6), xlabel="[Fe/H]")
    a.set_title(f"{g}  (N = {n}, χ²/N = {chi / n:.1f})", loc="left", fontsize=10, color=INK)
    ipk = np.argmax(tr["sfr"])
    a = ax[1, j]
    a.plot(tr["time"], tr["sfr"] / tr["sfr"][ipk], color=COLORS[g], lw=2)
    a.axvline(tr["time"][ipk], color=INK2, ls=":", lw=1)
    a.set(xlim=(0, 12), ylim=(0, 1.05), xlabel="time since start [Gyr]")
    a = ax[2, j]
    a.plot(tr["time"], tr["feh"], color=COLORS[g], lw=2)
    a.axhline(tr["feh"][ipk], color=INK2, ls=":", lw=1)
    for tq, fq in H21_QUOTES.get(g, []):
        a.plot(tq, fq, "x", color="#c0392b", ms=8, mew=2)
    a.set(xlim=(0, 12), ylim=(-2.5, 0.1), xlabel="time since start [Gyr]")
ax[0, 0].set_ylabel("[Si/Fe]"); ax[1, 0].set_ylabel("SFR / peak SFR"); ax[2, 0].set_ylabel("[Fe/H]")
for a in ax.flat:
    a.set_facecolor(SURF); a.grid(color=GRID, lw=0.6)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
ax[2, 0].plot([], [], "x", color="#c0392b", label="H21 §5.3 quoted (flexCE)"); ax[2, 0].legend(frameon=False, fontsize=8)
fig.suptitle("H21 best-fit models (Table 4) run with VICE 1.3.1 and H21-matched flexCE yields — layout of H21 Fig. 9",
             x=0.01, ha="left", color=INK)
fig.tight_layout()
out = Path(__file__).with_name("figures"); out.mkdir(exist_ok=True)
fig.savefig(out / "h21_fig9_vice.png", dpi=150, facecolor=SURF)
print("wrote", out / "h21_fig9_vice.png")
