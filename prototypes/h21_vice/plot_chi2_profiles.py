"""chi2 profiles of the full H21 grid, VICE vs flexCE, with H21's Table-4 values marked.

For each galaxy and each varied parameter: the lowest chi2 over all other parameters, minus the overall lowest.
Best-fit chi2/N is 4.7-8.5, i.e. the quoted errors are too small for a one-zone track, so the plotted quantity is
    Δχ²_scaled = Δχ² / (χ²_min / N)
(errors inflated so the best model has χ²/N = 1). Δχ²_scaled = 1 and 4 are drawn as rough 1σ / 2σ guides.

Usage (repo root, GCE venv): python prototypes/h21_vice/plot_chi2_profiles.py
Reads grid_runs/h21_grid_scores.csv; writes figures/h21_chi2_profiles.png.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import h21_vice as H

HERE = Path(__file__).parent
LABEL = {"sfe0": "SFE₀ [Gyr⁻¹]", "eta": "η (outflow)", "tau_b": "τ_b [Gyr]", "sigma_b": "σ_b [Gyr]", "f_b": "f_b (burst)"}
ENGINE = {"vice": ("#0b0b0b", "-", "VICE"), "flexce": ("#7b2cbf", "--", "flexCE")}
INK, GRID, SURF, RED = "#0b0b0b", "#e4e3df", "#fcfcfb", "#d62828"

df = pd.read_csv(HERE / "grid_runs/h21_grid_scores.csv")
fig, axes = plt.subplots(5, 4, figsize=(14, 14), facecolor=SURF)
for row, g in enumerate(H.H21_BEST):
    params = list(H.H21_GRID[g])
    for col in range(4):
        a = axes[row, col]
        if col >= len(params):
            a.axis("off")
            continue
        k = params[col]
        for e, (c, ls, lab) in ENGINE.items():
            s = df[(df.engine == e) & (df.galaxy == g)]
            best = s.chi2.min()
            scale = best / s.n.iloc[0]
            prof = s.groupby(k).chi2.min()
            a.plot(prof.index, (prof - best) / scale, color=c, ls=ls, marker="o", ms=3, lw=1.5, label=lab)
        a.axvline(H.H21_BEST[g][k], color=RED, lw=1.2, alpha=0.8, label="H21 Table 4")
        for y in (1, 4):
            a.axhline(y, color=INK, lw=0.6, ls=":")
        a.set_yscale("symlog", linthresh=1)
        a.set_ylim(-0.1, 1e3)
        a.set_xlabel(LABEL[k])
        if col == 0:
            a.set_ylabel(f"{g}\nΔχ² (scaled)")
        a.set_facecolor(SURF)
        a.grid(color=GRID, lw=0.6)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
axes[0, 0].legend(frameon=False, fontsize=8)
fig.suptitle("H21 grid: χ² profiles (min over other parameters), errors scaled so best χ²/N = 1; "
             "dotted: Δχ² = 1, 4", x=0.01, ha="left", color=INK)
fig.tight_layout()
out = HERE / "figures/h21_chi2_profiles.png"
fig.savefig(out, dpi=130, facecolor=SURF)
print("wrote", out)

# Numbers behind the figure: scaled Δχ² of H21's own model, per engine.
for g in H.H21_BEST:
    line = []
    for e in ENGINE:
        s = df[(df.engine == e) & (df.galaxy == g)]
        m = np.ones(len(s), bool)
        for k in H.H21_GRID[g]:
            m &= np.isclose(s[k], H.H21_BEST[g][k])
        best = s.chi2.min()
        line.append(f"{e} {float((s.chi2[m].iloc[0] - best) / (best / s.n.iloc[0])):6.1f}")
    print(f"{g:7s} scaled Δχ² of H21 model: " + "  ".join(line))
