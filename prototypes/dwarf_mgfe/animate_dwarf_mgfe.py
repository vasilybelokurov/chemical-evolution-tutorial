"""Animate the five dwarf [Mg/Fe]-[Fe/H] tracks growing from low to high [Fe/H].

A sweep in [Fe/H] (a monotonic proxy for time within each galaxy, not a common clock) reveals each
galaxy's median track with its bootstrap 16-84% band; a marker shows the
current end of each track. No grid, legend, star points or in-plot text
(colours: LMC blue, SMC orange, Sgr aqua, GSE yellow, Fornax pink). Tracks are the same as in plot_dwarf_mgfe_tracks.py (Fornax: binned
median for [Fe/H] >= -1.1, 10-star sliding windows below, drawn as one solid line).

Usage (from this folder): ~/Work/venvs/gce/bin/python animate_dwarf_mgfe.py   -> dwarf_mgfe_tracks.mp4 and .gif
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter

from plot_dwarf_mgfe import COLORS, FIGS, GALAXIES, GRID, INK, INK2, SURF, load, running_median
from plot_dwarf_mgfe_tracks import TAIL_EDGE, equal_count_median

X0, X1, DX, HOLD, FPS = -2.05, -0.15, 0.01, 40, 20
HERE = FIGS


def tracks():
    """Per galaxy: stars (x, y) and track arrays (c, m, lo, hi) sorted by c."""
    out = {}
    for key, name in GALAXIES:
        x, y = load(key)
        c, m, lo, hi, _ = running_median(x, y)
        if key == "fornax":
            t = x < TAIL_EDGE
            tc, tm, tlo, thi, _ = equal_count_median(x[t], y[t])
            c, m, lo, hi = (np.r_[a, b] for a, b in ((tc, c), (tm, m), (tlo, lo), (thi, hi)))
        out[key] = dict(name=name, x=x, y=y, c=c, m=m, lo=lo, hi=hi)
    return out


def upto(tr, xs):
    """Track points with c <= xs, plus a linearly interpolated end point at xs (smooth growth)."""
    c = tr["c"]
    if xs < c[0]:
        return None
    k = np.searchsorted(c, xs, side="right")
    seg = {q: tr[q][:k] for q in ("c", "m", "lo", "hi")}
    if k < len(c):
        f = (xs - c[k - 1]) / (c[k] - c[k - 1])
        for q in ("m", "lo", "hi"):
            seg[q] = np.r_[seg[q], tr[q][k - 1] + f * (tr[q][k] - tr[q][k - 1])]
        seg["c"] = np.r_[seg["c"], xs]
    return seg


def main():
    T = tracks()
    fig, ax = plt.subplots(figsize=(9, 5.9), facecolor=SURF)
    ax.set_facecolor(SURF)
    ax.axhline(0, color=INK2, lw=0.8, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=10)
    ax.set_xlim(-2.3, -0.1)
    ax.set_ylim(-0.35, 0.45)
    ax.set_xlabel("[Fe/H]", color=INK, fontsize=12)
    ax.set_ylabel("[Mg/Fe]", color=INK, fontsize=12)
    ax.set_title("Five dwarf galaxies building up their chemistry (APOGEE DR17, H21 selection)",
                 color=INK, fontsize=12, loc="left")

    art = {}
    for (key, _), col in zip(GALAXIES, COLORS):
        tr = T[key]
        art[key] = dict(
            ring=ax.plot([], [], color=SURF, lw=5, zorder=3)[0],
            line=ax.plot([], [], color=col, lw=2.4, zorder=4,
                         label=f"{tr['name']} (N = {len(tr['x'])})")[0],
            head=ax.plot([], [], "o", ms=9, color=col, mec=SURF, mew=2, zorder=5)[0],
            band=None)
    sweep = ax.axvline(X0, color=INK2, lw=1, ls=":", zorder=2)

    xs_all = np.r_[np.arange(X0, X1 + DX / 2, DX), np.full(HOLD, X1)]

    def frame(i):
        xs = xs_all[i]
        sweep.set_xdata([xs, xs])
        for key, a in art.items():
            tr = T[key]
            if a["band"] is not None:
                a["band"].remove()
                a["band"] = None
            seg = upto(tr, xs)
            if seg is None:
                for h in ("ring", "line", "head"):
                    a[h].set_data([], [])
                continue
            a["band"] = ax.fill_between(seg["c"], seg["lo"], seg["hi"], color=a["line"].get_color(),
                                        alpha=0.22, lw=0, zorder=2)
            a["ring"].set_data(seg["c"], seg["m"])
            a["line"].set_data(seg["c"], seg["m"])
            a["head"].set_data([seg["c"][-1]], [seg["m"][-1]])
        return []

    anim = FuncAnimation(fig, frame, frames=len(xs_all), interval=1000 / FPS)
    fig.tight_layout()
    mp4, gif = HERE / "dwarf_mgfe_tracks.mp4", HERE / "dwarf_mgfe_tracks.gif"
    anim.save(mp4, writer=FFMpegWriter(fps=FPS, bitrate=2400), dpi=150, savefig_kwargs={"facecolor": SURF})
    anim.save(gif, writer=PillowWriter(fps=FPS), dpi=80, savefig_kwargs={"facecolor": SURF})
    frame(len(xs_all) - 1)
    fig.savefig(HERE / "dwarf_mgfe_tracks_lastframe.png", dpi=150, facecolor=SURF)
    frame(int((-1.2 - X0) / DX))
    fig.savefig(HERE / "dwarf_mgfe_tracks_midframe.png", dpi=150, facecolor=SURF)
    print(f"wrote {mp4} and {gif}: {len(xs_all)} frames at {FPS} fps = {len(xs_all) / FPS:.1f} s")


if __name__ == "__main__":
    main()
