"""Score the H21 grid runs against the APOGEE stars and report the best fits (VICE, flexCE, H21 Table 4).

chi2 = chi2_track (every star's distance to the [Si/Fe]-[Fe/H] track, in units of its error in that direction)
     + chi2_endpoint (H21's penalty for evolving past the most metal-rich stars; not for GSE).
With --check-cerium-bump, the VICE scores are compared model by model with cerium_bump's grid.

Usage (repo root, GCE venv):
    python prototypes/h21_vice/score_grid.py [--check-cerium-bump]
Writes grid_runs/h21_grid_scores.csv.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

import h21_vice as H

RUNS = Path(__file__).with_name("grid_runs")
KEYS = ("sfe0", "eta", "tau_b", "sigma_b", "f_b")
VARIED = {g: tuple(H.H21_GRID[g]) for g in H.H21_GRID}


def stars(g):
    t = Table.read(H.ROOT / f"data_build/{g.lower()}_h21.ecsv")
    return np.c_[t["feh"], t["sife"]].astype(float), np.c_[t["feh_err"], t["sife_err"]].astype(float)


def score(engine):
    z = np.load(RUNS / f"h21_grid_{engine}.npz")
    rows, data = [], {g: stars(g) for g in H.H21_BEST}
    for i, g in enumerate(z["galaxy"]):
        pts, err = data[g]
        track = np.c_[z["feh"][i], z["sife"][i]]
        chi, n = H.chi2_track(pts, err, track)
        pen = H.chi2_endpoint(pts, err, track, z["time"], g, dt=0.03)
        rows.append(dict(engine=engine, galaxy=g, **{k: z[k][i] for k in KEYS}, n=n,
                         chi2_track=chi, chi2_endpoint=pen, chi2=chi + pen))
    return pd.DataFrame(rows)


def same(a, b):
    """Parameter match, treating NaN == NaN (unused parameters)."""
    return np.all([(np.isnan(a[k]) and np.isnan(b[k])) or np.isclose(a[k], b[k]) for k in KEYS])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-cerium-bump", action="store_true")
    args = ap.parse_args()
    df = pd.concat([score(e) for e in ("vice", "flexce") if (RUNS / f"h21_grid_{e}.npz").exists()])
    df.to_csv(RUNS / "h21_grid_scores.csv", index=False)

    if args.check_cerium_bump:
        cb = pd.read_csv(Path.home() / "Work/Code/cerium_bump/results/grids/h21_matched_dt0.03.csv")
        v = df[df.engine == "vice"].reset_index(drop=True)
        worst = 0.0
        for _, r in v.iterrows():
            m = cb[cb.galaxy == r.galaxy]
            for k in KEYS:
                m = m[np.isnan(m[k]) & np.isnan(r[k]) | np.isclose(m[k], r[k])]
            assert len(m) == 1, (r.galaxy, dict(r))
            worst = max(worst, abs(r.chi2_track / m.chi2_track.iloc[0] - 1),
                        abs(r.chi2_endpoint - m.chi2_overshoot.iloc[0]) / max(1.0, m.chi2_overshoot.iloc[0]))
        print(f"VICE vs cerium_bump, all {len(v)} models: max relative difference in chi2 = {worst:.1e}")

    print(f"\n{'galaxy':7s} {'source':7s} " + " ".join(f"{k:>7s}" for k in KEYS) + "   chi2/N   rank of H21 best")
    for g in H.H21_BEST:
        h21 = {k: (np.nan if H.H21_BEST[g][k] is None else H.H21_BEST[g][k]) for k in KEYS}
        print(f"{g:7s} {'H21':7s} " + " ".join(f"{h21[k]:7.3g}" for k in KEYS))
        for e in ("vice", "flexce"):
            s = df[(df.engine == e) & (df.galaxy == g)].sort_values("chi2").reset_index(drop=True)
            if s.empty:
                continue
            b = s.iloc[0]
            rank = int(np.flatnonzero([same(r, h21) for _, r in s.iterrows()])[0]) + 1
            dchi = s.chi2[rank - 1] - b.chi2
            print(f"{'':7s} {e:7s} " + " ".join(f"{b[k]:7.3g}" for k in KEYS)
                  + f"   {b.chi2 / b.n:6.2f}   {rank:4d}/{len(s)} (Δχ² = {dchi:.0f})")


if __name__ == "__main__":
    main()
