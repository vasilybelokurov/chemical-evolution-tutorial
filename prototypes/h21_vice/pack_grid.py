"""Pack the full H21 grid runs into compact tutorial data files (float32 tracks) and check nothing was lost.

Reads  prototypes/h21_vice/grid_runs/h21_grid_<engine>.npz  (from run_grid.py, float64)
       prototypes/h21_vice/grid_runs/h21_grid_scores.csv     (from score_grid.py)
Writes data_build/h21_grid_<engine>.npz  (same keys; tracks as float32; params float64; time float64)
       data_build/h21_grid_scores.csv    (chi2 per model and engine)
Check: re-scoring the float32 tracks reproduces every chi2 to < 1e-4 relative.

Usage (repo root, GCE venv): python prototypes/h21_vice/pack_grid.py
"""
import shutil

import numpy as np
import pandas as pd

import h21_vice as H
from score_grid import RUNS, stars

FIELDS = ("sfr", "mgas", "feh", "sife", "mgfe")
OUT = H.ROOT / "data_build"


def main():
    ref = pd.read_csv(RUNS / "h21_grid_scores.csv")
    data = {g: stars(g) for g in H.H21_BEST}
    for engine in ("vice", "flexce"):
        z = dict(np.load(RUNS / f"h21_grid_{engine}.npz"))
        for f in FIELDS:
            z[f] = z[f].astype(np.float32)
        path = OUT / f"h21_grid_{engine}.npz"
        np.savez_compressed(path, **z)
        # Re-score from the packed file and compare with the float64 scores.
        p = np.load(path)
        r = ref[ref.engine == engine].reset_index(drop=True)
        worst = 0.0
        for i, g in enumerate(p["galaxy"]):
            pts, err = data[g]
            track = np.c_[p["feh"][i], p["sife"][i]].astype(float)
            chi, _ = H.chi2_track(pts, err, track)
            chi += H.chi2_endpoint(pts, err, track, p["time"], g, dt=0.03)
            worst = max(worst, abs(chi / r.chi2[i] - 1))
        print(f"{engine}: wrote {path.name} ({path.stat().st_size / 1e6:.1f} MB); "
              f"max relative chi2 change after float32 = {worst:.1e}")
        assert worst < 1e-4
    shutil.copy(RUNS / "h21_grid_scores.csv", OUT / "h21_grid_scores.csv")
    print("copied h21_grid_scores.csv")


if __name__ == "__main__":
    main()
