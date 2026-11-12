"""Run the full H21 grid (3403 models, H21 Table 4) with VICE or with flexCE, in parallel.

Output: prototypes/h21_vice/grid_runs/h21_grid_<engine>.npz with
    galaxy (n,), sfe0, eta, tau_i, tau_b, sigma_b, f_b (n,; NaN where not used),
    time (m,) [Gyr], and (n, m) arrays sfr [Msun/yr], mgas [Msun], feh, sife, mgfe.

Usage (repo root):
    ~/Work/venvs/gce/bin/python    prototypes/h21_vice/run_grid.py vice   [--workers 12]   # ~10 min
    ~/Work/venvs/flexce/bin/python prototypes/h21_vice/run_grid.py flexce [--workers 12]   # ~40 min
VICE uses dt = 0.03 Gyr (as cerium_bump's grid); flexCE uses its fiducial 30 Myr.
"""
import argparse
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import h21_vice as H

KEYS = ("sfe0", "eta", "tau_i", "tau_b", "sigma_b", "f_b")
FIELDS = ("sfr", "mgas", "feh", "sife", "mgfe")


def work(job):
    engine, p = job
    warnings.filterwarnings("ignore")
    if engine == "vice":
        tr = H.run(p, dt=0.03)
    else:
        from gce_tools.flexce_h21 import run
        tr = run(p)
    return tr["time"], np.array([tr[k] for k in FIELDS])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("engine", choices=("vice", "flexce"))
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()
    jobs = [(g, p) for g in H.H21_BEST for p in H.grid(g)]
    t0 = time.time()
    out = []
    with Pool(args.workers) as pool:
        for i, res in enumerate(pool.imap(work, [(args.engine, p) for _, p in jobs], chunksize=4)):
            out.append(res)
            if (i + 1) % 250 == 0 or i + 1 == len(jobs):
                print(f"{i + 1}/{len(jobs)} models, {time.time() - t0:.0f} s", flush=True)
    times = out[0][0]
    assert all(np.array_equal(t, times) for t, _ in out), "models returned different time grids"
    data = np.stack([a for _, a in out])                         # (n_models, n_fields, n_times)
    params = {k: np.array([np.nan if p[k] is None else p[k] for _, p in jobs], dtype=float) for k in KEYS}
    path = Path(__file__).with_name("grid_runs") / f"h21_grid_{args.engine}.npz"
    path.parent.mkdir(exist_ok=True)
    np.savez_compressed(path, galaxy=np.array([g for g, _ in jobs]), time=times, **params,
                        **{f: data[:, i] for i, f in enumerate(FIELDS)})
    print(f"wrote {path} ({path.stat().st_size / 1e6:.1f} MB) in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
