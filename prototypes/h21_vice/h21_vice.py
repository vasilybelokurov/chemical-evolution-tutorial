"""Hasselquist et al. (2021, H21) one-zone models of five dwarf galaxies, run with VICE.

Ported and trimmed from cerium_bump (git 7c52f48) src/dwarf_vice/{model,yields,vice_runner,fitting}.py.
H21 ran these models with flexCE; here VICE runs them with the same (public flexCE) yields, built by
cerium_bump into data_build/h21_ssp_bundle.npz (SHA-256 in data/MANIFEST.json).

Model (H21 §5.1, Tables 3-4; time t in Gyr, masses in Msun):
    inflow   Mdot_in(t) = (M_i / tau_i) (t / tau_i) exp(-t / tau_i)            "delayed tau"
    SFE      SFE(t)     = SFE0 [1 + (f_b - 1) exp(-(t - tau_b)^2 / (2 sigma_b^2))]   Gaussian burst (f_b > 1) or dip (f_b < 1)
    SFR      = SFE(t) M_gas,   outflow = eta * SFR,   initial gas M_0 = 3e9, inflow mass M_i = 6e10.

Usage (from the repo root, GCE venv):
    python prototypes/h21_vice/h21_vice.py            # runs the five H21 best fits, prints a summary
"""
import os
import tempfile
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BUNDLE = ROOT / "data_build" / "h21_ssp_bundle.npz"
ELEMENTS = ("he", "c", "n", "o", "ne", "mg", "si", "s", "ca", "fe", "ni")
SOLAR_LOGEPS = {"fe": 7.47, "si": 7.54, "mg": 7.55}  # Lodders (2003) photosphere, as in flexCE
ATOMIC_MASS = {"fe": 56, "si": 28, "mg": 24}          # dominant isotopes

# H21 Table 3 (fixed) and Table 4 (best fit). f_b = 1 means no burst; sigma_b fixed for the MCs.
H21_BEST = {
    "LMC":    dict(sfe0=0.015, eta=2.0,   tau_i=2.0, tau_b=11.0, sigma_b=0.75, f_b=6.0),
    "SMC":    dict(sfe0=0.008, eta=10.0,  tau_i=2.0, tau_b=7.0,  sigma_b=0.75, f_b=4.0),
    "Sgr":    dict(sfe0=0.03,  eta=17.5,  tau_i=1.0, tau_b=5.0,  sigma_b=0.5,  f_b=0.01),
    "GSE":    dict(sfe0=0.14,  eta=6.0,   tau_i=2.5, tau_b=None, sigma_b=None, f_b=1.0),
    "Fornax": dict(sfe0=0.03,  eta=100.0, tau_i=1.0, tau_b=6.0,  sigma_b=5.0,  f_b=0.01),
}
M_GAS_0, M_INFLOW = 3e9, 6e10

# H21 Table 4 model grids: (first, last, step) for each varied parameter. Unvaried ones come from H21_BEST
# (Table 3 values). Strings avoid floating-point drift in the steps (0.006 + 4 * 0.001 must be exactly 0.010).
H21_GRID = {
    "LMC":    dict(sfe0=("0.010", "0.030", "0.005"), eta=("2", "10", "2"), tau_b=("8", "12", "1"), f_b=("4", "8", "1")),
    "SMC":    dict(sfe0=("0.006", "0.010", "0.001"), eta=("5", "25", "5"), tau_b=("6", "11", "1"), f_b=("2", "10", "2")),
    "Sgr":    dict(sfe0=("0.02", "0.06", "0.01"), eta=("12.5", "22.5", "2.5"), tau_b=("4", "8", "1"),
                   sigma_b=("0.25", "1.25", "0.25")),
    "GSE":    dict(sfe0=("0.08", "0.20", "0.01"), eta=("1", "11", "1")),
    "Fornax": dict(sfe0=("0.01", "0.05", "0.01"), eta=("20", "120", "20"), tau_b=("2", "8", "1"), sigma_b=("1", "6", "1")),
}


def grid(galaxy):
    """All parameter sets of the H21 grid for one galaxy (list of dicts, same keys as H21_BEST)."""
    from decimal import Decimal
    from itertools import product
    axes = {}
    for k, (lo, hi, step) in H21_GRID[galaxy].items():
        lo, hi, step = Decimal(lo), Decimal(hi), Decimal(step)
        axes[k] = [float(lo + i * step) for i in range(int((hi - lo) / step) + 1)]
    return [{**H21_BEST[galaxy], **dict(zip(axes, values))} for values in product(*axes.values())]


def inflow_rate(t, tau_i, mass=M_INFLOW):
    """Delayed-tau inflow in Msun/yr (VICE wants per yr; t and tau_i in Gyr)."""
    return mass / tau_i * (t / tau_i) * np.exp(-t / tau_i) / 1e9


def sfe(t, p):
    """Star-formation efficiency in Gyr^-1 with H21's Gaussian modulation."""
    if p["tau_b"] is None:
        return p["sfe0"]
    return p["sfe0"] * (1 + (p["f_b"] - 1) * np.exp(-0.5 * ((t - p["tau_b"]) / p["sigma_b"]) ** 2))


# ---- yields: metallicity-dependent net yields from the bundle, looked up on flexCE's Z grid -----------
@lru_cache(maxsize=1)
def bundle():
    with np.load(BUNDLE) as f:
        return {k: f[k] for k in f.files}


def z_index(Z):
    z = bundle()["z"]
    return min(len(z) - 1, int(np.searchsorted(z, Z, side="right")))   # flexCE's upper-bin lookup


class CCYield:            # CCSN and SN Ia: yield per unit mass formed, a function of Z only
    def __init__(self, key): self.key = key
    def __call__(self, Z): return float(bundle()[self.key][z_index(Z)])


class AGBYield:           # AGB: linear interpolation in progenitor mass m, divided by m (VICE convention)
    def __init__(self, key): self.key = key
    def __call__(self, m, Z):
        mass = bundle()["agb_mass"]
        if m < mass[0] or m > 8:
            return 0.0
        row = bundle()[self.key][z_index(Z)]
        j = min(len(mass) - 2, max(0, int(np.searchsorted(mass, m)) - 1))
        f = min(1.0, max(0.0, (m - mass[j]) / (mass[j + 1] - mass[j])))
        return float((row[j] * (1 - f) + row[j + 1] * f) / m)


def configure_yields():
    import vice
    for el in ELEMENTS:
        vice.yields.ccsne.settings[el] = CCYield(f"cc_{el}")
        vice.yields.sneia.settings[el] = CCYield(f"ia_{el}")
        vice.yields.agb.settings[el] = AGBYield(f"agb_{el}")


# ---- one VICE run ---------------------------------------------------------------------------------------
def run(p, dt=0.01, duration=12.0):
    """Run one H21 model in VICE. Returns arrays: time [Gyr], sfr [Msun/yr], mgas [Msun], feh, sife, mgfe."""
    import vice
    configure_yields()
    vice.mlr.setting = "pm1993"
    here = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        os.chdir(tmp)       # VICE splits output paths at spaces: run in a space-free temporary directory
        try:
            sz = vice.singlezone(name="h21", elements=ELEMENTS, mode="ifr",
                                 func=lambda t: inflow_rate(t, p["tau_i"]),
                                 tau_star=lambda t: 1 / sfe(t, p),    # VICE's depletion time = 1/SFE
                                 Mg0=M_GAS_0, eta=p["eta"], dt=dt, IMF="kroupa", m_lower=0.1, m_upper=100.0,
                                 RIa="exp", delay=0.15, tau_ia=1.5, recycling="continuous", postMS=0.0,
                                 Z_solar=0.014, schmidt=False, smoothing=0.0, Zin=0.0, enhancement=1.0)
            for channel in ("agb", "ccsne", "sneia"):
                for el in ELEMENTS:
                    getattr(sz.entrainment, channel)[el] = 1.0   # all ejecta enter the ISM
            h = sz.run(np.arange(round(duration / dt) + 1) * dt, capture=True, overwrite=True).history
            out = {k: np.asarray(h[c]) for k, c in [("time", "time"), ("sfr", "sfr"), ("mgas", "mgas")]}
            z = {el: np.asarray(h[f"z({el})"]) for el in ELEMENTS}
        finally:
            os.chdir(here)
    # Number ratios on the Lodders (2003) scale; H = 1 - He - tracked metals (mass fractions).
    hydrogen = 1 - z["he"] - sum(z[el] for el in ELEMENTS if el != "he")
    with np.errstate(divide="ignore", invalid="ignore"):
        logn = {el: np.log10(z[el] / ATOMIC_MASS[el]) for el in ATOMIC_MASS}
        out["feh"] = logn["fe"] - np.log10(hydrogen) - (SOLAR_LOGEPS["fe"] - 12)
        out["sife"] = logn["si"] - logn["fe"] - (SOLAR_LOGEPS["si"] - SOLAR_LOGEPS["fe"])
        out["mgfe"] = logn["mg"] - logn["fe"] - (SOLAR_LOGEPS["mg"] - SOLAR_LOGEPS["fe"])
    keep = out["time"] <= duration + 1e-6 * dt
    return {k: v[keep] for k, v in out.items()}


# ---- H21 chi^2: each star's distance to the track, divided by its error in that direction ----------------
def chi2_track(points, errors, track_xy):
    """points, errors: (N, 2) arrays of ([Fe/H], [Si/Fe]) and their errors; track_xy: (M, 2) polyline.

    For each star find the nearest point q on the piecewise-linear track. With residual d = p - q,
    the error along d is sigma_d^2 = sum(d_k^2 sigma_k^2) / |d|^2, so chi^2 = |d|^4 / sum(d_k^2 sigma_k^2).
    """
    ok = np.all(np.isfinite(points) & (points > -9000) & np.isfinite(errors) & (errors > 0), axis=1)
    p, e = points[ok], errors[ok]
    v = np.all(np.isfinite(track_xy), axis=1)
    seg = v[:-1] & v[1:]
    a, b = track_xy[:-1][seg], track_xy[1:][seg]
    ab = b - a
    ab2 = np.einsum("ij,ij->i", ab, ab)
    chi = np.empty(len(p))
    for s in range(0, len(p), 128):                       # chunks keep the (stars x segments) arrays small
        d = p[s:s + 128, None, :] - a                      # star minus segment start
        f = np.clip(np.divide(np.einsum("nsi,si->ns", d, ab), ab2, out=np.zeros(d.shape[:2]), where=ab2 > 0), 0, 1)
        q = a + f[:, :, None] * ab                         # nearest point on each segment
        d2 = np.sum((p[s:s + 128, None, :] - q) ** 2, axis=2)
        k = np.argmin(d2, axis=1)
        r = p[s:s + 128] - q[np.arange(len(k)), k]
        den = np.sum(r ** 2 * e[s:s + 128] ** 2, axis=1)
        chi[s:s + 128] = np.divide(d2[np.arange(len(k)), k] ** 2, den, out=np.zeros(len(k)), where=den > 0)
    return float(chi.sum()), int(ok.sum())


def chi2_endpoint(points, errors, track_xy, time, galaxy, dt, duration=12.0):
    """H21's extra term for models that evolve past the most metal-rich stars (H21 §5.1).

    H21 give the idea, not the formula; this is cerium_bump's reconstruction ('endpoint30_reference30myr_v1'):
    the end point is (max observed [Fe/H], median [Si/Fe] of the 30 most metal-rich stars), with the median
    errors of those 30 stars. Every model step beyond that [Fe/H] and earlier than the last 1 Gyr adds
    |d|^4 / sum(d_k^2 sigma_k^2), as in chi2_track, weighted by dt / 0.03 Gyr. Not applied to GSE.
    """
    if galaxy == "GSE":
        return 0.0
    ok = np.all(np.isfinite(points) & (points > -9000) & np.isfinite(errors) & (errors > 0), axis=1)
    p, e = points[ok], errors[ok]
    top = np.argsort(p[:, 0], kind="stable")[-30:]
    end, err = np.array([p[:, 0].max(), np.median(p[top, 1])]), np.median(e[top], axis=0)
    use = np.all(np.isfinite(track_xy), axis=1) & (track_xy[:, 0] > end[0]) & (time <= duration - 1 + 1e-9)
    d = track_xy[use] - end
    d2, den = np.sum(d ** 2, axis=1), np.sum(d ** 2 * err ** 2, axis=1)
    return float(np.sum(np.divide(d2 * d2, den, out=np.zeros_like(d2), where=den > 0)) * dt / 0.03)


if __name__ == "__main__":
    import time
    for g, p in H21_BEST.items():
        t0 = time.time()
        tr = run(p)
        i5 = np.argmin(np.abs(tr["time"] - 5))
        print(f"{g:6s} {time.time() - t0:5.1f} s  [Fe/H]_end={tr['feh'][-1]:+.3f}  [Fe/H](5 Gyr)={tr['feh'][i5]:+.3f}"
              f"  SFR peak at {tr['time'][np.argmax(tr['sfr'])]:.2f} Gyr")
