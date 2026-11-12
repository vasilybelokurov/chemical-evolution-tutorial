"""The tutorial's VICE port of the H21 models (prototypes/h21_vice) against cerium_bump and against H21's text.

Runs the five H21 best fits once at dt = 0.03 Gyr (about 1 min in total).
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "prototypes/h21_vice"))
import h21_vice as H  # noqa: E402

CB = Path.home() / "Work/Code/cerium_bump"
pytestmark = [pytest.mark.filterwarnings("ignore"),
              pytest.mark.skipif(not H.BUNDLE.exists(), reason="data_build/h21_ssp_bundle.npz missing")]


@pytest.fixture(scope="module")
def tracks():
    return {g: H.run(p, dt=0.03) for g, p in H.H21_BEST.items()}


def feh_at(tr, t):
    return float(np.interp(t, tr["time"], tr["feh"]))


def time_to(tr, feh):
    """First time the model reaches this [Fe/H] (tracks rise monotonically in these models)."""
    ok = np.isfinite(tr["feh"])
    return float(tr["time"][ok][np.argmax(tr["feh"][ok] >= feh)])


def cb_grid():
    pd = pytest.importorskip("pandas")
    path = CB / "results/grids/h21_matched_dt0.03.csv"
    if not path.exists():
        pytest.skip("cerium_bump grid not available")
    return pd.read_csv(path)


def cb_row(grid, g):
    p, m = H.H21_BEST[g], grid["galaxy"] == g
    for k in ("sfe0", "eta", "tau_b", "sigma_b", "f_b"):
        if p[k] is not None:
            m &= np.isclose(grid[k], p[k])
    assert m.sum() == 1
    return grid[m].iloc[0]


@pytest.mark.parametrize("g", list(H.H21_BEST))
def test_same_tracks_as_cerium_bump(tracks, g):
    """Same physics, same VICE: the port must reproduce cerium_bump's saved track to rounding error."""
    ref = np.load(CB / cb_row(cb_grid(), g)["track_path"])
    keep = ref["time"] <= 12 + 1e-8
    for k in ("time", "sfr", "mgas", "feh", "sife", "mgfe"):
        np.testing.assert_allclose(tracks[g][k], ref[k][keep], rtol=1e-12, atol=1e-12, equal_nan=True)


@pytest.mark.parametrize("g", list(H.H21_BEST))
def test_same_chi2_as_cerium_bump(tracks, g):
    from astropy.table import Table
    t = Table.read(ROOT / f"data_build/{g.lower()}_h21.ecsv")
    chi, n = H.chi2_track(np.c_[t["feh"], t["sife"]].astype(float), np.c_[t["feh_err"], t["sife_err"]].astype(float),
                          np.c_[tracks[g]["feh"], tracks[g]["sife"]])
    row = cb_row(cb_grid(), g)
    assert n == row["n_used"]
    assert chi == pytest.approx(row["chi2_track"], rel=1e-9)


# H21 §5.3 statements about their flexCE best fits. They quote round numbers, so allow 0.1 dex and 0.5 Gyr.
def test_h21_feh_at_5_gyr(tracks):
    """'by 5 Gyr into its evolution, GSE achieved [Fe/H] = -0.5, as compared to -1.2 and -1.5 for the LMC and SMC'"""
    for g, expected in [("GSE", -0.5), ("LMC", -1.2), ("SMC", -1.5)]:
        assert feh_at(tracks[g], 5.0) == pytest.approx(expected, abs=0.1), g


def test_h21_lmc_late_enrichment(tracks):
    """'it took the LMC 10 Gyr to enrich to [Fe/H] = -0.7'; 'enriching from -0.8 to -0.3 in the last 2-3 Gyr'"""
    lmc = tracks["LMC"]
    assert time_to(lmc, -0.7) == pytest.approx(10.0, abs=0.5)
    assert feh_at(lmc, 9.5) == pytest.approx(-0.8, abs=0.1)
    assert feh_at(lmc, 12.0) == pytest.approx(-0.3, abs=0.1)


def test_h21_smc_burst_4_gyr_before_lmc(tracks):
    """'the burst in the SMC ... occurred ~4 Gyr earlier than the burst in the LMC' (SFR peak times)"""
    peak = {g: tracks[g]["time"][np.argmax(tracks[g]["sfr"])] for g in ("LMC", "SMC")}
    assert peak["LMC"] - peak["SMC"] == pytest.approx(4.0, abs=0.5)


def test_h21_sgr_second_epoch(tracks):
    """'a second star formation epoch, beginning some 5-6 Gyr into its evolution'"""
    t, sfr = tracks["Sgr"]["time"], tracks["Sgr"]["sfr"]
    late = t > 3
    t_min = t[late][np.argmin(sfr[late])]
    after = t > t_min
    t_peak2 = t[after][np.argmax(sfr[after])]
    assert 4.5 <= t_min <= 6.0 and 5.0 <= t_peak2 <= 7.0


@pytest.mark.xfail(reason="known: VICE reaches -0.7 about 1 Gyr earlier than H21's flexCE GSE model", strict=True)
def test_h21_gse_time_to_minus_0p7(tracks):
    """'GSE took ~6 Gyr [Lian] to reach [Fe/H] = -0.7 instead of ~4 Gyr in the flexCE model'"""
    assert time_to(tracks["GSE"], -0.7) == pytest.approx(4.0, abs=0.5)
