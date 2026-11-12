"""The dwarf [Mg/Fe] prototype reproduces the numbers it was developed with (2026-11-11)."""
import sys
from pathlib import Path

import numpy as np
import pytest

PROTO = Path(__file__).resolve().parents[1] / "prototypes/dwarf_mgfe"
sys.path.insert(0, str(PROTO))
pytest.importorskip("astropy")
if not (PROTO.parents[1] / "data_build/fornax_h21.ecsv").exists():
    pytest.skip("data_build/ extracts missing", allow_module_level=True)

from plot_dwarf_mgfe import load, running_median  # noqa: E402
from plot_dwarf_mgfe_tracks import TAIL_EDGE, equal_count_median  # noqa: E402

# Stars after the quality cut, track [Fe/H] range, median [Mg/Fe] first -> last (from the development run).
FROZEN = {"lmc": (3900, -1.61, -0.41, 0.091, -0.005), "smc": (1138, -1.57, -0.67, 0.030, -0.096),
          "sgr": (935, -1.50, -0.20, 0.103, -0.126), "gse": (1811, -1.99, -0.64, 0.298, 0.003),
          "fornax": (183, -1.07, -0.77, -0.134, -0.185)}


@pytest.mark.parametrize("key", sorted(FROZEN))
def test_running_median(key):
    n, lo, hi, m_max, m_min = FROZEN[key]
    x, y = load(key)
    c, m, *_ = running_median(x, y)
    assert len(x) == n
    np.testing.assert_allclose([c.min(), c.max(), m.max(), m.min()], [lo, hi, m_max, m_min], atol=0.006)


def test_fornax_tail():
    x, y = load("fornax")
    tail = x < TAIL_EDGE
    tc, *_ = equal_count_median(x[tail], y[tail])
    assert tail.sum() == 46 and len(tc) == 19
    assert tc.min() == pytest.approx(-1.96, abs=0.006)
