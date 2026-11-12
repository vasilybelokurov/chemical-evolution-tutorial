"""The data extracts in data_build/ match data/MANIFEST.json (SHA-256 and row counts)."""
import hashlib
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = json.loads((ROOT / "data/MANIFEST.json").read_text())


@pytest.mark.parametrize("name", sorted(MANIFEST["files"]))
def test_checksum(name):
    path = ROOT / "data_build" / name
    if not path.exists():
        pytest.skip(f"{name} not in data_build/ (git-ignored; copy or download it)")
    assert hashlib.sha256(path.read_bytes()).hexdigest() == MANIFEST["files"][name]["sha256"]


def test_h21_counts():
    """Table-2 member counts per galaxy; H21 Table 1 gives 972 for GSE but Table 2 lists 1823 (known conflict)."""
    expected = {"lmc": 3909, "smc": 1146, "sgr": 946, "gse": 1823, "fornax": 192}
    for key, n in expected.items():
        assert MANIFEST["files"][f"{key}_h21.ecsv"]["rows"] == n
