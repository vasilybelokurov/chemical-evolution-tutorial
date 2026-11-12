"""VICE from the shared GCE toolchain runs a one-zone model with the expected qualitative behaviour."""
import pytest

vice = pytest.importorskip("vice")


def test_vice_singlezone(tmp_path):
    assert vice.__version__ == "1.3.1"
    sz = vice.singlezone(name=str(tmp_path / "sz"), elements=("fe", "mg"), dt=0.01)
    sz.run([i * 0.1 for i in range(101)], overwrite=True)
    h = vice.history(str(tmp_path / "sz"))
    feh, mgfe = h["[fe/h]"], h["[mg/fe]"]
    assert feh[-1] > feh[10]          # iron builds up
    assert mgfe[-1] < mgfe[10]        # SN Ia lower [Mg/Fe] with time
