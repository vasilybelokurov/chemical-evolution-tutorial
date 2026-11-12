# Chemical evolution of dwarf galaxies

A step-by-step Python tutorial on galactic chemical evolution, in six chapters:

| | Chapter |
|---|---|
| 0 | Five dwarf galaxies in APOGEE (LMC, SMC, Sagittarius, Gaia-Sausage/Enceladus, Fornax) |
| 1 | A one-zone model from scratch, checked against exact solutions |
| 2 | Yields, delay times and the [α/Fe] knee |
| 3 | VICE: the same model with fewer approximations; bursts |
| 4 | Fitting the five galaxies; a reproduction audit of Hasselquist et al. (2021) |
| 5 | What the abundance track cannot see: mass, metallicity distribution, gas |

Each chapter runs on its own in Google Colab (badge at its top) and downloads its data from the GitHub release
`data-v0.1`, checking every file against the SHA-256 in `data/MANIFEST.json`.

## Layout

| path | content |
|---|---|
| `notebooks/` | the chapters `NN_<slug>.ipynb` (committed without outputs) |
| `data/MANIFEST.json` | origin and SHA-256 of every data file |
| `prototypes/` | the development scripts and figures behind Chapters 0 and 3–4 |
| `scripts/colab_vice_check.py` | one-cell check that VICE 1.3.1 builds and runs in a fresh Colab runtime |
| `tests/` | data checksums, prototype reproduction, VICE smoke test |
| `STRUCTURE.md`, `AGENTS.md` | the plan (with its review log) and the working rules |

## Run locally

Python 3.13. VICE 1.3.1 is published only as source and its build fails with Cython ≥ 3.1, hence the first two lines:

```bash
python3.13 -m venv .venv && source .venv/bin/activate
pip install "cython>=3.0,<3.1" setuptools "numpy==2.5.3"
pip install --no-build-isolation vice==1.3.1
pip install -r requirements-dev.txt
pytest -q tests --nbmake notebooks           # runs every chapter (~5 min)
jupyter-book build . -W -n --keep-going
```

## Reference

Hasselquist, S., et al. 2021, ApJ, 923, 172, https://doi.org/10.3847/1538-4357/ac25f9
