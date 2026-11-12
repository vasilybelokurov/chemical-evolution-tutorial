# AGENTS.md

## Project

A Python tutorial, written as Jupyter notebooks and published from a public git repo, on galactic chemical evolution
(GCE). It builds one-zone models from first principles, then runs them with VICE, and compares them with the
[α/Fe]–[Fe/H] tracks of five dwarf galaxies (LMC, SMC, Sgr, GSE, Fornax) observed by APOGEE DR17. Students get the
data as pre-extracted files. **Public data only.**

## Rules

1. **Write the least code that works.**
2. **Always reproduce bugs with a failing test** before fixing them.
3. **Show a plot before claiming a result.**

## Workflow

- Chapters are `notebooks/NN_<slug>.ipynb`. Commit them without outputs. Each one must run on its own in Colab;
  VICE 1.3.1 installs there only as `pip install "cython>=3.0,<3.1" setuptools && pip install
  --no-build-isolation vice==1.3.1` (a plain install fails with Cython ≥ 3.1).
- Students have no access to the APOGEE allStar file. `scripts/build_data.py` writes the small extracts that the
  notebooks download, and each download is checked against the SHA-256 in `data/MANIFEST.json`. Published files are
  never changed in place: a changed selection gets a new file name. The allStar catalogue itself is never copied
  into this repository.
- Every claim in a notebook needs three things:
  - a numbered figure, with a **Figure N.k** caption stating what is plotted, the axes, units, colours and samples;
  - a discussion of what the figure shows;
  - an `assert` with a stated tolerance against a known result: a published value, an analytic solution, or a run
    of the model with known input.
- Comment the code for beginners. Where a line implements a formula, say which one; explain non-obvious numpy idioms.
- GCE toolchain: VICE 1.3.1 and the pinned flexCE yields come from `~/Work/Code/chemical_evolution` (venv
  `~/Work/venvs/gce`). Do not rebuild them here.
- Before pushing, run `pytest tests/ --nbmake notebooks/` and `jupyter-book build . -W -n --keep-going`.
- Every literature reference needs a working ADS, arXiv or DOI link, checked before it is committed.
