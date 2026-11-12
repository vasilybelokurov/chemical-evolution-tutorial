# STRUCTURE.md — revised after review (2026-11-12)

Plan of the tutorial **Chemical evolution of dwarf galaxies**. Reviewed by Fable and Codex (§ Review log); each point is
adopted or rejected there with a reason.

**Status (2026-11-12):** chapters 00–05 written as planned below, all asserts passing.

- **Audience:** students who know numpy and matplotlib, and what [X/Fe] means.
- **Exit level:**
  - write a one-zone GCE model and check it against analytic solutions;
  - explain the [α/Fe] knee with yields and SN Ia delay times;
  - run VICE and say which approximations it relaxes;
  - fit the [Si/Fe]–[Fe/H] tracks of real dwarf galaxies;
  - show with a mock that the fit recovers known input, and say what the fit cannot constrain.

## Principles

1. **Toy first, then data.** Every effect is first produced by a model with known input, then looked for in APOGEE.
   Every fit is first run on mock stars drawn from a known model.
2. **Same equations, approximations relaxed one at a time.** The numpy model of Chapter 01 uses instantaneous
   recycling and constant yields. Chapter 03 first runs VICE with the same choices (`recycling = r`, constant yields,
   exponential SN Ia DTD) and asserts agreement. It then switches on one realism item per plot (continuous recycling,
   metallicity-dependent yields, AGB, lifetimes) and reports the shift in dex.
3. **Doubt first.** These are not clean facts:
   - sample membership: GSE has 1823 rows in H21 Table 2 against 972 in Table 1, and 11 IDs are in both GSE and Sgr;
   - the sparse metal-poor tail of Fornax: 46 stars below [Fe/H] = −1.1, with foreground contamination unchecked;
   - [Fe/H] as a common clock: it is not;
   - H21's best fits: two of the five do not reproduce.
4. **Plot + assertion**, with the kind of each assertion named, as in gaia_nss_tutorial:
   - *identity*, exact by construction;
   - *limit*, against an analytic solution;
   - *recovery*, a known input recovered from mocks;
   - *reproduction*, against a published number;
   - *regression*, a frozen value that detects change and proves nothing physical.
5. **One quality mask** for plotting and fitting. It is defined once in Chapter 00 and imported everywhere.
6. **Runs without VICE where possible.** Chapters 04–05 work from the precomputed grids, so building VICE (3–5 min) is
   needed only in Chapter 03 and in the optional cells.

## Chapters

| # | Chapter | Content | Assert (kind) |
|---|---|---|---|
| 00 | The data | APOGEE DR17 members of five dwarfs (H21 Table 2); [Mg/Fe] and [Si/Fe] tracks, stating that the fits use Si as H21 did; the quality mask; running medians with bootstrap bands; membership doubts: GSE count, GSE/Sgr overlaps, Fornax tail | H21 Table 1 counts (reproduction; GSE a documented failure); medians (regression); conclusions unchanged without overlaps and the tail (sensitivity) |
| 01 | A one-zone model from scratch | gas, SFR = M_gas/τ★, inflow, outflow η, instantaneous recycling, the equilibrium abundance | mass conservation (identity); constant-SFR and exponential-inflow solutions of Weinberg, Andrews & Freudenburg 2017 (limit); timestep convergence |
| 02 | Yields and delays | CCSN vs SN Ia, DTD normalisation, plateau and knee. Constant yields, plus `h21_ssp_bundle.npz` read with numpy; no VICE | plateau = CCSN yield ratio (identity); [Fe/H] of the knee shifts by −Δlog τ★ at fixed DTD (controlled experiment, tolerance from the analytic solution) |
| 03 | VICE | matched-physics VICE equals Chapter 01, then one approximation relaxed per plot; time-dependent SFE (bursts); [Fe/H](t) of five models (the clock that is not common) | numpy vs matched VICE within a stated tolerance; shifts reported (regression) |
| 04 | Fitting five dwarfs | the H21 model and grid from the precomputed grids. First the score: distance to the track, endpoint penalty (reconstructed, with its last-Gyr exemption), χ²/N = 4.7–8.5, scaled Δχ². Then mock recovery, profiles and degeneracies, and an audit of the reproduction. Box: "is it the code?", VICE vs flexCE from `h21_grid_scores.csv` | mock recovery over an ensemble (recovery); re-scored grid = shipped scores to 1e-4 (regression); H21 model ranks ≤ 3 for Sgr, GSE and LMC (reproduction); SMC and Fornax ranks frozen and taught as a failed reproduction |
| 05 | What the track cannot see | mass scaling: the track is unchanged when M_i and M₀ are rescaled together. Where along the track the stars sit: two models with near-identical tracks but different mass per [Fe/H] bin. The gas left today. A mock "observed" mass distribution drawn from a known model, and a joint fit that recovers what abundances alone cannot. No real LMC mass data (author's decision, 2026-11-12) | mass-scaling identity (identity); mass per bin sums to ∫SFR dt (identity); joint-fit recovery on the mock (recovery) |
| later | Beyond Mg (element ranking, Ce); multizone (`vice.milkyway`) | needs more extracts and a Ce yield model | — |

## Open questions

1. **Real stellar mass per [Fe/H]: DECIDED 2026-11-12, left out.** No public selection-corrected measurement is used;
   Chapter 05 uses mocks only.
2. **Colab: PASSED 2026-11-12.** On real Colab the pinned recipe builds VICE in 459 s (7.7 min) and the smoke test
   matches macOS and CI exactly. That build time is the strongest argument for Principle 6, and for
   a prebuilt VICE wheel on the data release with the source build as fallback. Colab notebooks must not install
   the pinned `numpy==2.5.3` from requirements.txt.
3. **Python:** Colab now runs **3.13.18**, not 3.12. Target 3.13, as gaia_nss_tutorial does; VICE 1.3.1 builds and
   passes on 3.13 (macOS, CI, Colab). The local GCE venv stays on 3.12 until it is rebuilt.
4. **CI without the author's disk.** `tests/test_h21_vice.py` reads `~/Work/Code/cerium_bump`, and the checksum tests
   skip missing files. Once the data release exists, CI downloads it and fails on missing assets. The cerium_bump
   comparisons stay local-only, marked as such.

## Review log (2026-11-12)

Fable (`claude-fable-5-2`, Agent) and Codex (ChatGPT login, read-only, effort high) each reviewed the draft
independently. Each claim was checked against the code before being adopted.

| Point | Raised by | Verdict | Action |
|---|---|---|---|
| "One model throughout" is false: Ch 01 uses instantaneous recycling, the VICE model continuous recycling, Z-dependent yields and AGB (`h21_vice.py:115–125`) | both | confirmed; VICE accepts a number for `recycling` (checked in `help(vice.singlezone.recycling)`) | Principle 2 |
| Ch 04 asserts against H21 best fits that do not reproduce for SMC and Fornax | both | confirmed (Chapter 4) | Ch 04 asserts reworked |
| No parameter-recovery test | both | agreed | Principle 1; Ch 04 mocks |
| Mg vs Si mismatch: Ch 00 shows Mg, the fit uses Si (`score_grid.py:27`) | both | confirmed | Ch 00 shows both |
| Plotting and fitting use different masks (flags in `plot_dwarf_mgfe.py:31`; finite and error > 0 only in `h21_vice.py:152`) | Codex | confirmed | Principle 5 |
| χ² is a distance score, not a likelihood (χ²/N 4.7–8.5) | both | agreed | Ch 04 opens with the score |
| Mass/MDF findings need a home; the LMC data are not public | both | agreed; author decided to leave the real data out | Ch 05 with mocks only |
| "Beyond Mg" is circular and lacks a Ce model | both | agreed | moved to later |
| Ch 02 uses VICE and flexCE tables before VICE is introduced | Fable | agreed | Ch 02: constant yields + bundle |
| Weak asserts: Table 2 counts are circular; knee scaling unspecified | both | agreed | Table 1; −Δlog τ★ experiment |
| Bursts appear unprepared in Ch 04 | Fable | agreed | Ch 03 adds time-dependent SFE |
| VICE vs flexCE: Fable puts it in Ch 04, Codex in Ch 03 | both | Ch 04: it answers "why don't SMC/Fornax reproduce?" | box in Ch 04 |
| Code agreement does not prove the inputs cause the H21 mismatch | Codex | agreed: it rules out the code, not other causes | wording in Ch 04 |
| Membership sensitivity: overlaps, Fornax tail | both | agreed | Ch 00 sensitivity check |
| Colab "yes" overstates CI-only evidence; numpy pin; wheel | both | agreed | open question 2 |
| CI must fail on missing assets; cerium_bump tests skip | Codex, Fable | confirmed (`tests/test_data.py:15`) | open question 4; after the data release |
| Keep the population chapter short for this audience | Codex (its own counterargument) | agreed | Ch 05: one worked example + exercises |
| Python 3.12 | Fable | rejected: real Colab runs 3.13.18 | target 3.13 |
