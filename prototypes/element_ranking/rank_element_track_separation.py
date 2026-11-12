#!/usr/bin/env python3
"""Rank APOGEE DR17 elements by how different the five H21 dwarf tracks are in [X/Fe] vs [Fe/H].

Samples: results/samples/{lmc,smc,sgr,gse,fornax}_h21.ecsv (H21 Table-2 members; catalog_row indexes the allStar file).
Element quality: X_FE finite, > -9000 and X_FE_FLAG == 0 (element-specific; the samples already carry the stellar cuts).
Metric, per element, in [Fe/H] bins of width BIN between FEH_LO and FEH_HI with >= NMIN stars per galaxy:
  spread_dex   = max - min of the galaxy bin medians (dex)
  sep_scatter  = spread_dex / pooled within-galaxy robust scatter (1.4826 MAD) in that bin -> "how many scatters apart"
  chi2_nu      = chi^2/nu of the bin medians about their weighted mean, using median errors 1.2533*sigma/sqrt(N)
Element summary = median over bins with >= 2 galaxies; ranking by sep_scatter (visual track separation).
Outputs: results/element_track_separation/{all,without_<EXCLUDE>}/{ranking.csv,bins.csv}, plots/element_track_separation_top.png
Run (from project root): .venv/bin/python scripts/rank_element_track_separation.py
"""
from pathlib import Path
import os
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
ROOT = Path.cwd()
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'data/matplotlib_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from dwarf_vice.observations import config, GALAXIES, valid_abundance
from plot_h21_samples import COLORS

BIN, FEH_LO, FEH_HI, NMIN = 0.1, -2.5, 0.0, 15
EXCLUDE = [g for g in os.environ.get('EXCLUDE', '').split(',') if g]   # robustness: drop galaxies, e.g. EXCLUDE=GSE
OUT = ROOT / 'results/element_track_separation' / ('without_' + '_'.join(EXCLUDE) if EXCLUDE else 'all')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cfg = config(); path = Path(os.environ.get('APOGEE_ALLSTAR', cfg['catalog']))
    with fits.open(path, memmap=True) as h:
        d = h[1].data
        elems = [c[:-3] for c in d.columns.names if c.endswith('_FE') and c + '_FLAG' in d.columns.names and c != 'M_H']
        stars = []
        for g in [g for g in GALAXIES if g not in EXCLUDE]:
            s = Table.read(ROOT / f'results/samples/{g.lower()}_h21.ecsv')
            rows = np.asarray(s['catalog_row'])
            assert np.all(np.char.strip(np.asarray(d['APOGEE_ID'][rows], str)) == np.asarray(s['apogee_id'], str))
            t = pd.DataFrame(dict(galaxy=g, feh=np.asarray(s['feh'], float)))
            for e in elems:
                v = np.asarray(d[e + '_FE'][rows], float); f = np.asarray(d[e + '_FE_FLAG'][rows])
                t[e] = np.where(valid_abundance(v) & (f == 0), v, np.nan)
            stars.append(t)
    S = pd.concat(stars, ignore_index=True)
    S = S[valid_abundance(S.feh) & (S.feh >= FEH_LO) & (S.feh < FEH_HI)]
    S['bin'] = np.floor((S.feh - FEH_LO) / BIN).astype(int)
    rows = []
    for e in elems:
        for b, grp in S.dropna(subset=[e]).groupby('bin'):
            stats = []
            for g, x in grp.groupby('galaxy'):
                if len(x) >= NMIN:
                    v = x[e].to_numpy(); sig = 1.4826 * np.median(np.abs(v - np.median(v)))
                    stats.append((g, np.median(v), sig, 1.2533 * sig / np.sqrt(len(v)), len(v)))
            if len(stats) < 2:
                continue
            g_, med, sig, err, n = map(np.array, zip(*stats))
            w = 1 / np.maximum(err, 1e-3) ** 2; mean = np.sum(w * med) / np.sum(w)
            pooled = np.sqrt(np.sum((n - 1) * sig ** 2) / np.sum(n - 1))
            rows.append(dict(element=e, feh_mid=FEH_LO + (b + .5) * BIN, n_gal=len(stats), galaxies=','.join(g_),
                             spread_dex=med.max() - med.min(), pooled_scatter=pooled,
                             sep_scatter=(med.max() - med.min()) / max(pooled, 1e-3),
                             chi2_nu=np.sum(w * (med - mean) ** 2) / (len(stats) - 1)))
    B = pd.DataFrame(rows); B.to_csv(OUT / 'bins.csv', index=False)
    R = B.groupby('element').agg(n_bins=('feh_mid', 'size'), spread_dex=('spread_dex', 'median'),
                                 pooled_scatter=('pooled_scatter', 'median'), sep_scatter=('sep_scatter', 'median'),
                                 chi2_nu=('chi2_nu', 'median')).sort_values('sep_scatter', ascending=False)
    n_meas = {e: int(S[e].notna().sum()) for e in elems}
    R['n_measurements'] = R.index.map(n_meas)
    R.round(3).to_csv(OUT / 'ranking.csv')
    print(R.round(2).to_string())

    top = list(R.index[:6]) + list(R.index[-2:])
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharex=True, constrained_layout=True)
    for ax, e in zip(axes.flat, top):
        for g in [g for g in GALAXIES if g not in EXCLUDE]:
            x = S[(S.galaxy == g)].dropna(subset=[e])
            ax.scatter(x.feh, x[e], s=3, c=COLORS[g], alpha=.12, linewidths=0, rasterized=True)
            m = B[(B.element == e) & B.galaxies.str.contains(g)]
            med = x.groupby('bin')[e].agg(['median', 'size'])
            med = med[med['size'] >= NMIN]
            ax.plot(FEH_LO + (med.index + .5) * BIN, med['median'], color=COLORS[g], lw=2.2, label=g)
        r = R.loc[e]
        ax.set_title(f'[{e.title()}/Fe]: {r.sep_scatter:.2f} scatters, {r.spread_dex:.2f} dex', fontsize=10)
        lo, hi = np.nanpercentile(S[e], [1, 99]); ax.set_ylim(lo - .1, hi + .1)
        ax.set(xlabel='[Fe/H]', ylabel=f'[{e.title()}/Fe]'); ax.grid(alpha=.15)
    axes.flat[0].legend(fontsize=8, frameon=False)
    fig.suptitle('H21 dwarfs: median tracks (0.1 dex bins, N>=15); top 6 and bottom 2 by track separation', fontsize=11)
    if not EXCLUDE: fig.savefig(ROOT / 'plots/element_track_separation_top.png', dpi=130)
    plt.close(fig)


if __name__ == '__main__':
    main()
