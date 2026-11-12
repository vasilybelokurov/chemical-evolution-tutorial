# Data

`MANIFEST.json` lists every data file with its size, SHA-256 and origin. The files are attached to the GitHub release
`data-v0.1`; the notebooks download them on first use and check the checksum.

- `{lmc,smc,sgr,gse,fornax}_h21.ecsv`: the members of Hasselquist et al. (2021, Table 2) joined to APOGEE DR17
  `allStarLite-dr17-synspec_rev1`, one observation per star. Columns: Fe, Mg, Si and Ce abundances with errors and
  flags; Teff, logg, S/N, ASPCAPFLAG, EXTRATARG.
- `hasselquist2021_table2.ecsv`: the published membership table, unmodified.
- `h21_membership_summary.csv`: member counts against the paper's Table 1.
- `h21_ssp_bundle.npz`: the yield tables of H21's models (flexCE 2450aa29: Limongi & Chieffi CCSNe, Karakas 2010
  AGB, Iwamoto et al. W70 SNe Ia), per unit mass formed, against metallicity.
- `h21_grid_vice.npz`, `h21_grid_flexce.npz`: H21's 3403-model grid run with VICE 1.3.1 and with flexCE (float32
  tracks); `h21_grid_scores.csv`: their scores against the samples.

`membership_audit.json` is listed in the manifest but not released.
